"""This script takes an output file from the crosspoints script (X.crosspoints.csv) and a minimum bin size and identifies a list of bins and their genotypes for each RIL."""
from collections import OrderedDict
import math
import os

def _bins_batcher(input_path, output_path, min_bin_size, binmap_id=False):
    input_list = input_path
    if len(input_list) == 1:
        bins(input_path[0], output_path, min_bin_size, binmap_id)
    else:
        output_folder = output_path
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for file in input_list:
            print("Determining bins for:", file)
            out_file = os.path.join(output_folder, os.path.basename(file) + ".bins.csv")
            try:
                bins(file, out_file, min_bin_size, os.path.basename(file))
            except Exception as detail:
                print("FAILED! Details:")
                print(detail)

def bins(input_path, output_path, min_bin_size, binmap_id=False):
    line_cps = OrderedDict()
    chrom_len = 0
    with open(input_path) as input_file:
        for line in input_file:
            name, cps = line.split(",", 1)
            cps = [item.strip() for item in cps.split(",") if item.strip() != ""]
            for i in range(0, len(cps), 2):
                cps[i] = int(cps[i])
            line_cps[name] = cps
            if cps[-1] > chrom_len:
                chrom_len = cps[-1]

    with open(output_path, "w") as clear_output:
        pass

    cp_loc_count = {chrom_len: float('inf'), 0: float('inf')}
    for cps in line_cps.values():
        for i in range(0, len(cps), 2):
            cp_loc_count[cps[i]] = cp_loc_count.get(cps[i], 0) + 1
    all_cp = sorted(cp_loc_count.keys())
    print(all_cp)

    crosspoint_groups = [[all_cp[0]]]
    for cp in all_cp[1:]:
        if cp - crosspoint_groups[-1][-1] < min_bin_size:
            crosspoint_groups[-1].append(cp)
        else:
            crosspoint_groups.append([cp])
    print(crosspoint_groups)

    bin_bounds = []
    crosspoint_groups.reverse()

    for group in crosspoint_groups:
        group_len = group[-1] - group[0]
        expanded_group = []
        for cp in group:
            count = cp_loc_count[cp] if cp_loc_count[cp] != float('inf') else 1
            expanded_group.extend([cp] * count)

        print("\n\n" + "=" * 80)
        print(f"{len(expanded_group)} crosspoints, start={group[0]}, end={group[-1]}")
        print("." * 80)

        if group_len == 0:
            print("\nN=1")
            bin_bounds.append(group[0])
            print("F", bin_bounds[-1])
            continue

        max_new_cp = group_len // min_bin_size + 1

        if max_new_cp < 2:
            print("\nN=1")
            print("g", _bin_bound_visualize(expanded_group, group[0], group[-1]))
            bin_bounds.append(_crosspoint_avg(cp_loc_count, group, chrom_len))
            print("f", _bin_bound_visualize(bin_bounds[-1:], group[0], group[-1], aura=min_bin_size))
            continue

        solution_list = []
        for cp_count in range(max_new_cp, 0, -1):
            print(f"\nN={cp_count}\n")
            print("g", _bin_bound_visualize(expanded_group, group[0], group[-1]))

            start_cp_dist = group_len / float(cp_count)
            if min_bin_size > start_cp_dist:
                start_cp_dist = min_bin_size
            start_offset = group_len / 2.0 - (cp_count * start_cp_dist) / 2.0 + start_cp_dist / 2
            km_points = [int(j * start_cp_dist + group[0] + start_offset) for j in range(cp_count)]
            km_groups = [[] for _ in range(len(km_points))]

            nearest = 0
            for cp in group:
                while nearest + 1 < len(km_points) and abs(cp - km_points[nearest + 1]) < abs(cp - km_points[nearest]):
                    nearest += 1
                km_groups[nearest].append(cp)
            print("u", _bin_bound_visualize(km_points, group[0], group[-1], aura=min_bin_size))

            memo = set()
            adjustment_needed = True
            while adjustment_needed:
                km_points = [_crosspoint_avg(cp_loc_count, k, chrom_len) for k in km_groups]

                for k in range(1, len(km_points) - 1):
                    if math.isnan(km_points[k]):
                        next_not_nan = k + 1
                        while next_not_nan < len(km_points) and math.isnan(km_points[next_not_nan]):
                            next_not_nan += 1
                        if next_not_nan < len(km_points):
                            km_points[k] = (km_points[k - 1] + km_points[next_not_nan]) // 2

                get_overlaps = lambda: [
                    (min_bin_size - abs(km_points[k] - km_points[k - 1]), k)
                    for k in range(1, len(km_points))
                    if (min_bin_size - abs(km_points[k] - km_points[k - 1])) > 0
                ]
                overlap_adjustment_performed = False
                id_tuple = (tuple(len(g) for g in km_groups), tuple(km_points))

                if id_tuple in memo:
                    overlaps = get_overlaps()
                    while overlaps:
                        for overlap, k in overlaps:
                            left_shift = overlap - (overlap // 2)
                            right_shift = overlap - left_shift
                            if k == 1 and left_shift > km_points[k - 1] - group[0]:
                                left_shift = km_points[k - 1] - group[0]
                                right_shift = overlap - left_shift
                            elif k == len(km_points) - 1 and right_shift > group[-1] - km_points[k]:
                                right_shift = group[-1] - km_points[k]
                                left_shift = overlap - right_shift
                            km_points[k] += right_shift
                            km_points[k - 1] -= left_shift
                        overlaps = get_overlaps()
                    id_tuple = (tuple(len(g) for g in km_groups), tuple(km_points))

                for k in range(len(km_groups)):
                    if k > 0:
                        while km_groups[k] and abs(km_groups[k][0] - km_points[k - 1]) < abs(km_groups[k][0] - km_points[k]):
                            km_groups[k - 1].append(km_groups[k].pop(0))
                    if k < len(km_groups) - 1:
                        while km_groups[k] and abs(km_groups[k][-1] - km_points[k + 1]) < abs(km_groups[k][-1] - km_points[k]):
                            km_groups[k + 1].insert(0, km_groups[k].pop())

                if id_tuple in memo:
                    adjustment_needed = False
                    print("Done.")
                else:
                    memo.add(id_tuple)
                    print("adj", _bin_bound_visualize(km_points, group[0], group[-1], aura=min_bin_size))

            print("f", _bin_bound_visualize(km_points, group[0], group[-1], aura=min_bin_size))
            print("g", _bin_bound_visualize(expanded_group, group[0], group[-1]))

            dists = [[]]
            nearest = 0
            for cp in group:
                while nearest + 1 < len(km_points) and abs(cp - km_points[nearest + 1]) < abs(cp - km_points[nearest]):
                    nearest += 1
                    dists.append([])
                dists[-1].append(abs(cp - km_points[nearest]))
            average_average_group_dist = sum(sum(ds) / float(len(ds)) for ds in dists if ds)
            print("S", average_average_group_dist)

            solution_list.append((average_average_group_dist, km_points))

        _, best_bounds = sorted(solution_list)[0]
        bin_bounds += best_bounds

    bin_bounds.sort()
    print("\nBin Bounds:", bin_bounds)

    bin_genotypes = OrderedDict()
    for line in line_cps:
        bin_genotypes[line] = []
        line_cp_locs = list(enumerate((cp for i, cp in enumerate(line_cps[line]) if i % 2 == 0)))
        scan_lower_i, scan_lower_bound = line_cp_locs[0]
        scan_upper_i, scan_upper_bound = line_cp_locs[0]
        bin_iter = enumerate(bin_bounds)
        next(bin_iter)
        for bin_i, bin_upper_bound in bin_iter:
            bin_lower_bound = bin_bounds[bin_i - 1]
            while scan_lower_i + 1 < len(line_cp_locs) and line_cp_locs[scan_lower_i + 1][1] <= bin_lower_bound:
                scan_lower_i, scan_lower_bound = line_cp_locs[scan_lower_i + 1]
            while scan_upper_bound < bin_upper_bound and scan_upper_i + 1 < len(line_cp_locs):
                scan_upper_i, scan_upper_bound = line_cp_locs[scan_upper_i + 1]
            scan_range_contents = line_cps[line][scan_lower_i * 2:scan_upper_i * 2 + 1]
            bin_contents = [bin_lower_bound] + scan_range_contents[1:-1] + [bin_upper_bound]
            genotype_weights = {}
            for seg_i in range(1, len(bin_contents), 2):
                g = bin_contents[seg_i]
                genotype_weights[g] = genotype_weights.get(g, 0) + bin_contents[seg_i + 1] - bin_contents[seg_i - 1]
            bin_genotypes[line].append(max(genotype_weights.items(), key=lambda x: x[1])[0])

    bin_centers = [(bin_bounds[i - 1] + bin_bounds[i]) // 2 for i in range(1, len(bin_bounds))]

    for i in range(len(bin_centers) - 1, 0, -1):
        if all(bin_genotypes[line][i - 1] == bin_genotypes[line][i] for line in bin_genotypes):
            print("Combined Identical Bins.", bin_bounds[i - 1], ">", bin_bounds[i], "<", bin_bounds[i + 1])
            del bin_centers[i]
            del bin_bounds[i]
            bin_centers[i - 1] = (bin_bounds[i - 1] + bin_bounds[i]) // 2
            for line in bin_genotypes:
                del bin_genotypes[line][i]

    with open(output_path, "w", encoding="utf-8") as outfile:
        if binmap_id:
            outfile.write("##binmap id," + ",".join([binmap_id] * len(bin_centers)) + "\n")
        outfile.write("##bin start," + ",".join([str(i) for i in bin_bounds[:-1]]) + "\n")
        outfile.write("##bin end," + ",".join([str(i) for i in bin_bounds[1:]]) + "\n")
        outfile.write("bin center," + ",".join([str(i) for i in bin_centers]) + "\n")
        for line in bin_genotypes:
            outfile.write("%s,%s\n" % (line, ",".join(str(x) for x in bin_genotypes[line])))

def _crosspoint_avg(weights, cp_list, chrom_len):
    if len(cp_list) < 1:
        return float('nan')
    total_weight = sum(weights[k] for k in cp_list)
    if total_weight == 0:
        return 0
    avg = sum(k * weights[k] for k in cp_list) / float(total_weight)
    if math.isnan(avg) and (0 in cp_list or chrom_len in cp_list):
        avg = 0 if 0 in cp_list else chrom_len
    return int(avg)

def _bin_bound_visualize(cps, begin, end, bins=75, aura=0):
    empty_sym, aura_sym, aura_end_sym, aura_overlap_sym = (" ", "-", "|", "x")
    bin_size = (end - begin) / float(bins - 1)
    binned = [0] * bins
    aura_size = int(aura // bin_size // 2 if aura != 0 else 0)
    for cp in cps:
        bin_loc = int((cp - begin) // bin_size)
        if 0 <= bin_loc < len(binned):
            binned[bin_loc] += 1
    loc_map = [empty_sym] * bins
    for i in range(len(binned)):
        if binned[i] > 0:
            loc_map[i] = str(binned[i]) if binned[i] <= 10 else "G"
            for j in range(i - aura_size, i + aura_size + 1):
                if 0 <= j < len(loc_map):
                    if j == i - aura_size or j == i + aura_size:
                        loc_map[j] = aura_end_sym if loc_map[j] == empty_sym else aura_overlap_sym
                    elif j != i:
                        loc_map[j] = aura_sym if loc_map[j] == empty_sym else aura_overlap_sym
    return "".join(loc_map)

_cl_entry = _bins_batcher

def _parser(parser_add_func, name):
    p = parser_add_func(name, description=__doc__)
    p.add_argument("-i", "--input", metavar="PATH", dest='input_path', required=True, nargs='*',
                   help="Path to a crosspoints CSV or multiple.")
    p.add_argument("-o", "--output", metavar="PATH", dest='output_path', required=True,
                   help="Path for output.")
    p.add_argument("-l", "--min-bin-size", metavar="INT", dest='min_bin_size', required=True, type=int,
                   help="Minimum size of a bin.")
    p.add_argument("-n", "--binmap-id", metavar="ID", dest="binmap_id", default=False, type=str,
                   help="Optional binmap ID.")
    return p
