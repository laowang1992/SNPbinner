"""Sets up the main argparser and imports and sets up scripts in the program."""

import __init__ as snpbinner
import sys
import argparse

def main():
    program_dict = {}
    for module in snpbinner.__all__:
        program_dict[module] = getattr(snpbinner, module)

    main_parser = argparse.ArgumentParser(
        description=(
            "The snpbinner package can be run in three ways. "
            "It can be run as an executable directly from the command line, "
            "run as a python program with `python snpbinner`, or it can be imported "
            "as a python module `snpbinner` and used with other python scripts."
        )
    )

    prog_sub = main_parser.add_subparsers(
        title="Program",
        dest="program",
        help="Specifies which program in the snpbinner package should be run."
    )

    # Python 3: ensure `dest` is used and `required=True` (since Python 3.7)
    prog_sub.required = True

    program_run_dict = {
        name: program_dict[name]._cl_entry
        for name in program_dict
    }

    program_parser_dict = {
        name: program_dict[name]._parser(prog_sub.add_parser, name)
        for name in program_dict
    }

    args = main_parser.parse_args(sys.argv[1:])

    program_to_run = args.program
    arg_dict = vars(args)
    arg_dict_to_pass = {
        key: value for key, value in arg_dict.items()
        if key != "program" and value is not None
    }

    program_run_dict[program_to_run](**arg_dict_to_pass)

if __name__ == '__main__':
    main()
