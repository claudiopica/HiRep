#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os.path
import sys
from typing import Any, Set

no_compile = "NOCOMPILE="


def search_file(fname: str, opts: Set[Any]) -> bool:
    logic_string_file = []
    for line in open(fname):
        if no_compile in line:
            logic_string_line = []
            line_opts = set(line.split(no_compile)[1].strip().split(" "))
            for element in line_opts:
                if str(element).find("!") > -1:
                    new_element = set(str(element).split("!")[1:])
                    if new_element.issubset(opts):
                        logic_string_line.append(False)
                    else:
                        logic_string_line.append(True)
                else:
                    if set([element]).issubset(opts):
                        logic_string_line.append(True)
                    else:
                        logic_string_line.append(False)
            logic_string_file.append(all(logic_string_line))
    return any(logic_string_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", type=str, required=True, help="test file")
    parser.add_argument("-n", type=str, required=True, help="group size")
    parser.add_argument("-g", type=str, required=True, help="gauge group")
    parser.add_argument("-r", type=str, required=True, help="group representation")
    parser.add_argument(
        "-s", nargs="+", required=True, help="string of macros for the comparison"
    )

    args = parser.parse_args()
    file = f"{args.f}.c"

    macros = []
    for item in args.s:
        macros.extend(item.replace("-D", "").strip().split(" "))
    macros.append(args.g)
    macros.append(args.r)
    macros.append(f"NG={args.n}")

    if os.path.isfile(file):
        if search_file(file, set(macros)):
            exit(0)
        print(args.f)
    else:
        sys.stderr.write(f"Missing file {file}")
