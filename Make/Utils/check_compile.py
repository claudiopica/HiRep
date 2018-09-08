#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import argparse
import os.path
import sys

no_compile = 'NOCOMPILE='
def search_file(fname, opts):
    for line in open(fname):
        if no_compile in line:
            line_opts = set(line.split(no_compile)[1].strip().split(' '))
            if line_opts.issubset(opts):
                return True
    return False



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, required=True,
                    help='test file')
    parser.add_argument('-n', type=str, required=True,
                    help='group size')
    parser.add_argument('-g', type=str, required=True,
                    help='gauge group')
    parser.add_argument('-r', type=str, required=True,
                    help='group representation')
    parser.add_argument('-s', nargs='+', required=True,
                    help='string of macros for the comparison')


    args = parser.parse_args()
    file=args.f+".c"

    #sys.stderr.write(file+"AAA"+str(args.s[0])+"SSS\n")
    
    macros=[]
    for item in args.s:
        macros.extend(item.replace("-D","").strip().split(' '))
    macros.append(args.g)
    macros.append(args.r)
    macros.append("NG="+args.n)

    #sys.stderr.write(str(macros)+"\n")

    if(os.path.isfile(file)):
        if search_file(file, set(macros)):
            exit(0)
        print args.f
    else:
        sys.stderr.write("Missing file "+file)