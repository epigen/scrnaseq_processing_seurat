#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import argparse
import pysam

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="input fasta file")
    parser.add_argument(
        "--outputs", nargs="+", help="list of output files", required=True
    )
    args = parser.parse_args()
    NIDS = len(args.outputs)

    fasta = pysam.FastaFile(args.infile)

    outs = [open(f, "w+") for f in args.outputs]
    outidx = 0
    for name in fasta.references:
        seq = fasta.fetch(name)
        outs[outidx].write(">{}\n{}\n".format(name, seq))
        outidx += 1
        if outidx == NIDS:
            outidx = 0

    for out in outs:
        out.close()
