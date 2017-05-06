#!/usr/bin/evn python

import sys, argparse
from Bio import SeqIO

def get_parser():
    desc = "extract haplotigs and primary contigs"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('contig', type=str, help="contig name")
    parser.add_argument('inputfile', type=argparse.FileType('rU'), \
        default=sys.stdin, help="genome sequences")
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    infile = args.inputfile
    contig = args.contig

    _process( infile, contig )

def _process( infile, contig ):
    inhandle = open( infile, "rU" )
    fnm = contig+"_p.fa"
    pout = open( fnm, "w" )
    fnm = contig+"_h.fa"
    hout = open( fnm, "w" )

    for record in SeqIO.parse( inhandle, "fasta" ):
        if record.name.startswith(contig):
            if record.name.find("_") > 0: SeqIO.write(record,hout,"fasta")
            else: SeqIO.write(record, pout, "fasta")

    inhandle.close()
    pout.close()
    hout.close()

main()
