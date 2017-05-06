#!/usr/bin/evn python

import sys, argparse

get_range = lambda x: [min(x), max(x)]
get_size = lambda x: max(x)-min(x)

def build_dict(infile):

    coords = {}

    for line in open( infile ):
        row = []
        fields = line.split()
        ref_st, ref_end, qry_st, qry_end = map(int, fields[0:4])
        qry_chr = fields[14]
        row += [ref_st, ref_end] + sorted([qry_st, qry_end])

        if qry_chr not in coords: coords[qry_chr] = []
        coords[qry_chr].append( row )

    return coords

def get_parser():
    """ Parse input """
    desc = "Tool for a crude clustering of homologous \
             regions based on mummer coords file"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('coords', type=argparse.FileType("rU"), \
                        help="mummer show-coords output")
    parser.add_argument('ref', type=str, \
                        help="reference chromosome")
    parser.add_argument('output', nargs="?", type=argparse.FileType('wt'), default=sys.stdout, \
                        help="output file name, default stdout")
    parser.add_argument('-d', '--dist', type=int, default=1000, \
                        help="maximal distance for clustering, default 1000")
    return parser

def cluster_hom_regions( coords, dist ):

    while True:
        prev = len( coords )
        newcoords = clustering( coords, dist )
        curr = len( newcoords )
        if curr == prev or curr == 1: break
        coords = newcoords
    return newcoords

def write( rows, rchr, qchr, outf ):

    for row in rows:
        s1, s2 = get_ind_size( row )
        if s1 <= 1000 and s2 <= 1000: continue
        print("%s\t%s\t%s" % (rchr, qchr, "\t".join(map(str,row))), \
            file=outf)

def main():
    parser = get_parser()
    args = parser.parse_args()
    infile = args.coords
    rchr = args.ref
    dist = args.dist
    outf = args.output
    print("Homologous regions will be clustered if they are within %d bp." \
            % (5*dist), file=sys.stderr)

    coords = build_dict(infile)

    qry_chrs = list( coords.keys() )
    for qchr in qry_chrs:
        results = cluster_hom_regions( coords[qchr], dist )
        write( results, rchr, qchr, outf )

def get_ind_size( l ):

    s1, e1, s2, e2 = l
    return e1-s1, e2-s2
    
def gap_size( l ):

    st1, end1, st2, end2 = l
    sum_size = get_size(l)
    size1, size2 = get_ind_size( l )

    # handle overlay case
    if size1 + size2 > sum_size:
        return 0.5
    else:
        return st2 - end1 if st2 >= end1 else st1 - end2

def is_inclusive( l ):

    sum_size = get_size(l)
    size1, size2 = get_ind_size( l )

    if sum_size == size1: return 1
    elif sum_size == size2: return 2
    else: return 0

def handle_overlay( rl, ql, dist ):
    wh = is_inclusive( rl )
    if wh == 1:
        return get_range(rl), ql[0:2]
    elif wh == 2:
        return get_range(rl), ql[2:4]
    else:
        return get_range(rl), get_range(ql)

def clustering( coords, dist ):

    # initialize
    results = []
    prevr, prevq = coords[0][0:2], coords[0][2:4]

    # iterate over sequential pairs
    for i in range(1, len(coords) ):
        currr, currq = coords[i][0:2], coords[i][2:4]
        #print(prevr, prevq, currr, currq, file=sys.stderr)

        rl = prevr + currr
        ql = prevq + currq
        rgap = gap_size( rl )
        qgap = gap_size( ql )
        #print("gaps:", rgap, qgap, file=sys.stderr)

        if rgap == 0.5 and qgap <= 5*dist:
            prevr, prevq = handle_overlay( rl, ql, dist )
        elif qgap == 0.5 and rgap <= 5*dist:
            prevq, prevr = handle_overlay( ql, rl, dist )
        elif rgap + qgap <= 5*dist:
            prevr, prevq = get_range( rl ), get_range( ql )
        else:
            results.append( prevr + prevq )
            prevr, prevq = currr, currq

    results.append( prevr + prevq )
    return results

if __name__=="__main__":
    main()
