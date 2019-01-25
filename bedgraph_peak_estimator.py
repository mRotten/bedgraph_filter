# bedgraph_peak_estimator.py

import argparse

global curr_cluster
curr_cluster = {'members' : [],
                'cl_chr'  : None,
                'cl_start': 0,
                'cl_end'  : 0,
                'maxval'  : 0,
                'minval'  : 0}
global curr_int
curr_int = None

def set_curr_int(line):

    global curr_int

    ll = line.strip('\n').split('\t')

    curr_int = {'chrom': ll[0],
                'start': int(ll[1]),
                'end'  : int(ll[2]),
                'val'  : float(ll[3])}


def add_ci_to_cluster():
    
    global curr_cluster

    val = curr_int['val']
    curr_cluster['members'].append(curr_int)
    curr_cluster['cl_end'] = curr_int['end']
    curr_cluster['maxval'] = max(curr_cluster['maxval'], val)
    curr_cluster['minval'] = min(curr_cluster['minval'], val)


def is_cluster_member(line, min_gap):

    set_curr_int(line)
    
    if curr_int['chrom'] == curr_cluster['cl_chr'] and \
       curr_int['start'] <= curr_cluster['cl_end'] + min_gap:
        return True
    
    else:
        return False


def new_cluster():

    global curr_cluster
    global curr_int
    
    ci = curr_int
    curr_cluster = {'members'  : [ci, ],
                    'cl_chr'   : ci['chrom'],
                    'cl_ start': ci['start'],
                    'cl_end'   : ci['end'],
                    'maxval'   : ci['val'],
                    'minval'   : ci['val']}


def process_cluster(cutoff_frac, value_cutoff):

    global curr_cluster
    cc = curr_cluster

    if cutoff_frac < 1 and cutoff_frac > 0:
        fco = cc['maxval'] * cutoff_frac
    else:
        fco = cc['maxval']

    vco = value_cutoff

    # Filter for cutoff fraction.
    result1 = [x for x in cc['members'] if x['val'] > vco]

    # Filter on cutoff value.
    result2 = [x for x in result1 if x['val'] > fco]

    for r in result2:
        print("\t".join([r['chrom'],
                         str(r['start']),
                         str(r['end']),
                         str(r['val'])]))


def process_bedgraph(args):

    with open(args.path,'r') as infile:

        first_int = infile.readline()

        if "type=" in first_int:
            set_curr_int(infile.readline())
        else:
            set_curr_int(first_int)
        
        for line in infile:
            if is_cluster_member(line, args.mg):
                add_ci_to_cluster()
            else:
                process_cluster(args.fco, args.vco)
                new_cluster()


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('path', help="Path to bedgraph file.",
                        type=str)
    parser.add_argument('-mg', help="Minimum gap allowed in a "\
                        "cluster. Default is 10.",
                        type=int, default=10)
    parser.add_argument('-fco', help="Cutoff fraction. This "\
                        "determines how much of each interval is"\
                        "removed, given as a % of the local maximum. "\
                        "Default is 0.25 (25%).",
                        type=float, default=0.25)
    parser.add_argument('-vco', help="Value cutoff. Minimum value of "\
                        "interval to be put into the destination file. "\
                        "Default is 0.",
                        type=float, default=0)
    args = parser.parse_args()

    process_bedgraph(args)
    

if __name__ == "__main__":
    main()
