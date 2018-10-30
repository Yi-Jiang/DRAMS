import numpy as np

## To extract highly related sample pairs from GCTA results.
## Yi Jiang, October 2018

## Import modules
import pandas as pd
from struct import unpack, calcsize
import matplotlib.pyplot as plt
import sys, getopt

## Help page
def usage():
    print("")
    print("Usage: python %s --option=<argument>" %sys.argv[0])
    print("  --input_prefix=<STRING>    prefix of input file")
    print("  --output_prefix=<STRING>   prefix of output file")
    print("  --threshold=<FLOAT>        Threshold of genetic relatedness score to extract")
    print("                             highly related pairs [0.65]")
    #print("  --min_loci=<INT>           Minimum number of overlapped loci required to")
    #print("                             estimate sample relatedness scores [400]")
    print("  -h/--help                  Show this information")

## default options
min_loci = 400
threshold = 0.65

## deal with the options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "input_prefix=", "output_prefix=", "threshold=", "min_loci=", ] )
except getopt.GetoptError:
    print("get option error!")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ( "-h", "--help" ):
        usage()
        sys.exit(1)
    else:
        if opt in ( "--input_prefix", ):
            input_prefix = val
        if opt in ( "--output_prefix", ):
            output_prefix = val
        if opt in ( "--threshold", ):
            threshold = float(val)
        if opt in ( "--min_loci", ):
            min_loci = int(val)

## required options
try: input_prefix, output_prefix
except:
    print("\nMissing options!")
    usage()
    sys.exit(2)

def sum_n_vec(n):
    """
    return a vector giving the one less than the sum of the first i
    integers for i from 1 to n. Values are one less than the sum so
    that they can be used to index the elements of a vector storing
    the lower diagonal entries of a symmetric matrix that correspond
    to the diagonal entries
    
    Credit: Davis McCarthy (davismcc)
    https://github.com/davismcc
    Taken from:
    https://github.com/davismcc/GCTAtools/blob/master/GCTAtools.py
    """
    out = [int(0)] * n
    for i in range(n):
        out[i] = int(((i + 1) * (i + 2) / 2) - 1)
    return(out)

def ReadGRMBinN(prefix):
    """
    read the number of SNPs used to compute a GRM from GCTA-format binary file
    adapted from an R function on GCTA website
    
    Credit: Davis McCarthy (davismcc)
    https://github.com/davismcc
    Taken from:
    https://github.com/davismcc/GCTAtools/blob/master/GCTAtools.py
    """
    NFileName = prefix + ".grm.N.bin"
    entry_format = 'f' #N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    ## Open the binary file
    with open(NFileName, mode='rb') as f:
        #entry_count = os.fstat(f.fileno()).st_size / entry_size
        record = f.read(entry_size)
        N = unpack(entry_format, record)[0]
        N = int(N)
    return(N)

def ReadGRMBin(prefix, AllN = False):
    """
    read a GCTA binary GRM file storing relatedness values between individuals
    adapted from an R function on the GCTA website
    
    Credit: Davis McCarthy (davismcc)
    https://github.com/davismcc
    Taken from:
    https://github.com/davismcc/GCTAtools/blob/master/GCTAtools.py
    """
    BinFileName  = prefix + ".grm.bin"
    NFileName = prefix + ".grm.N.bin"
    IDFileName = prefix + ".grm.id"
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
    entry_format = 'f' # N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    ## Read IDs
    ids = pd.read_csv(IDFileName, sep = '\t', header = None)
    ids_vec = ids.iloc[:,1]
    n = len(ids.index)
    ids_diag = ['NA' for x in range(n)]
    n_off = int((n * (n + 1) / 2) - n)
    ids_off = ['NA' for x in range(n_off)]
    ## Generate ids for relatedness values by concatenating individual IDs
    ticker = 0
    for i in range(n):
        for j in range(i):
            if i == j:
                ids_diag[i] = str(ids_vec[i])
            else:
                ids_off[ticker] = str(ids_vec[i]) + '::' + str(ids_vec[j])
                ticker += 1
    ## Read relatedness values
    grm = np.fromfile(BinFileName, dtype = dt)
    ## Read number of markers values
    if AllN:
        N = np.fromfile(NFileName, dtype = dt)
    else:
        with open(NFileName, mode='rb') as f:
            record = f.read(entry_size)
            N = unpack(entry_format, record)[0]
            N = int(N)
    i = sum_n_vec(n)
    out = {'diag': grm[i], 'off': np.delete(grm, i), 'id': ids, 'id_off': ids_off, 'id_diag': ids_diag, 'N': N}
    return(out)

def plot_cor_distribution(d, output_prefix=output_prefix):
    """ 
    Create an histogram for distribution of sample relatedness scores
    """
    for k in d
        plt.hist(list(d[k][:,2]), bins='auto')
        plt.yscale('log')
        plt.title("Histogram of genetic relatedness scores")
        plt.savefig("%s.%s.relatedness.hist.pdf"%(output_prefix,k))

def extract_relate_pairs(d, t=threshold):
    """
    Extract highly related sample pairs based on the distribution of sample relatedness
    relatePairs = { Omics1::Omics2 -> [ [SampleID1, SampleID2, Relatedness], ... ]}
    """
    relatePairs = dict()
    for p in d:
        if d[p]>t:
            s = p.split("::")
            s1 = s[0].split("|")
            s2 = s[1].split("|")
            if len(s1)!=2:
                print("Warning: Sample ID \"%s\" is not formated like \"OmicsType|SampleID\". This sample will be ignored."%s[0])
                continue
            if len(s2)!=2:
                print("Warning: Sample ID \"%s\" is not formated like \"OmicsType|SampleID\". This sample will be ignored."%s[1])
                continue
            k = "%s::%s"%(s1[0],s2[0])
            if k not in relatePairs:
                relatePairs[k] = []
            relatePairs[k].append([ s1[1], s2[1], d[p] ])
    for k in relatePairs:
        relatePairs[k] = np.array(relatePairs[k])
    return relatePairs

def output_pairs(d, output_prefix=output_prefix):
    """
    Write highly related pairs and summary into files
    """
    f1 = open("%s.highlyrelatedpairs.txt"%output_prefix, "w")
    f2 = open("%s.highlyrelatedpairs.summary.txt"%output_prefix, "w")
    f1.write("Omics1\tSampleID1\tOmics2\tSampleID2\tRelatedness\tMatch\n")
    f2.write("Omics1\tOmics2\t#matchedPair\t#mismatchedPair\n")
    for k in d:
        nMatch = 0
        nMismatch = 0
        omics1,omics2 = k.split("::")
        for p in d[k]:
            if p[0]==p[1]:
                m = "Y"
                nMatch += 1
            else:
                m="N"
                nMismatch += 1
            f1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(omics1,p[0],omics2,p[1],p[2],m))
        f2.write("%s\t%s\t%s\t%s\n"%(omics1,omics2,nMatch,nMismatch))
    f1.close()
    f2.close()

if __name__ == "__main__":
    grmres = ReadGRMBin(input_prefix)
    grmres['diag'] = np.ndarray.tolist(grmres['diag'])
    grmres['off'] = np.ndarray.tolist(grmres['off'])
    relatedness = dict(zip(grmres['id_off'], grmres['off']))
    
    # Extract highly related sample pairs
    relatePairs = extract_relate_pairs(relatedness, threshold)
    
    # Plot distribution of sample relatedness scores
    plot_cor_distribution(relatePairs)
    
    # Write highly related pairs and summary into files
    output_pairs(relatePairs)
