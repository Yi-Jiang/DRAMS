## To extract highly related sample pairs from GCTA results.
## Yi Jiang, October 2018

## Import modules
import pandas as pd
import numpy as np
from struct import unpack, calcsize
import matplotlib.pyplot as plt
import sys, getopt

## Help page
def usage():
    print("")
    print("Usage: python %s --option=<argument>" %sys.argv[0])
    print("  --input_prefix=<STRING>    Prefix of input file")
    print("  --output_prefix=<STRING>   Prefix of output file")
    print("  --threshold=<FLOAT>        Threshold of genetic relatedness score to extract")
    print("                             highly related pairs [0.65]")
    print("  --min_loci=<INT>           Minimum number of overlapped loci required to")
    print("                             estimate sample relatedness scores [0]")
    print("  --plot                     Plot histograms of sample relatedness scores")
    print("  -h/--help                  Show this information")

## Default options
min_loci = 0
threshold = 0.65
plotHist = 0

## Deal with the options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "input_prefix=", "output_prefix=", "threshold=", "min_loci=", "plot", ] )
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
        if opt in ( "--plot", ):
            plotHist = 1

## Required options
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
    if AllN:
        out = {'diag': grm[i], 'off': np.delete(grm, i), 'id': ids, 'id_off': ids_off, 'id_diag': ids_diag, 'N_diag': N[i], 'N_off': np.delete(N, i)}
    else:
        out = {'diag': grm[i], 'off': np.delete(grm, i), 'id': ids, 'id_off': ids_off, 'id_diag': ids_diag, 'N': N}
    return(out)

def plot_cor_distribution(d, output_prefix=output_prefix):
    """ 
    Create an histogram for distribution of sample relatedness scores
    """
    for k in d:
        plt.hist(list(map(float, d[k][:,2])), bins="auto", color='skyblue', ec='skyblue')
        plt.yscale('log')
        plt.title("%s"%(k.replace("::"," vs. ")))
        plt.ylabel('Frequency')
        plt.xlabel('Genetic relatedness scores')
        plt.savefig("%s.%s.relatedness.hist.pdf"%(output_prefix,k.replace("::","-")))

def format_pairs(d):
    """
    Format sample pairs
    dd = { Omics1::Omics2 -> [ [SampleID1, SampleID2, Relatedness], ... ]}
    """
    dd = dict()
    for p in d:
        s = p.split("::")
        s1 = s[0].split("|")
        s2 = s[1].split("|")
        if len(s1)!=2:
            print("Warning: Sample ID \"%s\" is not formated like \"OmicsType|SampleID\". This sample will be ignored."%s[0])
            continue
        if len(s2)!=2:
            print("Warning: Sample ID \"%s\" is not formated like \"OmicsType|SampleID\". This sample will be ignored."%s[1])
            continue
        if s1[0]==s2[0]:
            continue
        k = "%s::%s"%(s1[0],s2[0])
        if k not in dd:
            dd[k] = []
        dd[k].append([ s1[1], s2[1], d[p] ])
    for k in dd:
        dd[k] = np.array(dd[k])
    return dd

def extract_relate_pairs(d, t, N, minN):
    """
    Extract highly related sample pairs based on the distribution of sample relatedness
    relatePairs = { Omics1::Omics2 -> [ [SampleID1, SampleID2, Relatedness], ... ]}
    """
    relatePairs = dict()
    for k in d:
        if k not in relatePairs:
            relatePairs[k] = []
        for i in range(len(d[k])):
            r = d[k][i]
            n = N[k][i]
            if float(n[2])>minN and float(r[2])>t:
                relatePairs[k].append(r)
    for k in relatePairs:
        relatePairs[k] = np.array(relatePairs[k])
    return relatePairs

def output_pairs(d, output_prefix=output_prefix):
    """
    Write highly related pairs and summary into files
    """
    f1 = open("%s.highlyrelatedpairs.txt"%output_prefix, "w")
    f2 = open("%s.highlyrelatedpairs.summary.txt"%output_prefix, "w")
    f3 = open("%s.highlyrelatedpairs.cytoscape.txt"%output_prefix, "w")
    f1.write("OmicsType1\tSampleID1\tOmicsType2\tSampleID2\tRelatedness\tMatch\n")
    f2.write("OmicsType1\tOmicsType2\t#matchedPair\t#mismatchedPair\n")
    f3.write("Source\tInteraction\tTarget\n")
    for k in d:
        nMatch = 0
        nMismatch = 0
        omics1,omics2 = k.split("::")
        for r in d[k]:
            if r[0]==r[1]:
                m = "Y"
                nMatch += 1
                f3.write("%s|%s\tMatch\t%s|%s\n"%(omics1,r[0],omics2,r[1]))
            else:
                m="N"
                nMismatch += 1
                f3.write("%s|%s\tMismatch\t%s|%s\n"%(omics1,r[0],omics2,r[1]))
            f1.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(omics1,r[0],omics2,r[1],r[2],m))
        f2.write("%s\t%s\t%s\t%s\n"%(omics1,omics2,nMatch,nMismatch))
    f1.close()
    f2.close()
    f3.close()

if __name__ == "__main__":
    # Read GCTA results
    grmres = ReadGRMBin(input_prefix, AllN=True)
    
    # Format results (relatedness scores)
    grmres['diag'] = np.ndarray.tolist(grmres['diag'])
    grmres['off'] = np.ndarray.tolist(grmres['off'])
    relatedness = format_pairs(dict(zip(grmres['id_off'], grmres['off'])))

    # Format results (Number of loci being used)
    grmres['N_diag'] = np.ndarray.tolist(grmres['N_diag'])
    grmres['N_off'] = np.ndarray.tolist(grmres['N_off'])
    Nloci = format_pairs(dict(zip(grmres['id_off'], grmres['N_off'])))
    
    # Plot distribution of sample relatedness scores
    if plotHist==1:
        plot_cor_distribution(relatedness)
    
    # Extract highly related sample pairs
    relatePairs = extract_relate_pairs(relatedness, threshold, Nloci, min_loci)
    
    # Write highly related pairs and summary into files
    output_pairs(relatePairs)
