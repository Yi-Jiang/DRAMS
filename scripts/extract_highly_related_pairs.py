import numpy as np
import pandas as pd
from struct import unpack, calcsize
import sys, re, getopt, gzip

def usage():
    print("")
    print("Usage: python %s --option=<argument>" %sys.argv[0])
    print("  --input_prefix=<STRING>    prefix of input file")
    print("  --output_prefix=<STRING>   prefix of output file")
    print("  --threshold=<FLOAT>        Threshold of genetic relatedness score to extract")
    print("                             highly related pairs. The threshold will be set")
    print("                             automatically if not set")
    print("  --min_loci=<INT>           Minimum number of overlapped loci required to")
    print("                             estimate sample relatedness scores [400]")
    print("  -h/--help                  Show this information")

try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "input_prefix=", "output_prefix=", "threshold=", "min_loci=", ] )
except getopt.GetoptError:
    print("get option error!")
    usage()
    sys.exit(2)

## default options
min_loci = 400

## deal with the options
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
                ids_off[ticker] = str(ids_vec[i]) + '_' + str(ids_vec[j])
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

def plot_correlogram(df,figsize=(20,20)):
    """ 
    Create an n x n matrix of scatter plots for every
    combination of numeric columns in a pandas dataframe

    Credit: Corey Chivers (cjbayesian)
    https://github.com/cjbayesian
    Taken from:
    https://gist.github.com/cjbayesian/f0f127cc57f26d968f10
    """
    cols = list(df.columns[df.dtypes=='float64'])
    n = len(cols)
    fig, ax = plt.subplots(n,n,figsize=figsize)
    for i,y in enumerate(cols):
        for j,x in enumerate(cols):
            if i != n-1:
                ax[i,j].xaxis.set_ticklabels([])
            if j != 0:
                    ax[i,j].yaxis.set_ticklabels([])
            if i != j:
                try:
                    tmp = df[[x,y]].copy()
                    tmp.dropna(inplace=True)
                    ax[i,j].plot(tmp[x].values,tmp[y].values,'.',markersize=0.5,alpha=0.5,color='black')
                except:
                    pass
            else:
                midx = df[x].min() + (df[x].max() - df[x].min())/2.0
                midy = df[y].min() + (df[y].max() - df[y].min())/2.0
                ax[i,j].text(midx, midy, y.replace(' ','\n'),
                             horizontalalignment='center',
                             verticalalignment='center')
                ax[i,j].set_ylim((df[y].min(),df[y].max()))
                ax[i,j].set_xlim((df[x].min(),df[x].max()))

def plot_cor_distribution(df,figsize=(20,20)):
    """ 
    Create an histogram for distribution of sample relatedness scores
    """
    

def find_threshold(relatedness,minloci=min_loci):
    """
    Extract highly related sample pairs based on the distribution of sample relatedness
    """
    

def extract_relate_pairs(relatedness,minloci=min_loci):
    """
    Extract highly related sample pairs based on the distribution of sample relatedness
    """
    

def summary_relate_pairs(relatedness):
    """
    Generate a summary table for highly related sample pairs
    """

d = ReadGRMBin("merge.braingvex-1000g")
d.keys()
# dict_keys(['diag', 'off', 'id', 'id_off', 'id_diag', 'N'])

find_threshold(d['off'])




