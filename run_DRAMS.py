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
