## To extract highly related data pairs from GCTA results.
## Yi Jiang, October 2018

## Import modules
import sys, re, getopt, gzip, copy, math
from collections import defaultdict
import numpy as np

## Help page
def usage():
    print("")
    print("Usage: python %s -option <argument>" %sys.argv[0])
    print("  --pair=<STRING>     A file of highly related data pairs.")
    print("  --prior=<STRING>    A file of omics priority.")
    print("                        Col1: Omics type. Multiple omics types separated by comma will be considered as one type.")
    print("                        Col2: Omics priority. Recommend to use proportion of sex-matched samples as priority (Integer started from 1).")
    print("  --nsex=<STRING>     A file of sample list with reported sex.")
    print("  --gsex=<STRING>     A file of sample list with genetic inferred sex.")
    print("  --cutoff=<STRING>   Relatedness score cutoff to extract highly related data pairs. [0.65]")
    print("  --coef=<STRING>     A list of coeffecients corresponding to \"intercept,a,b,c\" in Logistic Regression model. ")
    print("                        The values were separated by comma with no space. [0,24.14,0.89,3.39]")
    print("  --train=<STRING>    A list of data pairs with known switch directions being used as training set. If this parameter is specified, it will only train the logistic regression model.")
    print("  --stringent         To realign data ID stringently. In this mode, the data groups with less than three IDs or with no shared ID will be discarded, as it's hard to determine the IDs for these situations..")
    print("  --output=<STRING>   Prefix of output files.")
    print("  -h/--help           Show this information.")

## Default options
cutoff = 0.65
coef = [0,24.14,0.89,3.39]
train = 0
stringent = 0

## Deal with the options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "pair=", "prior=", "nsex=", "gsex=", "cutoff=", "coef=", "train=", "stringent", "output=", ] )
except getopt.GetoptError:
    print("get option error!")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ( "-h", "--help" ):
        usage()
        sys.exit(1)
    else:
        if opt in ( "--pair", ):
            pairFile = val
        if opt in ( "--prior", ):
            priorFile = val
        if opt in ( "--nsex", ):
            nsexFile = val
        if opt in ( "--gsex", ):
            gsexFile = val
        if opt in ( "--cutoff", ):
            cutoff = float(val)
        if opt in ( "--coef", ):
            coef = list(map(float, val.replace(" ","").split(",")))
        if opt in ( "--train", ):
            train = val
        if opt in ( "--stringent", ):
            stringent = 1
        if opt in ( "--output", ):
            output = val

### Required options
#try: pairFile, priorFile, nsexFile, gsexFile, output
#except:
#    print("\nMissing options!")
#    usage()
#    sys.exit(2)

## Initialize variables
omicsSurrogate = dict()  # omicsType -> mergedOmicsTypes(surrogateOmicsType)
omicsPriority = dict()   # surrogateOmicsType -> priorityScore
nsex = dict()      # sampleID -> reportedSex
gsex = dict()      # omicsType -> sampleID -> geneticSex
match = dict()       # omicsType -> sampleID -> [sampleID1, sampleID2, ...]
matchT = dict()      # omicsType -> sampleID -> [omicsType1, omicsType2, ...]
mismatch = dict()    # omicsType -> sampleID -> [sampleID1, sampleID2, ...]
mismatchT = dict()   # omicsType -> sampleID -> [omicsType1, omicsType2, ...]
indegree = dict()    # omicsType|sampleID  -> Number_of_arrow_pointed_to
idsInNet = dict()    # 1 -> {omicsType|sampleID -> sampleID, ... }, 2 -> {}, 3 -> {}, ... #Each subdict includes samples in a network
trueID = dict()        # Group -> [[omicsType, originalID, trueID, switchedOrNot, gSex, nSex_TrueID, SexMatchOrNot], ... ]
nSwitch = dict()   # omicsType -> Number_of_samples_switched_ID
idsWithConnection = []       # [omicsType|sampleID, ... ]
idsNotWithConnection = []    # [omicsType|sampleID, ... ]

def read_priority(infile):
    """
    Read omics priority.
    """
    print("\n## Reading omics types and priority ... ")
    nOmics = 0
    if re.search(r'\.gz$',infile):
        f = gzip.open(infile,'rt')
    else:
        f = open(infile,'r')
    while 1:
        l = f.readline()
        if not l:
            break
        s = l.rstrip("\n").split("\t")
        ss = s[0].split(",")
        for t in ss:
            t = t.strip()
            omicsSurrogate[t] = "|".join(ss)
            omicsPriority[omicsSurrogate[t]] = float(s[1])
            nSwitch[t] = 0
            nOmics += 1
    f.close()
    omicsSurrogateValues = set(omicsSurrogate.values())
    print("  A total of %s omics types were recogized: %s"%(nOmics,", ".join(omicsSurrogate.keys())))
    print("  A total of %s omics types after merged: %s"%(len(omicsPriority), ", ".join(omicsSurrogateValues)))
    print()

def read_reported_sex(infile):
    """
    Read sample list with reported sex information.
    """
    print("## Reading sample list with reported sex information ... ")
    if re.search(r'\.gz$',infile):
        f = gzip.open(infile,'rt')
    else:
        f = open(infile,'r')
    f.readline()
    n = 0
    male = 0
    female = 0
    unknownnsex = 0
    while 1:
        l = f.readline()
        if not l:
            break
        s = l.rstrip("\n").split("\t")
        if s[1]=="1" or s[1]=="M" or re.match("male",s[1],re.I):
            nsex[s[0]] = "1"
            male += 1
        elif s[1]=="2" or s[1]=="F" or re.match("female",s[1],re.I):
            nsex[s[0]] = "2"
            female += 1
        else:
            nsex[s[0]] = "0"
            unknownnsex += 1
        n += 1
    f.close()
    print("  Number of samples recognized: %s"%n)
    print("  Number of male samples: %s"%male)
    print("  Number of female samples: %s"%female)
    print("  Number of samples with unknown reported sex: %s"%unknownnsex)
    print()

def read_genetic_sex(infile):
    """
    Read sample list with genetic inferred sex information.
    """
    print("## Reading sample list with genetic inferred sex information ... ")
    SexMatched = dict()
    SexMismatched = dict()
    unknowngsex = dict()
    if re.search(r'\.gz$',infile):
        f = gzip.open(infile,'rt')
    else:
        f = open(infile,'r')
    f.readline()
    while 1:
        l = f.readline()
        if not l:
            break
        s = l.rstrip("\n").split("\t")
        if s[0] not in gsex:
            gsex[s[0]] = dict()
            SexMatched[s[0]] = 0
            SexMismatched[s[0]] = 0
            unknowngsex[s[0]] = 0
        if s[2]=="1" or s[2]=="M" or re.match("male",s[2],re.I):
            gsex[s[0]][s[1]] = "1"
            if gsex[s[0]][s[1]]==nsex[s[1]]:
                SexMatched[s[0]] += 1
            else:
                SexMismatched[s[0]] += 1
        elif s[2]=="2" or s[2]=="F" or re.match("female",s[2],re.I):
            gsex[s[0]][s[1]] = "2"
            if gsex[s[0]][s[1]]==nsex[s[1]]:
                SexMatched[s[0]] += 1
            else:
                SexMismatched[s[0]] += 1
        else:
            gsex[s[0]][s[1]] = "0"
            unknowngsex[s[0]] += 1
    f.close()
    print("  +" + "-"*91 + "+")
    print("  |" + "%15s"%("Omics Type") + " | " + "%25s"%("SexMatched (nSex==gSex)") + " | " + "%26s"%("SexMismatched (nSex!=gSex)") + " | " + "%15s"%("Unknown gSex") + " |")
    for i in omicsSurrogate:
        print("  |" + "%15s"%(i) + " | " + "%25s"%(SexMatched[i]) + " | " + "%26s"%(SexMismatched[i]) + " | " + "%15s"%(unknowngsex[i]) + " |")
    print("  +" + "-"*91 + "+")
    print()

def read_relate_pairs(infile):
    """
    Read highly related data pairs.
    """
    print("## Reading highly related data pairs ... ")
    for t in omicsSurrogate:
        match[t] = dict()
        matchT[t] = dict()
        mismatch[t] = dict()
        mismatchT[t] = dict()
    neti = 0
    nMatch = 0
    nMismatch = 0
    
    if re.search(r'\.gz$',infile):
        f = gzip.open(infile,'rt')
    else:
        f = open(infile,'r')
    f.readline()
    while 1:
        l = f.readline()
        if not l:
            break
        s = l.rstrip("\n").split("\t")
        t1,s1,t2,s2,relatedness = s[:5]
        k1 = t1+"|"+s1
        k2 = t2+"|"+s2
        
        # Filter unrelated data pairs
        if float(relatedness)<cutoff:
            continue
        
        # Initialize indegree for each sample
        if k1 not in indegree: indegree[k1] = 0
        if k2 not in indegree: indegree[k2] = 0
        
        # connect samples in each sub network
        flag = 0
        for i in idsInNet:
            subnet = idsInNet[i]
            if k1 in subnet:
                idsInNet[i][k2] = s2
                flag = 1
                break
            elif k2 in subnet:
                idsInNet[i][k1] = s1
                flag = 1
                break
        if flag==0:
            idsInNet[neti] = dict()
            idsInNet[neti][k2] = s2
            idsInNet[neti][k1] = s1
            neti += 1
        
        # Samples have at least one related sample
        idsWithConnection.append(k1)
        idsWithConnection.append(k2)
        
        # matched pairs & mismatched pairs
        if s1 not in match[t1]:
            match[t1][s1] = []
            matchT[t1][s1] = []
            mismatch[t1][s1] = []
            mismatchT[t1][s1] = []
        if s2 not in match[t2]:
            match[t2][s2] = []
            matchT[t2][s2] = []
            mismatch[t2][s2] = []
            mismatchT[t2][s2] = []
        if s[5]=="Y":
            match[t1][s1].append(s2)
            match[t2][s2].append(s1)
            matchT[t1][s1].append(t2)
            matchT[t2][s2].append(t1)
            nMatch += 1
        else:
            mismatch[t1][s1].append(s2)
            #mismatch[t2][s2].append(s1)
            mismatchT[t1][s1].append(t2)
            #mismatchT[t2][s2].append(t1)
            nMismatch += 1
    f.close()
    print("  Number of matched pairs: %s"%nMatch)
    print("  Number of mismatched pairs: %s"%nMismatch)
    print()

    ## For example: a sub network contains four nodes A, B, C, and D. If the first round read A & B, the second round read C & D, then they will be splited into two sub networks. We should merge then.
    keydict = dict()
    idsInNet2 = copy.deepcopy(idsInNet)
    for m in idsInNet2:
        for n in idsInNet2[m]:
            if n not in keydict:
                keydict[n] = m
            else:
                ## m and keydict[n] should be merged
                #import pprint
                #pprint.pprint(idsInNet[m])
                #pprint.pprint(idsInNet[keydict[n]])
                for i in idsInNet[keydict[n]]:
                    idsInNet[m][i] = idsInNet[keydict[n]][i]
                del idsInNet[keydict[n]]

def cal_numOmics(xt,x):
    """
    Number of samples in other omics match with the target sample.
    """
    omicsList = []
    for i in range(len(match[xt][x])):
        if omicsSurrogate[matchT[xt][x][i]]!=omicsSurrogate[xt]:
            omicsList.append(omicsSurrogate[matchT[xt][x][i]])
    return len(set(omicsList))

def train_logistic():
    """
    Training parameters in Logistic Regression.
    """
    print("## Reading training set ... ")
    trueset = []
    if re.search(r'\.gz$',train):
        f = gzip.open(train,'rt')
    else:
        f = open(train,'r')
    f.readline()
    while 1:
        l = f.readline()
        if not l:
            break
        s = l.rstrip("\n").split("\t")
        if s[2]=="->":
            trueset.append("%s|%s|%s|%s"%(s[0],s[1],s[3],s[4]))
        elif s[2]=="<-":
            trueset.append("%s|%s|%s|%s"%(s[3],s[4],s[0],s[1]))
        else:
            print("  Warning: Please use \"->\" and \"-<\" in the 3rd column to indicate the directions.")
            continue
    f.close()
    print("  Number of data pairs used in training set: %s"%(len(trueset)))
    
    print("## Preparing parameters for training ... ")
    ## Generate training set
    d = 1
    nOmicsMerge = len(omicsPriority)
    ft = open("%s.train"%output,'w')
    ft.write("difA\tdifB\tdifC\tdirection\n")
    for t in mismatch:
        for m in mismatch[t]:
            mt = t
            #print(mt,m,cal_numOmics(mt,m),nOmicsMerge)
            A1 = cal_numOmics(mt,m)/nOmicsMerge
            for i in range(len(mismatch[t][m])):
                n = mismatch[t][m][i]
                nt = mismatchT[t][m][i]
                #print(nt,n,cal_numOmics(nt,n),nOmicsMerge)
                A2 = cal_numOmics(nt,n)/nOmicsMerge
                B1 = 0
                if nsex[m]==gsex[mt][m]: B1 += 0.5   # nSex(target)==gSex(target)
                if nsex[m]==gsex[nt][n]: B1 += 0.5   # nSex(target)==gSex(source)
                B2 = 0
                if nsex[n]==gsex[nt][n]: B2 += 0.5   # nSex(target)==gSex(target)
                if nsex[n]==gsex[mt][m]: B2 += 0.5   # nSex(target)==gSex(source)
                C1 = omicsPriority[omicsSurrogate[mt]]
                C2 = omicsPriority[omicsSurrogate[nt]]
                difA = A1 - A2
                difB = B1 - B2
                difC = -(C1 - C2)
                if difC!=0: difC = difC/abs(difC)
                
                k1 = "|".join([m,mt,n,nt])
                k2 = "|".join([n,nt,m,mt])
                
                if k1 in trueset:
                    ft.write("%.2g\t%.2g\t%.2g\t%d\n"%(-difA*d,-difB*d,-difC*d,d))
                    d = -1 * d
                elif k2 in trueset:
                    ft.write("%.2g\t%.2g\t%.2g\t%d\n"%(difA*d,difB*d,difC*d,d))
                    d = -1 * d
                
    ft.close()
    
    print("## Training logistic regression model ... ")
    ## Model
    #import matplotlib.pyplot as plt
    import pandas as pd
    from patsy import dmatrices
    
    # Importing the dataset
    dataset = pd.read_table("%s.train"%output)
    # plt.hist(dataset.iloc[:,0].values)
    # plt.hist(dataset.iloc[:,1].values)
    # plt.hist(dataset.iloc[:,2].values)
    
    # read in the data & create matrices
    #y, X = dmatrices('leftward ~ difA + difB + difC', dataset, return_type = 'dataframe')
    X = dataset.iloc[:,:3].values
    y = dataset.iloc[:,3].values
    
    # Training by Logistic Regression
    from sklearn.linear_model import LogisticRegression
    model = LogisticRegression(fit_intercept = False, C = 1e9)
    model.fit(X, y)
    intercept = 0
    #intercept = model.intercept_
    coef = model.coef_
    coef = np.ndarray.tolist(coef)[0]
    print("## Logistic model:")
    print("Intercept: %s"%intercept)
    print("a: %s"%coef[0])
    print("b: %s"%coef[1])
    print("c: %s"%coef[2])
    print()

def judge_directions():
    """
    Judge possible switch directions and probabilities based on Logistic Regression model.
    Write a formated table for Cytoscape input. 
    """
    print("## Judging possible switch directions and probabilities ... ")
    
    nOmicsMerge = len(omicsPriority)
    fc = open("%s.cytoscapeInput"%output,'w')
    fc.write("Omics1\tSampleID1\tOmics_SampleID1\tOmics2\tSampleID2\tOmics_SampleID2\tDirection\tPossibility\n")
    nMis = 0
    nMisDirected = 0
    for t in mismatch:
        for m in mismatch[t]:
            mt = t
            mOmics = cal_numOmics(mt,m)/nOmicsMerge
            for i in range(len(mismatch[t][m])):
                n = mismatch[t][m][i]
                nt = mismatchT[t][m][i]
                nOmics = cal_numOmics(nt,n)/nOmicsMerge
                difA = mOmics - nOmics
                
                nsex_m,nsex_n,gsex_m,gsex_n = nsex[m],nsex[n],gsex[mt][m],gsex[nt][n]
                if nsex_m==0 or nsex_n==0 or gsex_m==0 or gsex_n==0:
                    difB = 0
                else:
                    Bm = 0
                    if nsex_m==gsex_m: Bm += 0.5   # nSex(target)==gSex(target)
                    if nsex_m==gsex_n: Bm += 0.5   # nSex(target)==gSex(source)
                    Bn = 0
                    if nsex_n==gsex_n: Bn += 0.5   # nSex(target)==gSex(target)
                    if nsex_n==gsex_m: Bn += 0.5   # nSex(target)==gSex(source)
                    difB = Bm - Bn
                
                #difC = omicsPriority[omicsSurrogate[mt]] - omicsPriority[omicsSurrogate[nt]]
                difC = -(omicsPriority[omicsSurrogate[mt]] - omicsPriority[omicsSurrogate[nt]])
                if difC!=0: difC = difC/abs(difC)
	 
                ## Scoring
                intercept,a,b,c = coef
                x = intercept + a*difA + b*difB + c*difC
                xx = math.e**x
                p = xx/(1 + xx)
                S = (p - 0.5)*2
                
                ## Judge direction for each pair
                if S>0:
                    fc.write("\t".join([nt,n,"%s|%s"%(nt,n),mt,m,"%s|%s"%(mt,m),"left2right","%.2g"%difA,"%.2g"%difB,"%.2g"%Bm,"%.2g"%Bn,"%.2g"%difC,"%.2g"%S]))
                    #fc.write("\t".join([nt,n,"%s|%s"%(nt,n),mt,m,"%s|%s"%(mt,m),"left2right","%.2g"%S]))
                    indegree["%s|%s"%(mt,m)] += S
                    nMisDirected += 1
                elif S<0:
                    fc.write("\t".join([mt,m,"%s|%s"%(mt,m),nt,n,"%s|%s"%(nt,n),"left2right","%.2g"%(-difA),"%.2g"%(-difB),"%.2g"%(-Bm),"%.2g"%(-Bn),"%.2g"%(-difC),"%.2g"%(-S)]))
                    #fc.write("\t".join([mt,m,"%s|%s"%(mt,m),nt,n,"%s|%s"%(nt,n),"left2right","%.2g"%(-S)]))
                    indegree["%s|%s"%(nt,n)] -= S
                    nMisDirected += 1
                else:
                    #fc.write("\t".join([mt,m,"%s|%s"%(mt,m),nt,n,"%s|%s"%(nt,n),"uncertain","%.2g"%difA,"%.2g"%difB,"%.2g"%difC,"%.2g"%S]))
                    fc.write("\t".join([mt,m,"%s|%s"%(mt,m),nt,n,"%s|%s"%(nt,n),"uncertain","%.2g"%S]))
                nMis += 1
                fc.write("\n")
    
    uniqlist = []
    for t in match:
        for m in match[t]:
            mt = t
            for i in range(len(match[t][m])):
                n = match[t][m][i]
                nt = matchT[t][m][i]
                pair0 = "||".join(sorted(["%s|%s"%(mt,m),"%s|%s"%(nt,n)]))
                if pair0 in uniqlist: continue
                fc.write("\t".join([mt,m,"%s|%s"%(mt,m),nt,n,"%s|%s"%(nt,n),"noSwitch",""]))
                fc.write("\t".join([mt,m,"%s|%s"%(mt,m),nt,n,"%s|%s"%(nt,n),"noSwitch",""]))
                fc.write("\n")
                uniqlist.append(pair0)
    
    fc.close()
    print("  For a total of %s mismatched data pairs, %s pairs with directions estimated."%(nMis,nMisDirected))
    print()

def sort_nodes():
    """
    Determine true sample IDs for samples in each small network by sorting nodes based on indegree and outdegree.
    """
    print("## Determine true sample IDs ... ")
    trueID[1] = []
    trueID[2] = []
    trueID[3] = []
    trueID[4] = []
    for m in idsInNet:
        nDifDataInGroup = len(set(idsInNet[m].values()))
        nDataInGroup = len(idsInNet[m])
        nOmicsInGroup = len(set(map(lambda x: x.split("|")[1],idsInNet[m].keys())))
        if nDifDataInGroup==1 or (stringent==1 and (nDataInGroup<3 or nDifDataInGroup>=nOmicsInGroup)):  ## (left): all samples matched in this sub network. (right): Stringent mode
            for n in idsInNet[m]:
                t,s = n.split("|")
                if gsex[t][s]==0 or nsex[s]==0:
                    sexM = "-"
                elif gsex[t][s]==nsex[s]:
                    sexM = "Y"
                else:
                    sexM = "N"
                trueID[2].append([t,s,s,"N",gsex[t][s],nsex[s],sexM])
            continue
        maxDegree = 0
        for n in idsInNet[m]:
            if indegree[n]>maxDegree:
                maxDegree = indegree[n]
                newID = idsInNet[m][n]
        if maxDegree==0:  # No certain directions for mismatch pairs in this network. Sample IDs in this network are uncertain.
            for n in idsInNet[m]:
                t,s = n.split("|")
                trueID[4].append([t,s,"-","ID_uncertain",gsex[t][s],"-","-"])
        else:
            for n in idsInNet[m]:
                t,s = n.split("|")
                if gsex[t][s]==0 or nsex[newID]==0:
                    sexM = "-"
                elif gsex[t][s]==nsex[newID]:
                    sexM = "Y"
                else:
                    sexM = "N"
                if s==newID:
                    trueID[3].append([t,s,newID,"N",gsex[t][s],nsex[newID],sexM])
                else:
                    trueID[1].append([t,s,newID,"Y",gsex[t][s],nsex[newID],sexM])
                    nSwitch[t] += 1
    noConnection = set(nsex.keys()) - set(idsWithConnection)
    print("  Number of samples switched IDs:")
    for t in nSwitch:
        print("    %s: %s"%(t,nSwitch[t]))
    print()
    
    print("## Writing output files ...")
    fd = open("%s.trueID"%output,'w')
    fd.write("OmicsType\tOriginalID\tTrueID\tSwitchedOrNot\tGeneticSex\tReportedSex_NewID\tSexMatchedOrNot\n")
    for i in trueID:
        for list0 in trueID[i]:
            fd.write("\t".join(list0))
            fd.write("\n")
    fd.close()
    print("  Finished.")
    print()

if __name__ == "__main__":
    read_priority(priorFile)
    read_reported_sex(nsexFile)
    read_genetic_sex(gsexFile)
    read_relate_pairs(pairFile)
    
    ## add train later 
    if train!=0:
        train_logistic()
        sys.exit(1)
    
    judge_directions()
    sort_nodes()


