#!/usr/bin/env python
import optparse
from Bio import SeqIO


def calcN(seqFile, N=50, more=False):
    seqlist = list()
    totallength = 0
    res = 0
    res_pos = None
    cumsuml = 0
    cumsumh = 0

    with open(seqFile, "rU") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            totallength += len(record.seq)
            seqlist.append(len(record.seq))

    seqlist.sort(reverse=True)

    for entry in seqlist:
        if cumsumh < totallength*(float(N)/100):
            cumsumh += entry
        else: #we are bigger than the treshold, so the last entry was the result
            res_pos = seqlist.index(entry)-1
            res = seqlist[res_pos]
            break

    for entry in seqlist[res_pos+1:]:
        cumsuml += entry

    if more:
        print("N"+str(N)+" = "+str(res) + ", Value found at position "+str(res_pos))
        print( "Number of sequences: " + str(len(seqlist)) + ", total sequence length: " + str(totallength)+\
               ", avg sequence length: " + str( totallength/len(seqlist) ))
        print("Cumulative sum less than N"+str(N)+": " + str(cumsuml)+ " with ratio: " + str(cumsuml/float(totallength)))
        print("Cumulative sum larger than N"+str(N)+": "  + str(cumsumh)+" with ratio: " +str(cumsumh/float(totallength)))

        medianpos = int(len(seqlist)/2)
        if len(seqlist)%2: #odd
            print("Median length is: " +str(  (seqlist[medianpos]+seqlist[medianpos-1])/2))
        else: #even
            print("Median length is: " +str( seqlist[int(len(seqlist)/2)] ))

    else:
        print(res)


def printHelp():
    print "This script expects the filename (or path) to the file to analyze. Currently, only fasta format is supported."
    print "As an optional parameter, you can specify the desired N value\n" \
          "after the file name, \"python calcN.py mySeq.fa 90\" would calculate the N90."
    print "If no value is specified, the N50 will be calculated."
    print "To print additional information, specify --more or -m"

def main():
    usage = "Usage: %prog [options] arg1 (arg2)\narg1 is expected to be the input file in fasta format"\
    +"\narg2 is optional, specifies the desires N Value (default = 50)"
    parser = optparse.OptionParser(usage=usage)
    parser.set_defaults(more=False)
    parser.add_option('--more', '-m', action='store_true', dest='more', help="show additional information")
    #parser.add_option('-m', action='store_true', dest='more')
    (options, args) = parser.parse_args()

    if len(args) > 0:
        if len(args) == 1:
            #run with n=50
            calcN(args[0], 50, options.more)
        elif len(args) == 2:
            #run with provided N
            if(float(args[1]) < 0 or float(args[1]) > 100):
                print ("N has to be between 0 and 100!\n")
                exit(2)
            calcN(args[0], args[1], options.more)
        else:
            print ("Problem with command line arguments!\n")
            printHelp()
    elif len(args) > 4:
        print ("Too many command line arguments!\n")
        printHelp()
        exit(2)
    else:
        printHelp()
        exit(2)


if __name__ == '__main__':
    main()
