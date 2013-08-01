#!/usr/bin/env python
import optparse, re
from Bio import SeqIO
from itertools import groupby


def findHomopolymers(seqFile, more=False, all=False, length=6):
    cumlen = 0
    num = 0
    results = []
    resulthist = {}
    max = 0
    maxresult = ""
    with open(seqFile, "rU") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if more:
                print "Processing " + str(record.id)
            lists = [list(g) for k, g in groupby(str(record.seq))]
            for x in lists:
                hplen = len(x)
                if(hplen>=length):
                    num += 1
                    results.append(str("Hp"+str(num)+": "+x[0]+" stretch on "+record.id+" at start="+str(cumlen)\
                                       +" and end="+str(cumlen+len(x))+" with length="+str(len(x))))
                    if resulthist.has_key(hplen):
                        resulthist[hplen] += 1
                    else:
                        resulthist[hplen] = 1
                    if(hplen > max): #new max found
                        max = hplen
                        maxresult = results[num-1]
                cumlen += len(x)
    if more:
        print "\nDone!"
        print "Longest stretch: "+maxresult
        print "\nDistribution of homopolymer lengths:"
        print "Length\tCount"
        print "-------------"
        keys = resulthist.keys()
        keys.sort()
        for key in keys:
            print str(key)+"\t"+str(resulthist[key])

    if all:
        for result in results:
            print result

    if not all or more:
        print "\nFound "+str(num)+" homopolymer streches longer or equal "+str(length)

def printHelp():
    print "This script finds homopolymer stretches in a given fasta file."
    print "This script expects the filename (or path) to the file to analyze. Currently, only fasta format is supported."
    print "As an optional parameter, you can specify the desired homopolymer length\n" \
          "after the file name, \"python findHpStretches.py mySeq.fa 8\" would find homopolymers longer or equal to 8bp."
    print "If no value is specified, the a 6bp cutoff will be used for computations."
    print "To print additional information, specify --more or -m"
    print "To print all found streches to stdout, specify --all or -a"

def main():
    usage = "Usage: %prog [options] arg1 (arg2)\narg1 is expected to be the input file in fasta format"\
    +"\narg2 is optional, specifies the desires homopolymer length (default = 6)"
    parser = optparse.OptionParser(usage=usage)
    parser.set_defaults(more=False)
    parser.set_defaults(all=False)
    parser.add_option('--more', '-m', action='store_true', dest='more', help="show additional information")
    parser.add_option('--all', '-a', action='store_true', dest='all', help="print all found stretches")
    #parser.add_option('-m', action='store_true', dest='more')
    (options, args) = parser.parse_args()

    if len(args) == 1:
        #run with n=50
        findHomopolymers(args[0], options.more, options.all)
    elif len(args) == 2:
        #run with provided number
        findHomopolymers(args[0], options.more, options.all, int(args[1]))
    else:
        print ("Problem with command line arguments!\n")
        printHelp()
        exit(2)



if __name__ == '__main__':
    main()
