#!/usr/bin/env python

import optparse, time, sys
from Bio import SeqIO


def processFile(seqFile, cutoff=500, more=False, stream=False):
    start_time = time.time()

    if(more):
        print("Processing File "+str(seqFile))
        print("Output will be written to "+outFile)
        print("Using cutoff length "+str(cutoff))

    input_seq_iterator = SeqIO.parse(open(seqFile, "rU"), "fasta")
    short_seq_iterator = (record for record in input_seq_iterator if len(record.seq) > float(cutoff))

    if stream:
        output_handle = sys.stdout
    else:
        splitname =  seqFile.split(".") #generate the modified file name
        splitname[-2] = splitname[-2]+"_filtered"
        outFile = ".".join(splitname)
        output_handle = open(outFile, "w")

    SeqIO.write(short_seq_iterator, output_handle, "fasta")

    if not stream:
        output_handle.close()

    if(more):
        print "Done! Processing took "+ str(time.time() - start_time)[:6]+ " seconds"
    elif not stream:
        print "Done!"



def printHelp():
    print "This script expects the filename (or path) to the file to process. Currently, only fasta format is supported."
    print "As an optional parameter, you can specify the desired cutoff length " \
          "after the file name, \"python filterShort.py mySeq.fa 90\" would filter sequences shorter than 90bp."
    print "If no value is specified, a treshold of 500 will be used."
    print "To print additional information, specify --more or -m"
    print "To print to stdout, specify --stream or -stream. If -s is set, no output file will be generated"

def main():
    usage = "Usage: %prog [options] <input file> <treshold>\nInput file should be in FASTA format"\
    +"\nThreshold is optional, specifies the desires cutoff Value (default = 500)"
    parser = optparse.OptionParser(usage=usage)
    parser.set_defaults(more=False)
    parser.set_defaults(stream=False)
    parser.add_option('--more', '-m', action='store_true', dest='more',
                      help="show additional information (this should not be used together with -s)")
    parser.add_option('--stream', '-s', action='store_true', dest='stream',
                      help="output to stdout instead of creating an output file")
    #parser.add_option('-m', action='store_true', dest='more')
    (options, args) = parser.parse_args()

    if len(args) > 0:
        if len(args) == 1:
            #run with threshold = 500
            processFile(args[0], 500, options.more, options.stream)
        elif len(args) == 2:
            #run with provided cutoff
            processFile(args[0], args[1], options.more, options.stream)
        else:
            print ("Problem with command line arguments!\n")
            printHelp()
    elif len(args) > 4:
        print ("Too many command line arguments! Ignoring some...\n")
        processFile(args[0], args[1], options.more, options.stream)
    else:
        printHelp()
        exit(2)


if __name__ == '__main__':
    main()