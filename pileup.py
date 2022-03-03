import sys
import getopt
import pysam
from pysam import AlignmentFile
import re
from collections import Counter


def phred_to_int(quality):
    return ord(quality) - 33


quality_threshold = queried_position = -1
try:
    opts, args = getopt.getopt(sys.argv[1:], "q:p:", ["quality-threshold=", "position="])
except getopt.GetoptError:
    sys.exit()
for opt, arg in opts:
    if opt in ("-q", "--quality-threshold"):
        quality_threshold = int(arg)
    elif opt in ("-p", "--position"):
        queried_position = int(arg)
if quality_threshold == -1:
    print("Missing quality threshold argument.")
    print("Usage: python pileup.py --quality-threshold=<number> --position=<number>")
    sys.exit(-1)
if queried_position == -1:
    print("Missing position argument.")
    print("Usage: python pileup.py --quality-threshold=<number> --position=<number>")
    sys.exit(-1)

pysam.index('./sample.bam')
bamfile = AlignmentFile('./sample.bam', 'rb')
all_alignments = bamfile.fetch()
pileup = dict()
for alignment in all_alignments:
    if alignment.cigarstring is None:
        continue
    position = alignment.reference_start
    cigar = re.findall("([0-9]+)([MIDNSHP])", alignment.cigarstring)
    sequence = alignment.query_sequence
    qual_string = alignment.qual
    offset = 0
    for aligned_chars, alignment_type in cigar:
        if alignment_type == 'H' or alignment_type == 'P':
            continue
        if alignment_type == 'N' or alignment_type == 'D':
            position += int(aligned_chars)
        elif alignment_type == 'S' or alignment_type == 'I':
            offset += int(aligned_chars)
        elif alignment_type == 'M':
            for _ in range(int(aligned_chars)):
                quality = phred_to_int(qual_string[offset])
                if position in pileup:
                    if quality >= quality_threshold:
                        pileup[position].append((alignment.qname, sequence[offset], quality))
                else:
                    if quality >= quality_threshold:
                        pileup[position] = [(alignment.qname, sequence[offset], quality)]
                position += 1
                offset += 1

read_list = pileup[queried_position]
count = Counter([read[1] for read in read_list])
print(count)
print("Number of reads in the queried position: " + str(len(read_list)))
for read in read_list:
    print("Query name: " + read[0] + "\t", "\tBase: " + read[1], "\tQuality: " + str(read[2]))
