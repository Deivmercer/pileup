import pysam
from pysam import AlignmentFile
import re


quality = 2

pysam.index('./sample.bam')
bamfile = AlignmentFile('./sample.bam', 'rb')
all_alignments = bamfile.fetch()
pileup = dict()
for alignment in all_alignments:
    position = alignment.reference_start
    cigar = re.findall("([0-9]+)([MIDNSHP])", alignment.cigarstring)
    sequence = alignment.query_sequence
    qual_string = alignment.qual
    offset = 0
    if 'S' in alignment.cigarstring:
        print("S")
    for aligned_chars, alignment_type in cigar:
        if alignment_type == 'N' or alignment_type == 'H' or alignment_type == 'D':
            position += int(aligned_chars)
        elif alignment_type == 'S':
            position += int(aligned_chars)
            offset += int(aligned_chars)
        elif alignment_type == 'I':
            offset += int(aligned_chars)
        elif alignment_type == 'M':
            for i in range(int(aligned_chars)):
                if (position + i) in pileup:
                    if ord(qual_string[offset]) - 33 >= quality:
                        pileup[position + i].append((alignment.qname, sequence[offset]))
                else:
                    if ord(qual_string[offset]) - 33 >= quality:
                        pileup[position + i] = [(alignment.qname, sequence[offset])]    
                offset += 1
print(pileup[289223])
