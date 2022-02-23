import pysam
from pysam import AlignmentFile
import re


pysam.index('./sample.bam')
bamfile = AlignmentFile('./sample.bam', 'rb')
all_alignments = bamfile.fetch()
pileup = dict()
for alignment in all_alignments:
    position = alignment.reference_start
    cigar = re.findall("([0-9]+)([MIDNSHP])", alignment.cigarstring)
    sequence = alignment.query
    offset = 0
    for aligned_chars, alignment_type in cigar:
        if alignment_type == 'N' or alignment_type == 'S' :
            position += int(aligned_chars)
        elif alignment_type == 'M':
            for i in range(int(aligned_chars)):
                if (position + i) in pileup:
                    pileup[position + i].append((alignment.qname, sequence[offset]))
                else:
                    pileup[position + i] = [(alignment.qname, sequence[offset])]    
                offset += 1
print(pileup[289223])
