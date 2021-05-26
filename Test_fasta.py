#!/usr/bin/env python3
from Bio import SeqIO
import sys

# $1 is the file to be tested; if it is a fasta file, print "TRUE"; if it is not a fasta file, print "FALSE"
FASTA=sys.argv[1]
FASTA=SeqIO.parse(FASTA,"fasta")
print(any(FASTA))
