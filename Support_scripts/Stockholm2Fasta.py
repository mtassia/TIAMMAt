#!/usr/bin/env python3

from Bio import AlignIO
from Bio import SeqIO
import sys

STOCKHOLM_FILE=sys.argv[1]
OUTPUT_FASTA=open(sys.argv[2],'w')

ALIGNMENT=AlignIO.read(STOCKHOLM_FILE,"stockholm")

for RECORD in ALIGNMENT:
	RECORD.letter_annotations={}
	RECORD.seq=RECORD.seq.ungap('-')
	RECORD.seq=RECORD.seq.ungap('.')
	SeqIO.write(RECORD,OUTPUT_FASTA,"fasta")

