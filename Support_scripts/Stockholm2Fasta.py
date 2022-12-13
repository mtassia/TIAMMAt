#!/usr/bin/env python3

from Bio import AlignIO
from Bio import SeqIO
import sys

if len(sys.argv) != 3:
	print("Input error.\nUsage is: Stockholm2Fasta.py [Stockholm_file] [Output_fasta]")
	sys.exit() 

STOCKHOLM_FILE=sys.argv[1]
OUTPUT_FASTA=open(sys.argv[2],'w')

ALIGNMENT=AlignIO.read(STOCKHOLM_FILE,"stockholm")

for RECORD in ALIGNMENT:
	RECORD.letter_annotations={}
	RECORD.seq=RECORD.seq.replace('-',"")
	RECORD.seq=RECORD.seq.replace('.',"")
	SeqIO.write(RECORD,OUTPUT_FASTA,"fasta")

