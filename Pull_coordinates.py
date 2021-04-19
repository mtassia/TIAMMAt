#!/usr/bin/env python2.7
#Input1=Coordinate_file, #Input2=Fasta_file

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import sys

#Obtain sequences headers and coordinates for extraction
COORD_FILE=open(sys.argv[1], 'r')
COORDINATE_LIST=[]
for COORD_LINE in COORD_FILE:
	SPLIT_LINE=COORD_LINE.split()
	count = 0
	TMP_LIST=[]
	for i in range(len(SPLIT_LINE)):
		if count == 0:
			TMP_LIST.append(SPLIT_LINE[i])
			count += 1
		elif count > 0:
			TMP_LIST.append(int(SPLIT_LINE[i]))
			count += 1
	COORDINATE_LIST.append(TMP_LIST)
COORD_FILE.close()

#Print a helpful message
NUM_COORD=len(COORDINATE_LIST)
print NUM_COORD,'coordinates have been entered.\n'

#For each sequence in fasta file, find matching coordinate data
COUNT=0
CUT_SEQ_RECORDS=[]
for SEQUENCE in SeqIO.parse(sys.argv[2],'fasta',IUPAC.ambiguous_dna):
	SEQ_NAME=SEQUENCE.id
	for COORDINATE in COORDINATE_LIST:	#COORDINATE is a list of [NAME, START, STOP]
		if COORDINATE[0] == SEQ_NAME:
			print 'Processing',SEQ_NAME
			SEQUENCE.description=SEQ_NAME+"_sites_"+str(COORDINATE[1])+"-"+str(COORDINATE[2]) # This line and the following rename sequence id and description
			SEQUENCE.id=SEQ_NAME+"_sites_"+str(COORDINATE[1])+"-"+str(COORDINATE[2])
			if COORDINATE[1] > COORDINATE[2]:	#If coordinates refer to reverse sequence
				if COORDINATE[2] == 1:
					STOP=COORDINATE[1]-1
					CUT_SEQ_RECORDS.append(SEQUENCE[STOP::-1])
				if COORDINATE[2] > 1:
					START=COORDINATE[1]-1
					STOP=COORDINATE[2]-1
					CUT_SEQ_RECORDS.append(SEQUENCE[START:STOP:-1])
			if COORDINATE[2] > COORDINATE[1]:	#If coordinates refer to forward sequence
				START=COORDINATE[1]-1
				STOP=COORDINATE[2]
				CUT_SEQ_RECORDS.append(SEQUENCE[START:STOP])
			COUNT+=1
print 'Total number of sequences processed:',COUNT

#Write output file
FILE_ID=sys.argv[2].split(".")
NEW_FILE_NAME=FILE_ID[0]+".regions_extracted.fasta"
SeqIO.write(CUT_SEQ_RECORDS,NEW_FILE_NAME,"fasta")
