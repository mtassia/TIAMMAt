#CODE DEPENDENCIES:
#	Python3.6+
#	Python3 modules: re, sys, svgwrite

#Initial module import and directory setup
import re
import sys
import svgwrite

#Following block imports necessary data from *besthits.tsv input into a nested list where each entry is a list with the following contents: ["SEQ_ID",SEQ_LENGTH,"DOMAIN_NAME",START,STOP]

#ATOM-SPECIFIC BLOCK FOR BUILDING CODE - LEAVE FOR FUTURE TESTING:
#import os
#print(os.getcwd())
#BEST_HITS=open("Hsapiens_genome_proteins.RHD_DNA_bind_present_after_revision.hmmscan_vs_revised_Pfam.domtblout.besthits.tsv", 'r')
#######################################

#LOAD BESTHITS FILE FROM COMMAND LINE ARGUMENT [1]
BEST_HITS=open(sys.argv[1],'r')

def convert_domtblout_values(string_list): #Create function to convert values of domain table to appropriate string/float elements
	string_list[2] = int(string_list[2])
	string_list[5] = int(string_list[5])
	string_list[9] = int(string_list[9])
	string_list[10] = int(string_list[10])
	string_list[15] = int(string_list[15])
	string_list[16] = int(string_list[16])
	string_list[17] = int(string_list[17])
	string_list[18] = int(string_list[18])
	string_list[19] = int(string_list[19])
	string_list[20] = int(string_list[20])
	string_list[6] = float(string_list[6])
	string_list[7] = float(string_list[7])
	string_list[8] = float(string_list[8])
	string_list[11] = float(string_list[11])
	string_list[12] = float(string_list[12])
	string_list[13] = float(string_list[13])
	string_list[14] = float(string_list[14])
	string_list[21] = float(string_list[21])

LINE_TO_LIST=[] #Initiate list object which will be filled with the necessary fields for domain diagram svg generation
for LINE in BEST_HITS: #read through best hits annotation and format lines appropriately to the hmmscan domain table fields
    SPLIT_LINE=LINE.split() #convert read string into a list object
    DOMTBL_LINE=re.split(r'\s+', LINE.rstrip())	#Initiates an object where fields 0-21 are as designated by HMMers domtblout and 22- are the final field(s)
    DESCRIPTION="" #Initiates a variable to create a single string corresponding to the 'description of target' column in the domtblout format
    for INDEX in range(22, len(DOMTBL_LINE)): #Indicies between 22: correspond to an improperly formated 'description of target' which needs to be converted to a single list index
        DESCRIPTION+=DOMTBL_LINE[INDEX]+" " #Import the 22: indices into a single string
    del DOMTBL_LINE[22:] #Remove the improperly formatted indices from DOMTBL_LINE
    DOMTBL_LINE.append(DESCRIPTION.rstrip()) #Replace the removed indices with the properly formatted index
    convert_domtblout_values(DOMTBL_LINE) #Convert integer and float indices to their appropriate classes
    LINE_TO_LIST.append(DOMTBL_LINE) #After cleaning datalines, append to this variable. LINE_TO_LIST will be a list-of-lists where each primary index is a line from the DOMTBLOUT
BEST_HITS.close()

ANNOTATION_LIST=[] #Initiate a list object to be used for svg diagram creation
for ANNOTATION in LINE_TO_LIST:
    ANNOTATION_DATA=[ANNOTATION[3],ANNOTATION[5],ANNOTATION[0],ANNOTATION[-4],ANNOTATION[-3]]
    ANNOTATION_LIST.append(ANNOTATION_DATA)
#######################################

#Remaineder of code below generates svg graphic from ANNOTATION_LIST list object

#ANNOTATION_LIST=[["SEQ1",300,"DOMAIN_1",50,150],["SEQ1",300,"DOMAIN_2",200,268],["SEQ100200300400500",400,"DOMAIN_1",20,80],["SEQ100200300400500",400,"DOMAIN_2",81,204],["SEQ100200300400500",400,"DOMAIN_200",290,399]] #Hard-coded data for code testing purposes; leave this here in case program breaks

SEQUENCE_IDS=[] #Create a list of unique sequence IDs such that number of indices in this list is the number of annotated sequences
SEQUENCE_LENGTHS=[] #Create a list to obtain sequence lengths which will be used to format the canvas
SEQLEN_DICT={} #Create a dictionary which will be filled with seq_IDs:seq_length pairs

for i in ANNOTATION_LIST: #Obtain basic sequence information
    if i[0] not in SEQUENCE_IDS:
            SEQUENCE_IDS.append(i[0]) # Add new sequence ID to list
            SEQLEN_DICT[str(i[0])]=int(i[1]) #Create a dictionary pair between new sequence ID and it's associated length
    if i[1] not in SEQUENCE_LENGTHS:
            SEQUENCE_LENGTHS.append(i[1]) #Create a list of sequence lengths to be used to format the svg document size

SEQUENCE_LENGTHS.sort(reverse=True) #Sort sequence lengths in descending order to find the longest sequence length
LONGEST_SEQ=SEQUENCE_LENGTHS[0] #Use the longest sequence to generate the width of the svg canvas
NUM_SEQS=len(SEQUENCE_IDS) #Use the number of sequences to generate the length of the svg canvas
LONGEST_ID_LEN=len(max(SEQUENCE_IDS, key=len)) #Create an integer variable which contains the length of the longest string - this is to be used to appropriately format the width of the canvas.
LINE_SCAFFOLD_START=(LONGEST_ID_LEN*10)+5

#CREATE SVG CANVAS
svgdoc=svgwrite.Drawing(
    filename="test_drawing.svg",
    size = (((LONGEST_ID_LEN*10)+LONGEST_SEQ+(LONGEST_SEQ*(1/20))),NUM_SEQS*30)) #At size 10 arial bold, the W glyph is 9.43px -- therefor the x for scaffold size is loaded as a function of (longest length sequence name * 10) + (length of longest sequence) + (Addition 5% of the length of the longest seq to account for trim). Y scaffold size is loaded as the number of sequences * 30 (15 above and below scaffold line) for appropriate spacing. Final 5 accounts for the spacing between text and the line object.

#Create TEXT AND LINE-SCAFFOLD PER SEQUENCE
for SEQi in range(len(SEQLEN_DICT)): #Iterate numerically through a series of numbers equal to the number of unique sequence
    svgdoc.add(svgdoc.text( #Create the text for the sequence scaffold
        str(SEQUENCE_IDS[SEQi]), #Load text
        insert=((LONGEST_ID_LEN*10),((SEQi*30)+15)), #Place cursor
        stroke='none', #No stroke on text
        fill=svgwrite.rgb(0,0,0), #Write letters in black with no stroke
        alignment_baseline="middle", #Positioning baseline is set to the middle of the written text
        text_anchor='end', #Text is written so that the last character is 5 pixels from the start of the line scaffold (below)
        font_size='10px',
        font_weight='bold',
        font_family='Arial'
    ))
    svgdoc.add(svgdoc.line( #Create the line for the sequence scaffold
        start=(LINE_SCAFFOLD_START,((SEQi*30)+15)), #Start of the line scaffold is equal to the LONGEST_ID_LEN * 10 + 5
        end=((LINE_SCAFFOLD_START)+SEQLEN_DICT[SEQUENCE_IDS[SEQi]],((SEQi*30)+15)),
        stroke='black',
        stroke_width='1px'
    ))
    svgdoc.add(svgdoc.text( #Create a text object which is the scaffold length and place it at the end of the line
        (str(SEQLEN_DICT[SEQUENCE_IDS[SEQi]])+'aa'), #Load sequence length
        insert=((LINE_SCAFFOLD_START)+SEQLEN_DICT[SEQUENCE_IDS[SEQi]]+5,((SEQi*30)+15)), #Insert cursor at the end of the line object + 5
        stroke='none',
        fill=svgwrite.rgb(0,0,0),
        alignment_baseline="middle",
        text_anchor='start',
        font_size='6px',
        font_family='Arial'
    ))

#Iterate through annotations per sequence to make domain overlays
COUNT=0 #Create a count variable which will help with y-axis formatting as the annotations are looped
for ID in SEQUENCE_IDS: #Loop through a list of the sequence IDs
    for ANNOTATION in ANNOTATION_LIST: #Loop through entries in the annotation list
            if ID == ANNOTATION[0]: #If the annotation line corresponds to the sequence in the list of sequences being iterated above, do the following
                DOMAIN_NAME=ANNOTATION[2] #Load the domain name
                DOMAIN_START=ANNOTATION[3] #Load the start coordinate for the domain
                DOMAIN_END=ANNOTATION[4] #Load the end coordinate for the domain
                DOMAIN_MIDPOINT=(LINE_SCAFFOLD_START+((DOMAIN_START+DOMAIN_END)/2)) #Create a variable with the X coordinate of the domain shape midpoint
                svgdoc.add(svgdoc.rect( #Create a rectangle overlaying the line that corresponds to each domain
                    insert=((LINE_SCAFFOLD_START+DOMAIN_START),((COUNT*30)+7.5)),
                    size=(((DOMAIN_END)-(DOMAIN_START)),15),
                    rx=2, #Make rectangles with rounded edges
                    fill="white",
                    stroke="black",
                    stroke_width="1px"
                ))
                svgdoc.add(svgdoc.text( #Create the lable for the domain ID within the domain shape
                    str(DOMAIN_NAME), #Load domain name text
                    insert=((DOMAIN_MIDPOINT),((COUNT*30)+15)), #Place cursor
                    stroke='none', #No stroke on text
                    fill=svgwrite.rgb(0,0,0), #Write letters in black with no stroke
                    alignment_baseline="middle", #Positioning baseline is set to the middle of the written text
                    text_anchor='middle', #Text is written from the middle position
                    font_size='6px',
                    font_family='Arial'
                ))
    COUNT+=1

#SAVE DRAWING
FILE_OUTPUT_NAME=(str(sys.argv[1])+".domain_diagram.svg")
svgdoc.saveas(FILE_OUTPUT_NAME)

#TEST OUTPUT - SAVE FOR FUTURE TESTING
#svgdoc.saveas("test_drawing.svg") #LINE FOR TESTING PURPOSES - CAN DELETE
