#!/usr/bin/env python3
# Kim Lauer (klauer2)
# Week 3, Parsing Fasta & Fastq Files
# Due 9/18/16

import argparse
import re
import gzip

# variables
peptide_count = 0
residue_count = 0
# used to temporarily store information for fastq files
t1 = ""
t2 = ""
temp1 = ""
temp2 = ""

# Parsing of Fasta files
def fasta(file_open, residue_count, peptide_count):
    # loop thru file line by line
    for line in file_open:
        # finds number of peptides
        if line.startswith('>'):
            peptide_count += 1
        # finds number of residues
        else: 
            residue_count = residue_count + len(line.rstrip())
    return(residue_count, peptide_count) 

# Parsing Fastq files
def fastq(file_open, residue_count, peptide_count):
    # ASSUMPTION - each entry is four lines 
    line_count = 0
    quality_count = 0
    error = False
    for line in file_open:
        line_count +=1
        # ASSUMPTION - first title line begins with @
        if(line.startswith('@') & (line_count%4==1)):
            # title line should match 3rd line.
            # Title stored as temporary variable
            t1 = re.search(r'(@)(.*)', line)
            temp1 = t1.group(2)
            peptide_count += 1
        # ASSUMPTION - second line is only the amino acid sequence
        elif (line_count%4==2):
            residue_count = residue_count + len(line.rstrip())
        # ASSUMPTION - third line is a repeat of title line starting with +
        elif (line_count%4==3):
            # this should be a repeated title line.  Error message if not.
            if(line.startswith('+')):
                t2 = re.search(r'\+(.*)', line)
                temp2 = t2.group(1)
                if (temp1 != temp2):
                    print("ERROR on line " + str(line_count) + ": title id's do not match")
                    error = True
            else:
                print("Error on line " + str(line_count) + ": repeat title line did not start with +")
                error = True
        # ASSUMPTION - fourth line is the quality codes for sequence data
        elif (line_count%4==0):
            quality_count = quality_count + len(line.rstrip())
            # error message if length of quality data does not match aa length
            if (residue_count != quality_count):
                print("Error: Residue Count on line " + str(line_count-2) + " and Quality Count on line " + str(line_count) + " do not match")
                error = True
                quality_count = residue_count
        # Error message if lines did not follow format descriped above
        else:
            print("Format of data incorrect on line " + str(line_count))
            error = True
        # if errors from above, prints warning
    if(error == True):
        print("Error(s) found in parsing file, please check data to ensure accuracy of results")
    return(residue_count, peptide_count)        

# parses the command line to obtain file name from user
parser = argparse.ArgumentParser()
parser.add_argument("file_name", help="Required - Enter file name with path")
args = parser.parse_args()
file_name = args.file_name

# checks if file is a zip file & opens file in text mode
if (file_name.endswith('.gz')):
    file_open = gzip.open(file_name, 'rt')
else:
    file_open = open(file_name)

# determines if program is a fasta or fastq file.  If neither, error message
if (file_name.endswith('.fa')):
    (residue_count, peptide_count) = fasta(file_open, 0, 0)
elif (file_name.endswith('.fastq')):
    (residue_count, peptide_count) = fastq(file_open, 0, 0)
else:
    print("Error, file name did not end with either .fa or .fastq")

        
# prints results
print("Total sequences found: " + str(peptide_count)) 
print("Total residues found: " + str(residue_count)) 

file_open.close()
