#!/usr/bin/env python3
#http://biopython-cn.readthedocs.io/zh_CN/latest/en/chr14.html

import subprocess
import sys
import Bio.Seq
import argparse

from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="FastafromBed & Motif search")

parser.add_argument("-f", "--fastagenome",
					help="fasta genome", required=True)
parser.add_argument("-b", "--bed",
					help="bedfile", required=True)
parser.add_argument("-o", "--fastaoutput",
					help="output from fastaFromBed in FASTA format", required=True)
parser.add_argument("-m", "--motifoutput",
					help="output from motif search", required=True)	
parser.add_argument("-s", "--motifsequence",
					help="motif sequence", required=True)	



args=parser.parse_args()

file1 = sys.argv[1] #fasta file
file2 = sys.argv[2] #bed file

fasta_write = open(args.fastaoutput, "w")

output = subprocess.run(['fastaFromBed', '-fi', args.fastagenome, '-bed', args.bed], stdout=subprocess.PIPE )
bytes = output.stdout
stdout=bytes.decode('utf-8')
#print(stdout)

fasta_write.write(stdout)
fasta_write.close()

outputfile = open(args.motifoutput, "w")
 
instances = [Seq(args.motifsequence)]	
m=motifs.create(instances)
myfile= args.fastaoutput

for line in SeqIO.parse(myfile, "fasta"):
 	#print(str(line.seq))
 	motif=Seq(str(line.seq))
 	for pos in m.instances.search(motif):
 		print(line.id, pos[0])
 		outputfile.write("{0}\t{1}\n".format(line.id, pos[0]))
 	
outputfile.close()



# 	test_seq=Seq(line, m.alphabet)
# 	# 
# instances = [Seq("TACAA"), Seq("AATGC")]
# 
# m = motifs.create(instances)
# r = m.reverse_complement(instances)
# 
# #m.alphabet
# #IUPACUnambiguousDNA()
# #IUPAC nucleotide ambiguity codes: W is either A or T, and V is A, C, or G
# # 
# test_seq=Seq("TACACTGCATTACAACCCAAGCATTA",m.alphabet)
# 
# for pos,seq in m.instances.search(test_seq):
# ...     print pos, seq