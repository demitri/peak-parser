#!/usr/bin/env python3
#http://biopython-cn.readthedocs.io/zh_CN/latest/en/chr14.html

import subprocess
import sys
import Bio.Seq

from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq

file1 = sys.argv[1] #fasta file
file2 = sys.argv[2] #bed file

output = subprocess.run(['fastaFromBed', '-fi', file1, '-bed', file2], stdout=subprocess.PIPE )
bytes = output.stdout
stdout=bytes.decode('utf-8')
#print(stdout)

fasta_write = open ("TAIR_fasta_output", "w")
fasta_write.write(stdout)
fasta_write.close()

#fasta_frombed = open ("TAIR_fasta_output", "r")

instances = [Seq("TACAA")]
m=motifs.create(instances)

myfile= "TAIR_fasta_output"

for line in SeqIO.parse(myfile, "fasta"):
	#print(str(line.seq))
	motif=Seq(str(line.seq))
	for pos in m.instances.search(motif):
		print(line.id, pos[0])
	
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