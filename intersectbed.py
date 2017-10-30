#! /usr/bin/env python3

import subprocess

bed_outputB = open ("../intersect_outputB.txt", "w")

outputB = subprocess.run(['intersectbed', '-a', 'chr1-5_GEM_events.narrowPeak_FUS3', '-b', 'chr1-5_GEM_events.narrowPeak_VRN1', '-wb', '-f', '0.3'], stdout=subprocess.PIPE )
bytesB = output.stdout
stdout_B=bytes.decode('utf-8')
bed_outputB.write(stdout_B)
bed_outputB.close()

bed_outputA = open ("../intersect_outputA.txt", "w")
outputA = subprocess.run(['intersectbed', '-a', 'chr1-5_GEM_events.narrowPeak_FUS3', '-b', 'chr1-5_GEM_events.narrowPeak_VRN1', '-wa', '-f', '0.3'], stdout=subprocess.PIPE )
bytesA = output.stdout
stdout_A=bytes.decode('utf-8')
bed_outputA.write(stdout_A)
bed_outputA.close()

#print(stdout)


