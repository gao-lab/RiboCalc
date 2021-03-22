import sys

tid_list = {}

file = open(sys.argv[1],'r')
for line in file:
    if line.startswith('ProtID'): continue
    line = line.rstrip('\n')
    arra = line.split('\t')
    tid = arra[1]
    tid_list[tid] = 0

bed_file = open('../raw_data/GENCODE/gencode.v24.annotation.exons.bed','r')
for line in bed_file:
    line = line.rstrip('\n')
    arra = line.split('\t')
    tid = arra[3].split('.')[0]
    if tid in tid_list: print(line)
