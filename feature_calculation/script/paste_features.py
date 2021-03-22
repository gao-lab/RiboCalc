import sys

output_dict={} #output line of features

#coding features
feat_file = open(sys.argv[1],'r') #../feature_data/coding_feature.temp.tab
p2t_idmapping = {} #protein ID to transcript ID
for line in feat_file:
    line = line.rstrip('\n')
    if line.startswith('ProtID'):
        header = line
    else:
        arra = line.split('\t')
        pid = arra[0]
        tid = arra[1]
        if pid in output_dict and tid.startswith('ENSTR'): continue
        output_dict[pid] = line
        if tid in p2t_idmapping:
            p2t_idmapping[tid].append(pid)
        else:
            p2t_idmapping[tid] = [pid]
feat_file.close()

#AGO binding list
ago_file = open(sys.argv[2],'r') #../feature_data/ago_feature.temp.list
ago_transcript = []
for line in ago_file:
    line = line.rstrip('\n')
    tid = line.split('.')[0]
    for pid in p2t_idmapping[tid]:
        ago_transcript.append(pid)
ago_file.close()


#CDS feature
cds_file = open(sys.argv[3],'r') #../feature_data/cds_feature.temp.tab
cds_list = []
for line in cds_file:
    line = line.rstrip('\n')
    arra = line.split('\t')
    if line.startswith('ID'):
        header_added = '\t'.join(arra[1::])
        header += '\t' + header_added
    else:
        pid = arra[0]
        if pid not in cds_list:
            cds_list.append(pid)
            output_line = '\t'.join(arra[1::])
            output_dict[pid] += '\t' + output_line
cds_file.close()


#MTDR feature
mtdr_file = open(sys.argv[4],'r') #../feature_data/mtdr_feature.temp.tab
mtdr_list = []
header += '\tMTDR\tago_number'
for line in mtdr_file:
    line = line.rstrip('\n')
    arra = line.split(' ')
    pid = arra[0]
    if pid in mtdr_list: continue #remove repeat lines
    mtdr_list.append(pid)
    mtdr = arra[1]
    if pid in ago_transcript:
        ago_num = "1"
    else:
        ago_num = "0"
    output_dict[pid] += '\t' + mtdr + '\t' + ago_num
mtdr_file.close()


#Initiation RNAfold MFE
init_file = open(sys.argv[5],'r') #../feature_data/init_feature.temp.tab
init_list = []
for line in init_file:
    line = line.rstrip('\n')
    arra = line.split('\t')
    pid = arra[0]
    if pid in init_list: continue #remove repeat lines
    init_list.append(pid)
    init_mfe = arra[1]
    if pid == "ID":
        header += '\t' + init_mfe
    else:
        output_dict[pid] += '\t' + init_mfe
init_file.close()


#TargetScan miRNA binding site number of 3'UTR
target_file = open(sys.argv[6],'r') #../feature_data/mirna_feature.temp.tab
targetN_dict = {}
for line in target_file:
    line = line.rstrip('\n')
    arra = line.split('\t')
    targetN_dict[arra[0]] = arra[1]
target_file.close()


#Translation initiation motifs
#Paste all the feature together and output the feature file
pwm_file = open(sys.argv[7],'r') #../feature_data/pwm_feature.temp.tab
for line in pwm_file:
    line = line.rstrip('\n')
    arra = line.split()
    pid = arra[0]
    if pid == "ID":
        header_added = '\t'.join(arra[1::])
        header += '\t' + header_added + '\ttarget_number'
        print(header)
    else:
        output_line = output_dict[pid]
        output_added = '\t'.join(arra[1::])
        if pid in targetN_dict:
            target_num = targetN_dict[pid]
        else:
            target_num = "0"
        output_line += '\t' + output_added + '\t' + target_num
        print(output_line)
pwm_file.close()
