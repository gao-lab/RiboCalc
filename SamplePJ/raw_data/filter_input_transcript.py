import sys
import gzip
import numpy as np


# record the protein measured protein abundance as y
myprot = open(sys.argv[1],'r') #'../raw_data/OCTOPOS/9606-iBAQ_HEK293_Geiger_2012_uniprot.txt'
prot_list = []
for line in myprot:
    if line.startswith('#'): continue
    line = line.rstrip('\n')
    arra = line.split('\t')
    prot_id = arra[0].split('.')[1]
    prot_list.append(prot_id)
myprot.close()


def UTR_length(field):
    info = field.split(':')
    my_type = info[0]
    [start,end] = info[1].split('-')
    my_len = int(end) - int(start) + 1
    return my_type, str(my_len)



# output the transcript structure and translation regulator features
myts = gzip.open('../raw_data/GENCODE/gencode.v24.pc_transcripts.fa.gz','rt')
utr5_file = open(output + '.5utr.fasta', 'w')
for line in myts:
    line = line.rstrip('\n')
    if line.startswith('>'):
        line = line.lstrip('>')
        arra = line.split('|')
        tid = arra[0].split('.')[0]
        if tid in  prot_list:
                    t_lenth = arra[6]
            length_dict = {"UTR5": "0", "CDS": "0", "UTR3": "0"}
            [mytype, mylen] = UTR_length(arra[7])
            length_dict[mytype] = mylen
            if len(arra) > 9:
                [mytype, mylen] = UTR_length(arra[8])
                length_dict[mytype] = mylen
            if len(arra) > 10:
                [mytype, mylen] = UTR_length(arra[9])
                length_dict[mytype] = mylen
            seq_file.write('>' + pid + '\n')
            cds_file.write('>' + pid + '\n')
            utr3_file.write('>' + pid + '\n')
            init_file.write('>' + pid + '\n')
    else:
        if flag == 1:
            tseq = line
            utr5_len = int(length_dict["UTR5"])
            cds_len = int(length_dict["CDS"])
            utr3_len = int(length_dict["UTR3"])            
            utr5_seq = tseq[0:utr5_len] # 5'UTR sequence
            cds_seq = tseq[utr5_len:utr5_len+cds_len] # CDS sequence
            ntail = cds_len % 3
            if ntail > 0: cds_seq = cds_seq[:-ntail]
            utr3_seq = tseq[-utr3_len::] # 3'UTR sequence
            utr5_gc = GC_content(utr5_seq)
            utr3_gc = GC_content(utr3_seq)
            if utr5_len >= 100:
                init_up = tseq[utr5_len-100:utr5_len]
            else:
                init_up = utr5_seq
            if cds_len + utr3_len >= 103:
                init_down = tseq[utr5_len:utr5_len+103]
            else:
                init_down = tseq[utr5_len::]
            init_seq = init_up + init_down # Upstream and downstream 100 nt sequence around start codon
            seq_file.write(tseq + '\n')
            cds_file.write(cds_seq + '\n')
            utr3_file.write(utr3_seq + '\n')
            init_file.write(init_seq + '\n')
            print('\t'.join(map(str,[pid,tid,gid,gname,t_lenth,cds_len,utr5_len,utr3_len,utr5_gc,utr3_gc,rna_tpm,regulator_val,ms_abundance])))
myts.close()
seq_file.close()
cds_file.close()
utr3_file.close()
init_file.close()
