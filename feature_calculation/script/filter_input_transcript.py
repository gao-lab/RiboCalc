import sys
import gzip
import numpy as np


# record the protein measured protein abundance as y
myprot = open(sys.argv[1],'r') #'../raw_data/OCTOPOS/9606-iBAQ_HEK293_Geiger_2012_uniprot.txt'
prot_abundance = {}
for line in myprot:
    if line.startswith('#'): continue
    line = line.rstrip('\n')
    arra = line.split('\t')
    prot_id = arra[1].split('.')[1]
    ms_abundance = float(arra[2])
    if ms_abundance > 0: prot_abundance[prot_id] = arra[2]
myprot.close()


# record the protein ORF locus from GENCODE annotation
protseq = gzip.open('../raw_data/GENCODE/gencode.v24.pc_translations.fa.gz','rt')
orf_list = []
id_mapping = {}
for line in protseq:
    if line.startswith('>'):
        line = line.lstrip('>')
        arra = line.split('|')
        prot_id = arra[0].split('.')[0]
        if prot_id in prot_abundance:
            orf_id = arra[5]
            orf_list.append(orf_id)
            id_mapping[orf_id] = prot_id
protseq.close()


# translation regulators list
regulator_list = ["MTOR","EIF4G3","EIF2B3","EIF2D","EIF4EBP2","EIF5AL1","NEURL1","NANOS1","EIF3A","EIF3F","EIF4G2","EIF3M","EEF1G",
    "EIF1AD","EIF2S3B","TSFM","DENR","C12orf65","EIF2B1","MTRF1","EIF2S1","EIF2B2","EIF5","CYFIP1","EIF2AK4","EIF3J","RPS27L","IREB2",
    "H3BNC9","CPEB1","GSPT1","EIF3C","TUFM","AARS","C1QBP","EIF5A","EIF4A1","FXR2","SHMT1","RARA","EIF1","IGF2BP1","MRPL58","MIF4GD",
    "EIF4A3","TYMS","CIRBP","EEF2","EIF3G","EIF3K","SAMD4B","RPS9","EIF2B4","EIF2AK2","MTIF2","PAIP2B","EIF2AK3","EIF5B","BOLL","EEF1B2",
    "EIF4E2","EIF2S2","EIF6","EEF1A2","AIRE","EIF3D","EIF3L","GTPBP1","DAZL","TRIM71","EIF1B","EIF4E3","ABTB1","EIF2A","GFM1","EIF5A2",
    "FXR1","EIF2B5","EIF4G1","IGF2BP2","EIF4A2","JAKMIP1","CPEB2","EIF4E","GATB","PAIP1","DHX29","GFM2","DHFR","ETF1","PAIP2","PURA",
    "EIF4EBP3","RPS14","LARP1","CPEB4","EIF4E1B","ABCF1","GTPBP2","EEF1A1","HBS1L","MTRF1L","EIF3B","EIF2AK1","IGF2BP3","EIF4H","EIF4EBP1",
    "COPS5","PABPC1","EIF3E","EIF3H","AGO2","EEF1D","EIF1AX","EIF2S3","RBM3","GSPT2","MCTS1","FMR1","RPL10","EIF1AY","DAZ1","DAZ3","DAZ4"]
regulator_line = "\t".join(regulator_list)
regulator_tpm = {}


# record the transcript TPM calculated by stringtie
mytpm = open(sys.argv[2],'r') #'../raw_data/OCTOPOS/SRR500121.tab'
gene_tpm = {}
for line in mytpm:
    if line.startswith('Gene'): continue
    line = line.rstrip('\n')
    arra = line.split('\t')
    gene_id = arra[0]
    gene_name = arra[1]
    rna_tpm = float(arra[-1])
    if gene_name in regulator_list: regulator_tpm[gene_name] = rna_tpm
    if gene_name == "HPRT1": scale_factor = rna_tpm/109.3212
    if rna_tpm > 0: gene_tpm[gene_id] = rna_tpm
mytpm.close()


# scale the TPM value as log2(TPM+1)
reg_val_list = []
for reg_gene in regulator_list:
    if reg_gene in regulator_tpm:
        rna_tpm = regulator_tpm[reg_gene]
    else:
        rna_tpm = 0
    scaled_tpm = np.log2(rna_tpm/scale_factor+1)
    reg_val_list.append(str(scaled_tpm))
regulator_val = '\t'.join(reg_val_list)


def UTR_length(field):
    info = field.split(':')
    my_type = info[0]
    [start,end] = info[1].split('-')
    my_len = int(end) - int(start) + 1
    return my_type, str(my_len)


def GC_content(seq):
    if len(seq)==0: seq = "A"
    numGC = seq.count("C") + seq.count("G")
    return numGC*1.0/len(seq)


# output the transcript structure and translation regulator features
print('\t'.join(["ProtID","TsID","GeneID","GeneName","Length","CDS_length","X5UTR_length","X3UTR_length",
    "X5UTR_GC","X3UTR_GC","RNA_TPM",regulator_line,"ProtAbundance"]))
myts = gzip.open('../raw_data/GENCODE/gencode.v24.pc_transcripts.fa.gz','rt')
output = sys.argv[3]
seq_file = open(output + '.transcript.fasta', 'w')
cds_file = open(output + '.cds.fasta', 'w')
utr3_file = open(output + '.3utr.fasta', 'w')
init_file = open(output + '.init.fasta', 'w')
for line in myts:
    line = line.rstrip('\n')
    if line.startswith('>'):
        flag = 0
        line = line.lstrip('>')
        arra = line.split('|')
        gene_id = arra[1]
        orf_id = arra[4]
        if orf_id in orf_list and gene_id in gene_tpm:
            flag = 1
            tid = arra[0].split('.')[0]
            gid = gene_id.split('.')[0]
            pid = id_mapping[orf_id]
            gname = arra[5]
            rna_tpm = np.log2(gene_tpm[gene_id]/scale_factor+1)
            rna_tpm = str(rna_tpm)
            ms_abundance = prot_abundance[pid]
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
