import sys
import subprocess
from Bio import SeqIO
from itertools import product


emboss_dir = sys.argv[1]

def ENC(outfile):
    cmd = emboss_dir + "/EMBOSS-6.6.0/emboss/chips -seqall " + outfile + ".tmp.fa -outfile " + outfile + ".chips"
    exitstatus = subprocess.call(cmd, shell=True)
    cmd = "tail -1 " + outfile + ".chips | awk '{print $3}'"
    outtext_enc_tmp = subprocess.check_output(cmd, shell=True)
    outtext_enc = outtext_enc_tmp.decode().rstrip('\n')
    ##add CAI here
    return outtext_enc

def CAI(outfile):
    cmd = emboss_dir + "/EMBOSS-6.6.0/emboss/cai -seqall " + outfile + ".tmp.fa -cfile Ehuman.cut -outfile " + outfile + ".chips"
    #cmd = emboss_dir + "/EMBOSS-6.6.0/emboss/cai -seqall " + outfile + "tmp.fa -cfile Eyeast_cai.cut -outfile " + outfile + ".chips"
    exitstatus = subprocess.call(cmd, shell=True)
    cmd = "awk '{print $NF}' " + outfile + ".chips"
    outtext_cai_tmp = subprocess.check_output(cmd, shell=True)
    outtext_cai = outtext_cai_tmp.decode().rstrip('\n')
    return outtext_cai


def get_codon_frequency(seq):
    codon_freq = {}
    for codon in codon_list:
        codon_freq[codon] = 0
    n = 0
    for i in range(0,len(seq)-2,3):
        codon = seq[i:i+3]
        if codon in codon_freq: codon_freq[codon] = codon_freq[codon] + 1
        n = n + 1
    freq_list = []
    if n==0: n=1
    for codon in codon_list:
        freq = float(codon_freq[codon])/n
        freq_list.append(str(freq))
    line = "\t".join(freq_list)
    return line


'''use itertools to generate the header line of codon name'''
rna_file = open(sys.argv[2],'r')
outfile = sys.argv[3]
out = open(outfile + '.tab','w')
codon_list = []
for i in list(product('ATGC', repeat=3)):
    codon = "".join(i)
    codon_list.append(codon)
codons = "\t".join(codon_list)
header = "ID\tENC\tCAI\t" + codons + "\n"
out.write(header)


seqs = {}
for seq_record in SeqIO.parse(rna_file, "fasta"):
    '''reads FASTA'''
    seqid = seq_record.id
    orfseq = str(seq_record.seq)
    fmrna = open(outfile + '.tmp.fa',"w")
    fmrna.write(">id\n%s\n" % orfseq)
    fmrna.close()
    '''ENC and CAI'''
    enc = ENC(outfile)
    cai = CAI(outfile)
    #print(enc + '\t' + cai)
    frequency = get_codon_frequency(orfseq)
    out.write(seqid + '\t' + enc + '\t' + cai + '\t' + frequency + '\n')
rna_file.close()
out.close()

cmd = 'rm -f ' + outfile + '.tmp.fa ' + outfile + '.chips'
subprocess.call(cmd, shell=True)
