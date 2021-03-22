###############################################################
#
# RiboCalc feature calculation and RiboTPM prediction
#
# By Yujian Kang - kangyj@mail.cbi.oku.edu.cn
#
# Take OCTOPOS (Trösemeier et al.) data as an exanple
#
###############################################################


# Input data
# File1 - predicted value: ../raw_data/OCTOPOS/9606-iBAQ_HEK293_Geiger_2012_uniprot.txt
# File2 - genes' RNA abundance file calculated by stringtie: ../raw_data/OCTOPOS/SRR500121.tab
## replace File1 and File2 path to your own files when doing independent testing
prot_abundance='../raw_data/OCTOPOS/9606-iBAQ_HEK293_Geiger_2012_uniprot.txt'
stringtie_output='../raw_data/OCTOPOS/SRR500121.tab' #see stringtie.sh


# Download gene model data from GENCODE
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.pc_transcripts.fa.gz -P ../raw_data/GENCODE
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.pc_translations.fa.gz -P ../raw_data/GENCODE


# Selected expressed RNA for prediction
## replace File1 and File2 path to your own files when doing independent testing
## this script calculate trancript structure and translation regulator features and output sequence files
python filter_input_transcript.py $prot_abundance $stringtie_output ../feature_data/testing > ../feature_data/coding_feature.temp.tab


# Feature ago_number calculation
## humanAGO_*_CLIP.bed are AGO binding sites curated by ANNOLNC (Hou et al.)
par_clip_bed='../raw_data/ANNOLNC/humanAGO_PAR_CLIP.bed'
hits_clip_bed='../raw_data/ANNOLNC/humanAGO_HITS_CLIP.bed'

python transcript_exon_bed.py ../feature_data/coding_feature.temp.tab > ../feature_data/coding_feature.temp.bed
bedtools intersect -a ../feature_data/coding_feature.temp.bed -b $par_clip_bed > ../feature_data/PAR_CLIP_intersect.temp.bed
bedtools intersect -a ../feature_data/coding_feature.temp.bed -b $hits_clip_bed > ../feature_data/HITS_CLIP_intersect.temp.bed

cat ../feature_data/PAR_CLIP_intersect.temp.bed ../feature_data/HITS_CLIP_intersect.temp.bed |cut -f4|sort|uniq > ../feature_data/ago_feature.temp.list
rm -f ../feature_data/coding_feature.temp.bed ../feature_data/PAR_CLIP_intersect.temp.bed ../feature_data/HITS_CLIP_intersect.temp.bed


# ORF feature calculation: MTDR/ENC/CAI/codon frequency
## CAI and ENC were calculated by EMBOSS (Rice et al.) package
## EMBOSS installation or you can calculate CAI and ENC online
tool_bin_path='../bin'
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz -P $tool_bin_path
tar -xvf $tool_bin_path/EMBOSS-6.6.0.tar.gz
cd $tool_bin_path/EMBOSS-6.6.0
./configure
make
cd -
## ORF/CDS feature calculation
python transcript_cds_features.py $tool_bin_path ../feature_data/testing.cds.fasta ../feature_data/cds_feature.temp.tab

## MTDR score file were calculated by MTDRcalculator (Dana and Tuller) (Select organism: H.sapiens-HEK293)
## MTDRcalculator was download at https://www.cs.tau.ac.il/~tamirtul/MTDR/
mtdr_output_file='../feature_data/mtdr_feature.temp.tab'

# Feature init_fold: translation initiation site MFE calculated by RNAfold (Lorenz et al.)
RNAfold < ../feature_data/testing.init.fasta > ../feature_data/testing_init.RNAfold.out
rm -f *_ss.ps
perl RNAfold2MFE.pl ../feature_data/testing_init.RNAfold.out ../feature_data/init_feature.temp.tab


# Feature pwm_*: translation initiation motif identification by FIMO (Grant et al.)
## the translation initiation motifs were reported by TITER (Zhang et al.)
fimo -o ../feature_data/fimo/ ../raw_data/TITER/titer.meme.motif ../feature_data/testing.transcript.fasta
mv ../feature_data/fimo/fimo.txt ../feature_data/testing_pwm.fimo.out
rm -rf ../feature_data/fimo
perl fimo2feature.pl ../feature_data/testing.init.fasta ../feature_data/testing_pwm.fimo.out > ../feature_data/pwm_feature.temp.tab


# Feature target_number: human mature miRNAs target sites number
## scan 3’UTR region against MiRbase () by TargetScan (Grimson et al.)
## raw output processing: count target site number of each query excluding context+ score below -0.2 and "too_close" sites
targetscan_output_file='../feature_data/mirna_feature.temp.tab'


# Paste all the feature data together and do the test
coding_feature_file='../feature_data/coding_feature.temp.tab'
ago_number_file='../feature_data/ago_feature.temp.list'
cds_feature_file='../feature_data/cds_feature.temp.tab'
initiation_mfe_file='../feature_data/init_feature.temp.tab'
titer_pwm_file='../feature_data/pwm_feature.temp.tab'

python paste_features.py $coding_feature_file $ago_number_file $cds_feature_file $mtdr_output_file $initiation_mfe_file \
 $targetscan_output_file $titer_pwm_file > ../feature_data/OCTOPOS_feature.testing.tab


# Clean up the temporary data
# rm -f ../feature_data/*.temp.* ../feature_data/*.fasta


