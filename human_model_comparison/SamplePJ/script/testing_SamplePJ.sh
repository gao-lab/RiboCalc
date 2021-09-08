# Dependencies
#
# See environment of SamplePJ's model at
# 	https://github.com/pjsample/human_5utr_modeling
# 
# The model is downloaded from 
#	https://github.com/pjsample/human_5utr_modeling/tree/master/modeling/saved_models/main_MRL_model.hdf5
#


# Generate data
# Ensembl_5UTR.raw.fa was fetched from Ensembl biomart with transcript stable IDs
perl filter_seq.pl ../../../RiboCalc/raw_data/coding_features.train.data ../raw_data/Ensembl_5UTR.raw.fa ../raw_data/coding_features.train.RiboTPM ../raw_data/coding_features.train.fasta
perl filter_seq.pl ../../../RiboCalc/raw_data/coding_features.test.data ../raw_data/Ensembl_5UTR.raw.fa ../raw_data/coding_features.test.RiboTPM ../raw_data/coding_features.test.fasta

# Run testing
python ScanOurData_ave.py #Figure S7 (../result/SampleModel_testing.png)
