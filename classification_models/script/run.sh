## This script is to build PLS classification model to assess characteristics (feature selection) of coding/CDCT/noncoding transcripts

# Dependencies:
# Please install the following packages
# https://github.com/JHwarlock/rklib
# https://github.com/JHwarlock/rblib


## We build two models
# coding vs noncoding
myprefix="../feature_data/coding_noncoding"
# coding vs CDCT
myprefix="../feature_data/coding_CDCT"


# step 1 number of latent varialbles optimazing
plsregress -t category -c1v2 --nlvs=10 --nlvs-max=10 $myprefix.train.feat -s $myprefix.train.info --nlvs-optimize

# step 2 feature importance assessment
plsregress -t category -c1v2 --nlvs=6 --nlvs-max=6 $myprefix.train.feat -s $myprefix.train.info --mcuve --mc-nlvs=6

# step 3 independent testing
plsregress -t category -c1v2 --nlvs=6 --nlvs-max=6 $myprefix.train.feat -s $myprefix.train.info --ex_predict --ex_sampleinfo $myprefix.test.info --ex_matrix $myprefix.test.feat
