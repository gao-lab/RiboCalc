## This script is to build PLS classification model to assess characteristict of coding/CDCT/noncoding transcripts

# Dependencies:
# Pliease install the following packages
# https://github.com/zju3351689/rklib
# https://github.com/zju3351689/rblib


## We build two models
# coding vs noncoding
myprefix="../feature_data/coding_noncoding"
# coding vs CDCT
myprefix="../feature_data/coding_CDCT"


# step 1 number of latent varialbles optimazing
~/bin/Mutilanalysis/plsregress -t category -c1v2 --nlvs=10 --nlvs-max=10 $myprefix.train.feat -s $myprefix.train.info --nlvs-optimize

# step 2 feature importance assessment
~/bin/Mutilanalysis/plsregress -t category -c1v2 --nlvs=6 --nlvs-max=6 $myprefix.train.feat -s $myprefix.train.info --mcuve --mc-nlvs=6

# step 3 independent testing
~/bin/Mutilanalysis/plsregress -t category -c1v2 --nlvs=6 --nlvs-max=6 $myprefix.train.feat -s $myprefix.train.info --ex_predict --ex_sampleinfo $myprefix.test.info --ex_matrix $myprefix.test.feat
