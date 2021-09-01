# s1
~/bin/Mutilanalysis/plsregress  -t category -c1v2 --nlvs=10 --nlvs-max=10 feat.tab -s sample1.dat.random.sel --nlvs-optimize

# s2
~/bin/Mutilanalysis/plsregress  -t category -c1v2 --nlvs=6 --nlvs-max=6  feat.tab  -s sample1.dat.random.sel --mcuve --mc-nlvs=6

# s3
~/bin/Mutilanalysis/plsregress  -t category -c1v2 --nlvs=6 --nlvs-max=6  feat.tab  -s sample1.dat.random.sel --ex_predict --ex_sampleinfo sample1.dat.nosel --ex_matrix dat.nosel.T
