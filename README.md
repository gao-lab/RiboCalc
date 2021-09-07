RiboCalc
====

* 2020-04-30 22:00, Yu-jian Kang
* This is a repository of RiboCalc data and scripts.

1 Download
----
Dependencies:<br>
Model building: R packages - caret, glmnet<br>
Feature caculation scripts: Python packages - biopython, numpy<br>
Feature caculation tools: fimo, RNAfold, TargetScan (miRbase), bedtools, MTDRcalculator, EMBOSS<br>
Human model comparison: Python packages - pdb, pandas, numpy, sklearn, scipy, matplotlib, seaborn


	tom@linux$ git clone git@github.com:gao-lab/RiboCalc.git
  
2 Feature calculation
----
Taking OCTOPOS raw data as an example (Figure S8):<br>
The scripts for feature calculation are at feature_calculation/script<br>
The public datasets are provided in feature_calculation/raw_data<br>
The calculated feature value of OCTOPOS is at feature_calculation/feature_data/

	tom@linux$ cd feature_calculation/script
	tom@linux$ sh test_OCTOPOS.sh

3 RiboCalc model building
----
RiboCalc model building(Figure 2a-c, Figure 3a):<br>
The RiboCalc model data is RiboCalc/RiboCalc.RData<br>
To build RiboCalc without RNA expression, see RiboCalc/remove_RNAexpression.r<br>
The prediction result for TE is TE_testing_result.tab

	tom@linux$ cd RiboCalc
	tom@linux$ Rscript RiboCalc_build_model.r

4 Cell specific model
----
Build cell specific models for the 5 cell lines (Table 1)

	tom@linux$ cd cell_specific_model
	tom@linux$ Rscript cell_specific_model.r

**5 Human model comparison**
---
Performance comparison with LiJJ's human model and SamplePJ's model
	
	tom@linux$ cd human_model_comparison/LiJJ
	tom@linux$ Rscript testing_LiJJ_human.r
	tom@linux$ cd ../SamplePJ/script
	tom@linux$ sh testing_SamplePJ.sh

6 RiboCalc yeast
----
RiboCalc performance testing in yeast (Figure 2d-e, Table 3)

	tom@linux$ cd RiboCalc_yeast
	tom@linux$ Rscript RiboCalc_yeast_model.r

7 Ribo-lncRNA Tesinting
----
RiboCalc prediction of lncRNAs binding with ribosomes reported by previous studies (Figure 3b-c)

	tom@linux$ cd ribo-lncRNA
	tom@linux$ Rscript ribo_lncRNA.r


Contact
----
>If you have any questions about RiboCalc, please mail ribocalc@mail.cbi.pku.edu.cn.

