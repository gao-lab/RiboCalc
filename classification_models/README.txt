RiboCalc is a quantitative model for coding ability prediction. To assess
characteristics for each class of transcripts (coding, CDCTs, and noncoding),
we use feature selection strategy to look into the distinguishing features
between coding vs noncoding RNAs and coding vs CDCTs. Accordingly, we built
two classification models (coding vs CDCT and coding vs noncoding).

We used Partial least squares (PLS) method to construct the models and
select features according to variable importance in projection (>1.5) and
stability of coefficient (>3). The "coding vs noncoding" model's accuracy
was 94.4% in cross validation (n=10,000) and 95.7% in independent testing
(n=18,370). The most important feature for classification was expression
abundance (result/coding_noncoding.selected.features). The sequence features
such as length and CDS length were also distinct between coding and noncoding
transcripts. In the "coding vs CDCT" model, the accuracy was 83.3% in cross
validation (n = 5,000) and 85.0% in independent testing (n=5,804). The most
important 2 features to classify coding and CDCT transcripts were also "RNA
TPM" and "AGO binding" (result/coding_CDCT.selected.features). Different from
"coding vs noncoding", the motifs around start codon (Translation Initiation
motifs, TI motifs) become important, indicating a more suitable sequence
for translation initiation of coding transcripts.

See script/run.sh for details of methods.
