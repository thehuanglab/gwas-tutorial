# this package can be downloaded from https://github.com/hailianghuang/gwas-tutorial
# contact: Hailiang Huang (hhuang@broadinstitute.org)
# install PLINK2 from https://www.cog-genomics.org/plink2
# install R from https://www.r-project.org
# install the 'qqman' package in R.

#data inspection
head dat/gwas_example.ped | cut -c 1-100 
head dat/gwas_example.map

#convert to binary files
plink --file dat/gwas_example --make-bed --out gwas_example
ls -l 
head gwas_example.bim 
head gwas_example.fam

ls -rlth gwas_example.bed
ls -rlth dat/gwas_example.ped

#calculate missing rate
plink --bfile gwas_example --missing --out gwas_example
head gwas_example.imiss
awk '$6 > 0.02' gwas_example.imiss

#step 4: calculate heterozygosity rate
plink --bfile gwas_example --het --out gwas_example
head gwas_example.het

#Call the R code to generate a list of IDs failing the missing rate and heterozygosity test
Rscript gwas_remove_imiss_het.R

#remove the individuals in the list
plink --bfile gwas_example --remove gwas_example.imiss-vs-het.remove  --make-bed --out gwas_example_qc1

#variants failing the missing rate test
head gwas_example.lmiss
awk 'NR>1 && $5 > 0.05 {print $2} ' gwas_example.lmiss > gwas_example.lmiss.exclude
plink --bfile gwas_example_qc1 --make-bed --exclude gwas_example.lmiss.exclude --out gwas_example_qc2

#variants failing the differential missing rate test
plink --bfile gwas_example_qc2 --test-missing --out gwas_example_qc2
head gwas_example_qc2.missing
awk '$3-$4 > 0.05 || $3-$4 < -0.05 {print $2}' gwas_example_qc2.missing > gwas_example_qc2.missing.exclude
plink --bfile gwas_example_qc2 --exclude gwas_example_qc2.missing.exclude --make-bed --out gwas_example_qc3


#variants failing the HWE test
plink --bfile gwas_example_qc3 --hardy --out gwas_example_qc3 
head gwas_example_qc3.hwe
awk '$3=="UNAFF" && $9 < 0.000001 {print $2}' gwas_example_qc3.hwe > gwas_example_qc3.hwe.exclude
plink --bfile gwas_example_qc3 --exclude gwas_example_qc3.hwe.exclude --make-bed --out gwas_example_qc4

#LD pruning for removing related individuals
plink --bfile gwas_example_qc4 --indep-pairwise 200 5 0.2 --out gwas_example_qc4
head gwas_example_qc4.prune.in
plink --bfile gwas_example_qc4 --extract gwas_example_qc4.prune.in --make-bed --out gwas_example_qc4_indep
plink --bfile gwas_example_qc4_indep --genome --out gwas_example_qc4_indep
head gwas_example_qc4_indep.genome
awk 'NR>1&& $10 > 0.2 {print $1"\t"$2}' gwas_example_qc4_indep.genome > gwas_example_qc4_indep.genome.remove
plink --bfile gwas_example_qc4 --remove gwas_example_qc4_indep.genome.remove --make-bed --out gwas_example_final

wc -l gwas_example.bim 
wc -l gwas_example.fam

wc -l gwas_example_final.bim 
wc -l gwas_example_final.fam

#generate PC for population structure (using the new dataset without related individuals)
plink --bfile gwas_example_final --indep-pairwise 200 5 0.2 --out gwas_example_final
plink --bfile gwas_example_final --extract gwas_example_final.prune.in --make-bed --out gwas_example_final_indep
plink --bfile gwas_example_final_indep --pca --out gwas_example_final_pca

column -t gwas_example_final_pca.eigenvec | less -S

#GWAS!
plink --bfile gwas_example_final --logistic hide-covar --covar gwas_example_final_pca.eigenvec --covar-number 1-10 --out gwas_example_final  --ci 0.95

head gwas_example_final.assoc.logistic

#make figure 
Rscript gwas_make_figure.R
