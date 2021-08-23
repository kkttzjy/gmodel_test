# gmodel_test
The R codes with simulation experiments and real data analysis examples for the two-sample test with g-modeling: two p-value estimation approaches are include -- simple bootstrap (KS) and accelerated bootstrap (ASY). "deconv_new" is the R code for our two-sample test.

## Simulation
The "ZIF simulation" file shows an example code of the similation for the zero-inflated Poisson case. 

## Real data analysis: scRNA-seq data set
The "counts.txt.zip" is scRNA-seq data set used. The "real data ASY.R" and "real data KS.R" are the example code of our two-sample test with two p-value estimation procedures for the single-ccell RNA-seq data set. "real data ASY.R" uses the accelerated bootstrap p-value estimation procedure while "real data KS.R" uses the simple bootstrap p-value estimation procedure

## Real data-based simulation
The "null case.R", "alternative case 1.R" and "alternative case 2.R" are the example codes for the binomial simulation based on the surgical real data set. "null case.R" simulates data under the null hypothesis where the two groups of sample are generated from the same distribution. "alternative case 1.R" and "alternative case 2.R" simulate data under the alternative hypothesis where the two groups of sample are generated from two different distributions.
