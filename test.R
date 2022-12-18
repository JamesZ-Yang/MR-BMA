# set working directory as your path
setwd("C:/Users/James/demo_AMD")

#download necessary packages
require(knitr)
require(markdown)
require(ggplot2)
require(combinat)
require(hash)
require(corpcor)

#download example sources
source("summary_mvMR_SSS.R")
source("summary_mvMR_BF.R")

#replication 1: all genetic instrumental variables included (n=148, runtime: 14 minutes):
knit("run-BMA-AMD-n148.Rmd")
markdownToHTML('run-BMA-AMD-n148.md', 'run-BMA-AMD-n148.html', options=c("use_xhml"))

#replication 2: after excluding outliers and influential points (n=145, runtime: 11 minutes):
knit("run-BMA-AMD-n145.Rmd")
markdownToHTML('run-BMA-AMD-n145.md', 'run-BMA-AMD-n145.html', options=c("use_xhml"))

#mvMRInput object
amd_nmr_input=new("mvMRInput", betaX = as.matrix(betaX_ivw), betaY = as.matrix(amd_beta_ivw), snps=rs, exposure=rf, outcome = "amd")
BMA_output=summarymvMR_SSS(amd_nmr_input,kmin=1,kmax=12, prior_prob=0.1, max_iter=100)

#combination of risk factors (best.model) output of MR-BMA
best.model.out = sss.report.best.model(BMA_output, top = 10, write.out = TRUE, csv.file.name="amd_best_10models_n145")
best.model.out

#model-averaged results output
mr.bma.out = sss.report.mr.bma(BMA_output, top = 10, write.out = TRUE, csv.file.name="amd_mr_bma_n145")
mr.bma.out

#permutation procedure to calculate empirical p-values (hours of run time)
#nrepeat 10-100 to evaluate the run time; nrepeat = 100,000 for stable results
permute_bma = create.permutations(BMA_output, nrepeat = 1000, save.matrix=TRUE, file.name = "permutation_mrBMA")

empirical.p = calculate.p(BMA_output, permute_bma)
