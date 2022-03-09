library(ggplot2)
library(pROC)

set.seed(1000)
## ---------------------------------------------------------
## read-in dataset
## ---------------------------------------------------------
remove(list=ls()); gc()

dat <- readRDS("/final_data.rds")

## remove "Other" votes and alternate specifications for the outcomes
dat <- subset(dat, se=-c(percentother_2012_votes, rcl_2012, sens_2014, sens2_2014, sens3_2014))

## make all names lower cases
names(dat) <- tolower(names(dat))

## drop PCs if desired
dat <- subset(dat, se=-c(pca3, pca4, pca5, pca6, pca7, pca8, pca9, pca10)) 

## drop missing outcomes from analysis
dat <- subset(dat, !is.na(lag_2014)) 
## dat <- within(dat,                   # xt: standardize large numbers
## {
##     area_water <- as.vector(scale(area_water))
##     area_land  <- as.vector(scale(area_land))
## })
dat <- na.omit(dat)                  # xt: 3103 to 3094

## -------------------------------------------------------
## Create macro variables for L and seL
## -------------------------------------------------------
seL <- grep("sel$", names(dat), value = TRUE)              # xt: sd of log
muL <- setdiff(grep("l$", names(dat), value = TRUE), seL)  # xt: mu of log
vnL <- sub('l$', '', muL)                                  # xt: var names
P <- length(vnL)
    
## -------------------------------------------------------
## define train/test sample and case control ratio
## -------------------------------------------------------
## ratio of train-to-test data split
trn <- 0.7; tst <- 0.3
## ratio of legalized-to-nonlegalized to sample
zrs <- 1.6; ons <- .8 

## -------------------------------------------------------
## define objects to store results
## -------------------------------------------------------
var.lst <- list() # data to store variable importance
prd.lst <- list() # data to store county-level predictions
roc.lst <- list() # data to store AUC stuff
acc.lst <- list() # data to store accurary assessment

## ---------------------------------------------------------
## start simulations
## ---------------------------------------------------------
B <- 1e3
for(i in 1:B)
{
    ## bootstrap, training and testing DO NOT share counties.
    tmp <- with(dat,
    {
        i0 <- which(lag_2014 == 0); n0 <- length(i0) # 0s
        i1 <- which(lag_2014 == 1); n1 <- length(i1) # 1s
        ## divide (a) training and (b) testing, then permute
        a0 <- sample(i0, n0 * trn); b0 <- sample(setdiff(i0, a0))
        a1 <- sample(i1, n1 * trn); b1 <- sample(setdiff(i1, a1))
        ## bootstrap, number of 1s as the basic number
        a <- c(rep(a0, len=n1 * zrs * trn), rep(a1, len=n1 * ons * trn))
        b <- c(rep(b0, len=n1 * zrs * tst), rep(b1, len=n1 * ons * tst))
        list(a=sort(a), b=sort(b))
    })
    tmp <- rbind(cbind(dat='a', dat[tmp$a, ]), cbind(dat='b', dat[tmp$b, ]))
    N <- nrow(tmp)
    ## check
    with(tmp, intersect(qname[dat=='a'], qname[dat=='b'])) # must be empty
    with(tmp, table(dat, lag_2014)) # row should be trn:tst, col should be zro:ons

    ## ------------------------------------------------
    ## introduce variability to NSDUH RDAS estimates
    ## ------------------------------------------------
    ## xt: expand standard deviation to 3 times
    .x <- rnorm(N * P, unlist(tmp[, muL]), unlist(tmp[, seL]))
    .x <- matrix(1 / (1 + exp(-.x)), N, P, dimnames=list(NULL, vnL))
    tmp <- cbind(tmp[, !grepl("^(twelve|eighteen)", names(tmp))], .x)

    ## ---------------------------
    ## Perform training testing split
    ## ----------------------------
    A <- subset(tmp, dat == 'a', -dat)
    B <- subset(tmp, dat == 'b', -dat)

    ## ----------------------------
    ## fit a GLM model, predict testing data
    ## ----------------------------
    mld <- glm(lag_2014 ~ ., data    = subset(A, se=-qname), family = 'binomial')
    pht <- predict(mld,      newdata = subset(B, se=-qname), type   = 'response')

    ## ---------------------------------------
    ## store county-level predictions
    ## ---------------------------------------
    ## xt: fraction of 0s as threshold
    res <- cbind(sim=i, B, pht=pht, yht=0 + (pht > mean(B$lag_2014 < 1))) 
    prd.lst[[i]] <- res

    ## -------------------------------------
    ## store predictor statistics
    ## -------------------------------------
    var.lst[[i]] <- data.frame(sim=i, variables=names(coef(mld))[-1], Z=summary(mld)$coef[-1, 3], row.names=NULL)

    ## -------------------------------------
    ## store prediction accuracy results
    ## -------------------------------------
    tbl <- with(res, table(lag_2014, yht))
    acc.lst[[i]] <- data.frame(
        sim=i,
        TP  = tbl[2, 2], TN = tbl[1, 1], FP = tbl[1, 2], FN = tbl[2, 1],
        TPF = tbl[2, 2] / sum(tbl[2, ]), TNF = tbl[1, 1] / sum(tbl[1, ]),
        FPF = tbl[1, 2] / sum(tbl[1, ]), FNF = tbl[2, 1] / sum(tbl[2, ]),
        acc = sum(diag(tbl)) / sum(tbl))
    ## --------------------------
    ## create a ROC curve
    ## --------------------------
    ROC <- with(res, roc(lag_2014, pht, quiet = TRUE))
    roc.lst[[i]] <- with(ROC, data.frame(sim=i, TPR=rev(sensitivities), FPR=rev(1 - specificities)))
}
var2 <- do.call(rbind, var.lst)    # predictor z-scores
acc2 <- do.call(rbind, acc.lst)    # accurary performance
roc2 <- do.call(rbind, roc.lst)    # roc
prd2 <- do.call(rbind, prd.lst)    # truth and prediction, per county, per sim

## ROC for all simulations combined (overall)
rocA <- with(prd2, roc(lag_2014, pht, quiet = TRUE))
aucA <- round(auc(rocA), 3)
rocA <- with(rocA, data.frame(TPR=rev(sensitivities), FPR=rev(1-specificities)))

## ---------------------------------------
## plot all ROC curves with the average
## ---------------------------------------
setwd("C:/Users/montg/Dropbox/Ph.D Work/Dissertation/Data/Aim 1/results/")
g <- ggplot() + theme_bw()
g <- g + geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") # ref
g <- g + geom_line(aes(x=FPR, y=TPR, group=sim), roc2, alpha=0.3) # per sim
g <- g + geom_line(aes(x=FPR, y=TPR), rocA, color="red", size=2)    # overall
g <- g + geom_text(aes(x = 0.75, y=0.4, label=paste0("AUC: ", aucA)), color="red", size=10)
g <- g + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1))
ggsave("roc_glm.png", g)

## -------------------------------------------
## summarize predictive power of variables
## -------------------------------------------
round(do.call(rbind, with(var2, by(Z, variables, summary)))[, 2:5], 3)

#                         1st Qu. Median       Mean 3rd Qu.
# area_land                 0.733  1.346  320390.08   1.932
# area_water                0.111  0.467  228688.61   0.829
# eighteen_plus_pysmiprev   0.527  1.351  441991.16   1.986
# eighteen_plus_pystprev   -1.107 -0.326 -194993.85   0.332
# pca1                     -1.649 -0.715 -195687.63   0.555
# pca2                     -1.588 -0.802 -216793.23  -0.007
# percentdem_2012_votes    -1.619 -0.912 -598545.37  -0.202
# percentrep_2012_votes    -1.782 -1.100 -612526.57  -0.392
# twelve_plus_pmalcprev    -0.004  0.618  268324.83   1.314
# twelve_plus_pmcigprev    -1.567 -0.829 -362955.85  -0.010
# twelve_plus_pmmjprev      2.203  2.871 1202101.03   3.326
# twelve_plus_pyaudprev    -0.085  0.536  392370.13   1.279
# twelve_plus_pycocprev     1.636  2.239  333135.16   2.685
# twelve_plus_pysudprev    -1.320 -0.481   35547.05   0.260

## ------------------------------------------
## Summarize Classification Accuracies
## ------------------------------------------
summary(acc2)
#       sim               TP              TN              FP               FN              TPF        
#  Min.   :   1.0   Min.   : 7.00   Min.   :33.00   Min.   : 0.000   Min.   : 0.000   Min.   :0.3182  
#  1st Qu.: 250.8   1st Qu.:16.00   1st Qu.:40.00   1st Qu.: 2.000   1st Qu.: 3.000   1st Qu.:0.7273  
#  Median : 500.5   Median :17.00   Median :41.00   Median : 3.000   Median : 5.000   Median :0.7727  
#  Mean   : 500.5   Mean   :17.05   Mean   :40.97   Mean   : 3.026   Mean   : 4.945   Mean   :0.7752  
#  3rd Qu.: 750.2   3rd Qu.:19.00   3rd Qu.:42.00   3rd Qu.: 4.000   3rd Qu.: 6.000   3rd Qu.:0.8636  
#  Max.   :1000.0   Max.   :22.00   Max.   :44.00   Max.   :11.000   Max.   :15.000   Max.   :1.0000  
#       TNF              FPF               FNF              acc        
#  Min.   :0.7500   Min.   :0.00000   Min.   :0.0000   Min.   :0.7576  
#  1st Qu.:0.9091   1st Qu.:0.04545   1st Qu.:0.1364   1st Qu.:0.8485  
#  Median :0.9318   Median :0.06818   Median :0.2273   Median :0.8788  
#  Mean   :0.9312   Mean   :0.06877   Mean   :0.2248   Mean   :0.8792  
#  3rd Qu.:0.9545   3rd Qu.:0.09091   3rd Qu.:0.2727   3rd Qu.:0.9091  
#  Max.   :1.0000   Max.   :0.25000   Max.   :0.6818   Max.   :0.9848 

#slightly lower accuracy but much better balance of sens and spec

## ------------------------------------------------
## ensemble predictions per county using different weights
## ------------------------------------------------
WTS <- c('TPF', 'TNF', 'FNF', 'acc')
mrg <- merge(prd2, acc2, by="sim")[, c(names(prd2), WTS)]
ens <- by(mrg, mrg$qname, function(x)
{
    wts <- t(colMeans(x[, 'pht'] * x[, WTS])) # esemble voted probability
    data.frame(
        county = x[1, 'qname'], N = nrow(x), Y = x[1, 'lag_2014'],
        YHT = mean(x[, 'yht']), # mean of per-simulation hard-call
        PHT = mean(x[, 'pht']), # mean of per-simulation preditive probability
        wts)
})
ens <- data.frame(do.call(rbind, ens), row.names=NULL)
## check performance
CRI <- c("YHT", "PHT", WTS) # criteria
round(cor(ens[, 'Y'], ens[, CRI]), 3)

## make hard-calls
THD <- zrs / (zrs + ons) # xt: fraction of 0s as threshold
hdc <- cbind(ens[, !names(ens) %in% CRI], (ens[, CRI] > THD) + 0)

with(hdc, table(Y, TNF))
##    TNF
## Y      0    1
##   0 2751  251
##   1    3   89

with(hdc, table(Y, FNF))
##    FNF
## Y      0
##   0 3002
##   1   92

with(hdc, table(Y, acc))
##    acc
## Y      0    1
##   0 2746  256
##   1    3   89

with(hdc, table(Y, PHT))
##    PHT
## Y      0    1
##   0 2732  270
##   1    1   91

with(hdc, table(Y, YHT))
##    YHT
## Y      0    1
##   0 2734  268
##   1    1   91

## Counties with false positive predictions
subset(ens, Y == 0 & TNF>THD)

## Output for SAS mapping
write.csv(ens, file='C:\\Users\\montg\\Dropbox\\Ph.D Work\\Dissertation\\Data\\Aim 1\\Processed/ens_lag.csv')