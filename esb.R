## Leave one (state) out analysis
## Authors:
## - Montgomery, Barrett [Michigan State University, montg270@msu.edu]
## - Vsevolozhskaya, Olga [University of Kentucky, ovs222@uky.edu]
## - Tong, Xiaoran [National Institute of Health, xiaoran.tong@nih.gov]
rm(list=ls(all.names=TRUE))
library(ggplot2)
library(pROC)
library(reshape2)
source("hlp.R")
options(width=103)

## helper to extract estimates of a linear model
elm <- function(mdl, ...)
{
    bta <- summary(mdl)$coef
    bta <- data.frame(rownames(bta), bta)
    colnames(bta) <- c("bta", "est", "std", "zvl", "pvl")
    rownames(bta) <- NULL
    cbind(..., bta)
}

## ---------------------------------------------------------
## read-in dataset
## ---------------------------------------------------------
dat <- readRDS("final_data.rds") # 3124 x 85
names(dat) <- tolower(names(dat))
## drop other votes, alternate outcomes, lesser PCs, and SE
dat <- subset(dat, , -c(percentother_2012_votes)) # rm
dat <- subset(dat, , -c(rcl_2012, sens_2014, sens2_2014, sens3_2014))
dat <- subset(dat, , -c(pca3:pca10))
dat <- dat[, -grep("(lower|upper|se)$", names(dat))]
names(dat) <- sub("estimate$", "", names(dat))        # keep estimates
dat <- na.omit(dat)                                   # 3094 x 72
rownames(dat) <- dat[, "qname"]                       # county, state
dat <- within(dat, state <- gsub("^.+, ", "", qname)) # state

## ---------------------------------------------------------
## model specificaion
## ---------------------------------------------------------
NSD <- c("twelve_plus_pmalcprev", "twelve_plus_pmcigprev",
         "twelve_plus_pmmjprev", "twelve_plus_pyaudprev",
         "twelve_plus_pycocprev", "twelve_plus_pysudprev",
         "eighteen_plus_pysmiprev", "eighteen_plus_pystprev") # NSDUH vars
LSE <- paste0(NSD, "sel")                                     # log sd
LMU <- paste0(NSD, "l")                                       # log mu
RHS <- c(NSD, "percentdem_2012_votes", "percentrep_2012_votes", "area_water",
         "area_land", "pca1", "pca2")
LHS <- "lag_2014"
FRM <- paste(LHS, "~", paste(RHS, collapse = " + ")) |> as.formula()
IGS <- paste(LHS, "~", 1) |> as.formula()

## ---------------------------------------------------------
## start rotation
## ---------------------------------------------------------
res <- list()
ZRS <- 1.6; ONS <- 0.8 # legalize to nonlegalize
THD <- ONS / (ZRS+ONS) # threshold to call positive
NSM <- 1e3             # assemble size
RSD <- 1000
out <- "esb"
cache <- mkdir(out, "res")

## check cached results
rds <- FP(out, "res.rds")
if(file.exists(rds))
{
    res <- rpt <- readRDS(rds)
    prj <- rpt$prj
    cat(sprintf("%20s: %5d x %3d (<-\"%s\")\n", loo, NROW(prj), NCOL(prj), rds))
}
trn <- dat

## (a) ensemble
esb <- local(
{
    ret <- list(); j <- 1L; while(j <= NSM) # build NSM weak models
    {
        itr <- SP("%03X", j)
        TZS <- sum(trn[, LHS] == 0) / sum(trn[, LHS] == 1)
        ## bootstrap sub-sampling of training data with replacement
        idx <- c(sample(which(trn[, LHS]==0), sum(trn[, LHS]==1) * ZRS, replace=ONS > 1),
                 sample(which(trn[, LHS]==1), sum(trn[, LHS]==1) * ONS, replace=ONS > 1))
        ## bootstrap sub-sampling of validate data with replacement
        jdx <- c(sample(which(trn[, LHS]==0), sum(trn[, LHS]==1) * ZRS, replace=ONS > 1),
                 sample(which(trn[, LHS]==1), sum(trn[, LHS]==1) * ONS, replace=ONS > 1))
        use <- trn[idx, ] # sub-training
        vld <- trn[jdx, ] # sub-validate
        ## random realization of NSDUH RADS estimates for sub-train and sub-validate
        use[, NSD] <- sgm(rnorm(nrow(use) * length(NSD), unlist(use[, LMU]), unlist(use[, LSE])))
        vld[, NSD] <- sgm(rnorm(nrow(vld) * length(NSD), unlist(vld[, LMU]), unlist(vld[, LSE])))
        mdl <- glm(FRM, 'binomial', use)
        bta <- elm(mdl, mtd="esb", itr)
        ## model quality (internal)
        qua <- cfx(vld[, LHS], predict(mdl, newdata=vld, type='response') > THD)
        ret[[itr]] <- list(qua=qua, bta=bta)
        j <- j + 1L
    }
    ## bta <- do.call(rbind, lapply(ret, `[[`, "bta"))
    bta <- do.call(rbind, lapply(ret, `[[`, "bta"))
    data.frame(zrs=ZRS, ons=ONS, nsm=NSM, bta)
})

## (b) logistic regression
lgr <- local(
{
    mdl <- glm(FRM, 'binomial', trn)
    bta <- elm(mdl, mtd="lgr", itr="NUL")
    data.frame(zrs="NUL", ons="NUL", nsm="NUL", bta)
})

## (c) informed guess
igs <- local(
{
    mtd <- "igs"
    mdl <- glm(IGS, 'binomial', trn)
    bta <- elm(mdl, mtd, itr="NUL")
    data.frame(zrs="NUL", ons="NUL", nsm="NUL", bta)
})

## combine
bta <- data.frame(rbind(esb, lgr, igs), row.names=NULL)

## median betas
mbt <- aggregate(cbind(est, ort=0, zvl, pvl) ~ mtd + bta, bta, median)
mbt <- within(mbt, ort <- exp(est))

## save
mbt <- mbt[with(mbt, order(mtd, bta)), ]
write.tsv(mbt, file.path(out, "2023_10_10_mbt.tsv"))
