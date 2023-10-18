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
NSM <- 2e3             # assemble size
RSD <- 1000
out <- "loo"
cache <- mkdir(out, "res")
for(loo in unique(dat$state))
{
    rds <- paste0(FP(cache, loo), ".rds")
    if(file.exists(rds))
    {
        res[[loo]] <- rpt <- readRDS(rds)
        prj <- rpt$prj
        cat(sprintf("%20s: %5d x %3d (<-\"%s\")\n", loo, NROW(prj), NCOL(prj), rds))
        next
    }
    ## divide training and testing - leave one state out.\
    trn <- subset(dat, state != loo)
    tst <- subset(dat, state == loo)

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
            bta <- elm(mdl, loo, mtd="esb", itr)
            ## model quality (internal)
            qua <- cfx(vld[, LHS], predict(mdl, newdata=vld, type='response') > THD)
            ## estimated probability on testing data
            tpr <- predict(mdl, newdata=tst, type='response')
            ret[[itr]] <- list(qua=qua, tpr=tpr, bta=bta)
            j <- j + 1L
        }
        ## bta <- do.call(rbind, lapply(ret, `[[`, "bta"))
        qua <- do.call(rbind, lapply(ret, `[[`, "qua"))
        tpr <- do.call(cbind, lapply(ret, `[[`, "tpr")) # pr(lag=1) of testing data
        tbi <- 0 + (tpr > THD)                          # binary vote?
        bta <- do.call(rbind, lapply(ret, `[[`, "bta"))
        ## weights
        wgt <- within(data.frame(qua), {FPF <- 1-FPF; FNF <- 1-FNF}) |> as.matrix()
        wgt <- cbind(wgt, NUL=1)
        wgt[is.na(wgt)] <- 0.0
        wgt <- sweep(wgt, 2, colSums(wgt), `/`)
        ## projection: weighted enssemble by estimated prob or binary vote
        epr <- melt(tpr %*% wgt, varnames=c("qnm", "wgt"), value.name = "pht")
        ebi <- melt(tbi %*% wgt, varnames=c("qnm", "wgt"), value.name = "pht")
        prj <- rbind(cbind(mtd="epr", epr), cbind(mtd="ebi", ebi))
        prj <- within(prj, yht <- 0 + (pht > mean(trn[, LHS])))
        ## prj <- within(prj, yht <- 0 + (pht > THD))
        prj <- data.frame(loo, prj)
        list(prj=prj, bta=bta)
    })

    ## (b) logistic regression
    lgr <- local(
    {
        mdl <- glm(FRM, 'binomial', trn)
        bta <- elm(mdl, loo, mtd="lgr", itr="NUL")
        prj <- cbind(NUL=predict(mdl, newdata=tst, type="response"))
        prj <- melt(prj, varnames=c("qnm", "wgt"), value.name = "pht")
        prj <- within(prj, yht <- 0 + (pht > mean(trn[, LHS])))
        prj <- data.frame(loo, mtd="lgr", prj)
        list(prj=prj, bta=bta)
    })

    ## (c) informed guess
    igs <- local(
    {
        mtd <- "igs"
        mdl <- glm(IGS, 'binomial', trn)
        bta <- elm(mdl, loo, mtd, itr="NUL")
        prj <- cbind(NUL=predict(mdl, newdata=tst, type="response"))
        prj <- melt(prj, varnames=c("qnm", "wgt"), value.name = "pht")
        prj <- within(prj, yht <- 0 + (pht > mean(trn[, LHS])))
        prj <- data.frame(loo, mtd, prj)
        list(prj=prj, bta=bta)
    })

    ## combined prediction on testing data
    lst <- list(esb, lgr, igs)
    prj <- do.call(rbind, lapply(lst, `[[`, "prj"))
    prj <- within(prj, ref <- tst[qnm, LHS]) # truth
    bta <- do.call(rbind, lapply(lst, `[[`, "bta"))

    ## pack up
    res[[loo]] <- rpt <- list(prj=prj, bta=bta, loo=loo)
    saveRDS(rpt, rds)
    cat(sprintf("%20s: %5d x %3d (<-\"%s\")\n", loo, NROW(prj), NCOL(prj), rds))
}
## combine
prj <- do.call(rbind, lapply(res, `[[`, "prj"))
bta <- do.call(rbind, lapply(res, `[[`, "bta"))
cfg <- data.frame(loo, ZRS, ONS, THD, NSM)

## performance tally
## LOO mean
M14 <- with(dat,
            (sum(lag_2014) - tapply(lag_2014, state, sum))/
            (length(lag_2014) - tapply(lag_2014, state, length)))
tly <- lapply(split(prj, ~ mtd + wgt), with,
{
    ## yht <- 0 + (pht > 1 - M14[loo])
    ## yht <- 0 + (pht > THD)
    if(length(yht) > 0)
        data.frame(mtd=mtd[1], wgt=wgt[1], t(cfx(yht , ref)))
    else
        NULL
})
tly <- data.frame(do.call(rbind, tly), row.names=NULL)
## append AUC
roc <- lapply(split(prj, ~ mtd + wgt), with,
{
    if(length(yht) > 0) pROC::roc(ref, pht, quiet=TRUE) else NULL
})
roc <- roc[!sapply(roc, is.null)]
tly$AUC <- sapply(roc, pROC::auc)
## median betas
mbt <- aggregate(cbind(est, ort=0, zvl, pvl) ~ mtd + bta, bta, median)
mbt <- within(mbt, ort <- exp(est))

## save
tly <- tly[with(tly, order(mtd, wgt)), ]
mbt <- mbt[with(mbt, order(mtd, bta)), ]
write.tsv(prj, file.path(out, "2023_08_08_prj.tsv"))
write.tsv(tly, file.path(out, "2023_08_08_tly.tsv"))
write.tsv(mbt, file.path(out, "2023_08_08_mbt.tsv"))
