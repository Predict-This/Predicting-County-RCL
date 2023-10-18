## TCOM Tool related helpers
HLP <- new.env()

## default values
HLP$NAS <- c("", "NA", "N/A", "NULL", NA) # na-strings

#' Function Composition
#'
#' Given f(x, a) and g(., b), compose function f o g: g(f(x, a), b).
#'
#' @details
#' The composed function
#' * supports partial arguments assignment like [base::pipeOp()].
#' * ensures environment sensitive function (e.g., [ls()]) acting normally.
#' * encloses syntax of cascated component function call in its body.
#'
#' @usage
#' f(a) %*% g(b)
#'
#' @param f(x, a), g(., b) component functions with partial arguments.
#' @return composite function g(f(x, a), b) of `x`
#'
#' @examples
#' ## component functions for testing
#' ff <- function(n, a=1.0) 1:(n*a)
#' gg <- function(x, b=0.5) x^b
#' hh <- function(x, c=2.0, d=0.5) d * log(x, c)
#' 
#' example #1, rotate `c` in hh()
#' (tmp <- ff(a=1) %.% gg(b=2) %.% hh(x=3, d=1) %.% round(digits=3))
#' tmp(9)
#'
#' example #2: rotate `x` in hh(); anonymous function usage.
#' sapply(1:4, ff(a=1) %.% gg(b=2) %.% hh(d=1, c=2) %.% round(digits=3))
#'
#' example #3: one of the component function is environment sensitive.
#' (ls() %.% grep(pattern="^[a-z]", value=TRUE) %.% toupper)(environment())
#' should not show `f`, `g` and `ret` (and `.dp`) since the composition was
#' "supposedly" done outside of `%.%`'s environment.
HLP$"%.%" <- function(f, g, .dp=0)
{
    f <- substitute(f)
    g <- substitute(g)

    ## treatment of f()
    if("%.%" == all.names(f)[1])
    {
        ## f() is itslef a %.% expression, expand it by a deeper call.
        f <- eval(as.call(append(as.list(f), c(.dp=.dp+1))))
    }
    else
    {
        ## f() is  the deepest function, (i.e.,  the left most in the
        ## %.% chain), insert ... as its first argument.
        f <- as.call(append(as.list(f), quote(...), 1L))
    }

    ## expands g() to g(f(...)) by treating f(...) as the 1st argument
    g <- as.call(append(as.list(g), f, 1L))

    ## pack up and return
    if(.dp > 0)
    {
        ret <- g # a deeper call of %.% return expanded expression
    }
    else
    {
        ## the top call of %.% build a function enclosing g(f(...))
        ret <- function(...) {}
        body(ret) <- g
        ## treat the function as if it was defined from outside.
        environment(ret) <- parent.frame()
    }
    ret
}

#' assign if not already exists
HLP$"%:-%" <- function(x, y) 
{
    Var <- deparse(substitute(x))
    if (!exists(Var, parent.frame()) || length(x) == 0L)
    {
        assign(Var, y, parent.frame())
    }
    invisible(NULL)
}

#' y if x is none.
HLP$"%||%" <- function(x, y) if(length(x)) x else y

#' x if not y.
HLP$"%&!%" <- function(x, y) if(y) NULL else x

#' x if y.
HLP$"%&&%" <- function(x, y) if(y) x else NULL


#' make and return deep directories without warning
HLP$mkdir <- function(...)
{
    d <- file.path(...)
    dir.create(d, FALSE, TRUE)
    d
}

#' cache evaluation 
#'
#' Idealy, a cached expression should only be evaluated once and the future call
#' return the cached result.
#'
#' @param .rds the R dataset (*.rds) to store the cache
#' @param .rxp the R expression to evaluate.
#' @param over overwrite existing cache? (def = N)
#' @param here evaluate {.rxp} locally? (def=N)
#' @param pack automatically pack up new / changed objects?
#'
#' Bewere not  to write R expression  relying on external objects  named {over},
#' {here}, or {pack}.
#' @examples
#' r <- cache("ex1.rds",
#' {
#'     a <- rnorm(9)
#'     b <- rnorm(9)
#'     a * b
#' }, over=0)
HLP$cache <- function(.rds, .rxp, over=0, here=0, pack=0)
{
    out <- parent.frame() # outside
    env <- new.env()      # sandbox
    if(file.exists(.rds) && !over) {
        ret <- readRDS(.rds)
    } else {
        ret <- eval(substitute(.rxp), env, out)
        if(pack)
            ret <- as.list(env)
        saveRDS(ret, .rds)
    }
    if(here == 0)
    {
        for(n in names(ret))
            assign(n, ret[[n]], out)
    }
    invisible(ret)
}

#' flood objects from a list into an environment
#'
#' Unpack objects from a list into an environment.
#'
#' @param pck the container of newer objects.
#' @param env the environment (def=<caller>).
HLP$unpack <- function(pck, env=NULL, vbs=0)
{
    env %:-% parent.frame()
    enm <- environmentName(env)
    pnm <- as.character(substitute(pck))
    if(vbs)
        PL("unpack obj from \"%s\" into env [%s]:", pnm, enm)
    for(. in names(pck))
    {
        env[[.]] <- pck[[.]]
        if(vbs)
            PL("- %s: %7d x %5d, %16s", ., NROW(env[[.]]), NCOL(env[[.]]),
               format(as.integer(object.size(env[[.]])), big.mark=","))
    }
    invisible(NULL)
}
HLP$upk <- HLP$unpack

#' Repack objects from an environment into a list.
#'
#' @param pck the container of older objects.
#' @param env the environment (def=<caller>).
HLP$repack <- function(pck, env=NULL, vbs=0)
{
    env %:-% parent.frame()
    enm <- environmentName(env)
    pnm <- as.character(substitute(pck))
    ret <- list()
    if(vbs)
        PL("repack obj from env [%s] into \"%s\":", enm, pnm)
    for(. in names(pck))
    {
        ret[[.]] <- env[[.]]
        if(vbs)
            PL("- %s: %7d x %5d, %16s", ., NROW(env[[.]]), NCOL(env[[.]]),
               format(as.integer(object.size(env[[.]])), big.mark=","))
    }
    ret
}
HLP$rpk <- HLP$repack

#' Write TSV file.
#'
#' A wrapper of  [write.table()] that always uses "\t" as  separator, by default
#' uses no quotation and blank for NA.
HLP$write.tsv <- function(x, f, quote=FALSE, na="", ...)
{
    #' wrapper of write.table(...)
    #' - sep="\t",
    #' - row.names=FALSE,
    #' - na="",
    #' - quote=FALSE
    write.table(x, f, quote=quote, sep="\t", na=na, row.names=FALSE, ...)
}
HLP$saveTSV <- HLP$Pwrite.tsv

#' Read TSV file.
#'
#' A wrapper of [read.delim()], use no quotation and treat blank as NA.
HLP$read.tsv <- function(f, ...)
{
    #' wrapper of read.delim(...)
    #' - na.strings="",
    #' - quote=NULL,
    #' - check.names=FALSE
    read.delim(f, na.strings="", quote=NULL, check.names=FALSE, ...)
}
HLP$readTSV <- HLP$read.tsv

#' split {x} by {g}, apply {f}, and unsplit.
HLP$xgf <- function(x, g, f, ...) unsplit(lapply(split(x, g), f, ...), g)

#' split {x} by {g}, tabluate unique values.
#'
#' @param na TRUE to count NA (def=1), or a string to name the count.
HLP$xgt <- function(x, g, na=TRUE)
{
    x <- as.factor(x)
    y <- levels(x)
    r <- outer(x, y, `==`)
    colnames(r) <- y
    if(!is.null(na) && na != FALSE)
    {
        if(na == TRUE)
            na <- NA
        r <- cbind(r, matrix(is.na(x), dimnames=list(NULL, na)))
    }
    rowsum(0 + r, g, na.rm=TRUE)
}

#' lengh of unique values in {x}
HLP$lux <- function(x, na.rm=FALSE)
{
    #' length of unique values
    if(na.rm)
        x <- x[!is.na(x)]
    length(unique(x))
}

#' difference between values in {x} padded with an initial value.
#'
#' A wrapper of [diff()] which pads an  initial value (def=0) to the results, so
#' it has the same length with {x}.
HLP$xdf <- function(x, ini=0) c(ini, diff(x))

#' emulated printf
HLP$PF <- function(fmt, ...)
{
    msg <- if(missing(fmt)) "" else sprintf(fmt, ...)
    cat(msg)
    invisible(c(msg, ..., recursive=FALSE))
}

#' emulated printf and with new lines
HLP$PL <- function(fmt, ...)
{
    msg <- if(missing(fmt)) "" else sprintf(fmt, ...)
    cat(msg, "\n", sep="")
    invisible(c(msg, ..., recursive=FALSE))
}
HLP$SP <- sprintf

#' emulated file path
HLP$FP <- file.path

#' praint a horizontal line
HLP$HL <- function(ch="-", bg=ch, ed=ch, wd=NULL, sep="")
{
    if(is.null(wd))
        wd <- options()$width
    if(is.double(wd) && 0.0 < wd && wd <= 1.0)
        wd <- options()$width * wd
    ml <- wd - nchar(bg) - nchar(ed) - 2 * nchar(sep)
    md <- rep(ch, times=ceiling(ml / nchar(ch)))
    md <- substr(paste0(md, collapse = ""), 1, ml)
    ln <- paste0(bg, md, ed, collapse = sep)
    cat(ln, "\n")
    invisible(NULL)
}

#' concatenate two strings by ":".
#'
#' wrapper of `paste(a, b, sep=":")`.
HLP$`%:%` <- function(a, b) paste(a, b, sep=":")

#' helper: two class confusion matrices
HLP$cfx <- function(ref, est) # confusion related
{
    ## confusion of case
    r <- addmargins(table(factor(0+ref, 0:1), factor(0+est, 0:1)))
    ## ret <- rbind(TPC=c(2, 2), TNF=c(1, 1), FPC=c(1, 2), FNC=c(2, 1))
    ## F1S = TP / (TP + 0.5 (FP + FN))
    ret <- c(
        TPF=r[2, 2] / r[2, 3], TNF=r[1, 1] / r[1, 3],
        FPF=r[1, 2] / r[1, 3], FNF=r[2, 1] / r[2, 3],
        PRC=r[2, 2] / (r[2, 2] + r[1, 2]), # precision=TPC/(TPC+FPC)
        RCL=r[2, 2] / (r[2, 2] + r[2, 1]), # recall   =TPC/(TPC+FNG)
        ACC = (r[1, 1] + r[2, 2]) / r[3, 3])

    ## F1 score for cases
    F1S <- r[2, 2] / (r[2, 2] + 0.5 * (r[1, 2] + r[2, 1]))

    ## confusion of ctrl
    s <- addmargins(table(factor(1-ref, 0:1), factor(1-est, 0:1)))
    F2S <- s[2, 2] / (s[2, 2] + 0.5 * (s[1, 2] + s[2, 1]))

    ## Macro Average Score
    FAS <- F1S * 0.5          + F2S * 0.5          # macro avg
    FWS <- F1S * mean(ref==1) + F2S * mean(ref==0) # macro avg weighted
    c(ret, F1S=F1S, FAS=FAS, FWS=FWS)
}

#' sigmoid function
#'
#' wrapper of [stats::binomial()].
HLP$sgm <- binomial()$linkinv

#' convert YYYY-MM-DD to yyyy-qN (year-quater).
HLP$y2q <- function(x)
{
    YMD <- "^(....)-(..)-(..)$"
    m <- sub(YMD, "\\2", x)
    y <- sub(YMD, "\\1", x)
    q <- (as.integer(m) - 1) %/% 3L + 1L
    sprintf("%s Q%d", y, q)
    ifelse(is.na(x), NA, sprintf("%s Q%d", y, q))
}

## tools for data clean up 
#' word check
#'
#' A initial quality summary of a data as words in bytes.
#' 
#' For each field, summerize the size  of words in bytes, missing rate, uniques,
#' followed by a preview of few non-missing words.
#'
#' Ill-formated data may have abnormally long  words for non-text fields such as
#' age, sex, and uid, due to  line shifting; highly missing and highly identical
#' fields may be considered to be dropped.
#'
#' Caution: the few values previewed here are not from the same records, because
#' the purpose of wor-check is to showcase some unique values in each field. Use
#' [head()] instead to preview a few compelte records.
#'
#' @param dat the table to check, typically a R [data.frame].
#' @param few how many elements to preview (def=4)
#' @param len maximum word character length to preview (def=10).
#' @param nas na-strings, (def={"", "NA", "N/A", "NULL"}).
HLP$wck <- function(dat, few=4, len=10, nas=NAS)
{
    rpt <- sapply(dat, \(x)
    {
        ## v <- iconv(a, "latin1", "UTF-8")
        v <- as.character(x)
        v <- v[!(is.na(v) | toupper(enc2native(v)) %in% nas)] # non-NA
        Encoding(v) <- "bytes"
        ## prop of missing
        pms <- c(PMS = round((length(x) - length(v)) / length(x), 4))
        ## word sizes
        ## QSZ <- c(LMN=0, LQ1=1, LMD=2, LQ3=3, LMX=4) / 4
        QSZ <- c(LMN=0, LMD=2, LMX=4) / 4
        wsz <- quantile(nchar(v, "bytes"), QSZ)
        names(wsz) <- names(QSZ)
        ## number of unique values
        nux <- c(NUX = length(unique(v)))
        ## preview a few unique values, truncated
        val <- substr(v, 1, len) |> unique() |> rep(length.out=few)
        names(val) <- sprintf("x%02d", 1:few)
        c(wsz, nux, pms, val)
    }) |> t()
    data.frame(rpt)
}

#' identifiable principle components
#'
#' 
#' wrap R's prcomp() to ensure maximum positive span, so the
#' principal components become identifiable.
#'
#' @param x data matrix to apply IPC.
#' @param mxp ensure positive maximum span? (def=1)
#' @param ... arguments to pass to prcomp().
HLP$ipc <- function(x, mxp=1, ...)
{
    pca <- prcomp(x, ...)
    pcs <- pca$x
    ldv <- pca$rotation
    sdv <- pca$sdev
    cnt <- pca$center
    if(mxp)
    {
        for(i in seq(ncol(pcs)))
            ldv[, i] <- ldv[, i] * sign(max(pcs[, i])^2 - min(pcs[, i])^2)
        pcs <- scale(x, cnt, FALSE) %*% ldv
    }
    structure(list(pcs=pcs, ldv=ldv, sdv=sdv, cnt=cnt), class=c("ipc", "list"))
}

#' identifiable class labels
#'
#' A wrapper of R's [reorder()] to ensure consistant order of factor levels with
#' principle components or any weighted coordinates.
#'
#' The order of  a level is dertermined by the  *aggregated proximity* of points
#' in that level to the origin of the space defined by `pcs`, with dimensions in
#' `pcs` weighted by `sdv`.
#'
#' By default, *proximity* to the origin is measured by Euclidian distance while
#' *aggregation* is done by passing [FUN=median()] to [reorder()].
#'
#' @param lbl {n} points labeled by factor levels.
#' @param pcs principle components or coordinates of label point.
#' @param sdv explanable standard deviations or dimension weights.
#' @param FUN function to calculated an aggregated proximity.
HLP$icl <- function(lbl, pcs, sdv, FUN=median)
{
    dst <- sqrt((pcs^2 %*% (1/sdv^2)))
    lbl <- reorder(lbl, dst, FUN)
    lbl
}
HLP$ilv <- HLP$icl

#' wrapper to calculate percentile instead of quantile
HLP$percentile <- function(x, probs=seq(0, 1, 0.01))
{
    quantile(x, probs, na.rm=TRUE)
}

if("TCOM:HLP" %in% search())
    detach("TCOM:HLP")
attach(HLP, name="TCOM:HLP", warn.conflicts = TRUE)
rm(HLP)
