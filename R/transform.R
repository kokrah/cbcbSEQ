#' Quantile transform counts.
#'
#' Perform normalize.quantiles on the count matrix, likely used with log2CPM().
#'
#' @param x raw counts
#' @param verbose print head of transformed matrix if TRUE (default=FALSE)
#' @return matrix of quantile normalized counts
#' @export
qNorm <- function(x, verbose=FALSE) {
    rn <- rownames(x)
    cn <- colnames(x)

    y <- normalize.quantiles(as.matrix(x), copy=TRUE)

    rownames(y) <- rn
    colnames(y) <- cn

    if (verbose) {
        print(head(y))
    }

    return(y)
}

#' Compute log2(quantile counts per mil reads) and library size for each sample.
#'
#' Performing a log2(cpm(quantile(counts))) of a count matrix is a quick way to get data ready for
#' visualization and analysis so that we may see trends more easily.  Unfortunately, that format is
#' inappropriate for any of the differential expression tools, ergo voomMod().
#'
#' @param qcounts quantile normalized counts
#' @param lib.size default is colsums(qcounts)
#' @return list containing log2(quantile counts per mil reads) and library sizes
#' @export
log2CPM <- function(qcounts, lib.size=NULL){
    if (is.null(lib.size)) {
        lib.size <- colSums(qcounts)
    }
    y <- t(log2(t(qcounts + 0.5)/(lib.size + 1) * 1e+06))
    retlist <- list(
        "y" = y,
        "lib.size" = lib.size)
    return(retlist)
}

#' Call modified voom on ComBatLoc adjusted data.
#'
#' This is a copy of limma's voom() with a couple changes to allow the input of log2(cpm()) data.
#'
#' @param x ComBatLoc adjusted data
#' @param design reparameterized design matrix
#' @param lib.size library sizes of quantile counts obtained from calling logPCM()
#' @return Elist object including a plot slot containing a ggplot2 version of voom()'s plot.
#' @export
voomMod <- function(x, design, lib.size) {
    ## minor modification to voom code:
    ## 1. Allows the function to take in data on log2CMP scale.
    ## Note to whomever reads this:  I have a separate version of this which detects the scale of
    ## the data as well as checks to see if it is cpm'd and reacts accordingly, this allowing it to
    ## be used with arbitrary data.  I didn't include those changes, but they are trivial if()
    ## statements which may be left as an exercise to the reader.
    out <- list()

    y <- as.matrix(x)
    fit <- lmFit(y, design)

    if (is.null(fit$Amean)) {
        fit$Amean <- rowMeans(y, na.rm = TRUE)
    }
    sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma)
    allzero <- rowSums(y) == 0
    if (any(allzero)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
    }
    l <- lowess(sx, sy, f = 0.5)
    ##if (plot) {
    ##    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
    ##         pch = 16, cex = 0.25)
    ##    title("voom: Mean-variance trend")
    ##    lines(l, col = "red")
    ##}
    f <- approxfun(l, rule = 2)
    mean_var_df <- data.frame(mean=sx, var=sy)
    mean_var_plot <- ggplot2::ggplot(mean_var_df, ggplot2::aes_string(x="mean", y="var")) +
        ggplot2::geom_point() +
        ggplot2::xlab("Log2(count size + 0.5)") +
        ggplot2::ylab("Square root of the standard deviation.") +
        ## stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE, show_guide=FALSE) +
        ggplot2::stat_density2d(geom="tile", ggplot2::aes_string(fill="..density..^0.25"),
                                contour=FALSE, show.legend=FALSE) +
        ggplot2::scale_fill_gradientn(colours=grDevices::colorRampPalette(c("white","black"))(256)) +
        ggplot2::geom_smooth(method="loess") +
        ggplot2::stat_function(fun=f, colour="red") +
        ggplot2::theme(legend.position="none")

    if (is.null(fit$rank)) {
        warning("Some samples cannot be balanced across the experimental design.")
    }

    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, j, drop = FALSE])
    } else {
        fitted.values <- fit$coef %*% t(fit$design)
    }

    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    out$E <- y
    out$weights <- w
    out$design <- design
    out$lib.size <- lib.size
    out[["plot"]] <- mean_var_plot

    new("EList", out)
}


#' Call ComBat function with minor adjustments
#'
#' It is worth noting that sva's ComBat has some of these modifications now.
#'
#' @param dat Genomic measure matrix (dimensions probe x sample)
#' @param batch Batch covariate (multiple batches allowed)
#' @param mod bio factor of interest
#' @param numCovs The column numbers of the variables in mod to be treated as continuous variables
#' @param noScale do not adjust scale for batch (use default(TRUE) for now)
#' @param prior.plots (Optional) TRUE gives prior plots. Used only when noScale is FALSE
#' @return list containing adjusted data(bayesdata), shift adjustments(gamma.star), and scale adjustments(delta.star)
#' @export
combatMod <- function(dat, batch, mod, numCovs=NULL, noScale=TRUE, prior.plots = FALSE) {

    ## Always use parametric (for now)
    par.prior <-TRUE

    mod <- cbind(mod, batch)
    check <- apply(mod, 2, function(x) all(x == 1))
    mod <- as.matrix(mod[, !check])
    colnames(mod)[ncol(mod)] <- "Batch"
    if (sum(check) > 0 & !is.null(numCovs)) {
        numCovs <- numCovs - 1
    }
    ## design <- survJamda::design.mat(mod, numCov = numCovs)
    design <- survJamda::design.mat(mod)
    batches <- survJamda::list.batch(mod)
    n.batch <- length(batches)
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    NAs <- any(is.na(dat))
    B.hat <- NULL
    ## This is taken from sva's github repository in helper.R
    Beta.NA <- function(y, X) {
        des <- X[!is.na(y), ]
        y1 <- y[!is.na(y) ]
        B <- solve(t(des)%*%des)%*%t(des)%*%y1
        B
    }
    var.pooled <- NULL
    if (NAs) {
        message(pasteo("Found ", sum(is.na(dat)), " missing data values."))
    }

    message("Standardizing Data across genes.")
    if (!NAs) {
        B.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(dat))
    } else {
        B.hat <- apply(dat, 1, Beta.NA, design)
    }
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]

    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, n.array)
    } else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var, na.rm = T)
    }

    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, n.array)))

  #----------------------------------------------------------------------------------#
  # --------------------- Add no scale code here !!--------------------------------- #
  #----------------------------------------------------------------------------------#

    if (isTRUE(noScale)) {

        m.data <- dat - stand.mean
        mse <- ((dat - t(design %*% B.hat))^2) %*% rep(1/(n.array-ncol(design)), n.array)
        hld <- NULL
        bayesdata <- dat

    for (k in 1:n.batch) {
        message(paste0("Fitting 'shrunk' batch ", k, " effects."))
        sel <- batches[[k]]
        gammaMLE <- rowMeans(m.data[, sel])
        mprior <- mean(gammaMLE, na.rm=TRUE)
        vprior <- var(gammaMLE, na.rm=TRUE)
        prop <- vprior / (mse/(length(sel)) + vprior)
        gammaPost <- prop*gammaMLE + (1 - prop)*mprior

        for(i in sel) {
            bayesdata[,i] <- bayesdata[,i] - gammaPost
        }

        stats <- data.frame(gammaPost=gammaPost, gammaMLE=gammaMLE, prop=prop)
        hld[[paste("Batch", k, sep=".")]] <- list(
            "stats" = stats,
            "indices" = sel,
            "mprior" = mprior,
            "vprior" = vprior)
    } ## End for k in 1:n.batch

        message("Adjusting data for batch effects.")
        return(bayesdata)

    } else {
        ## Continue with orignal method
        message("Fitting L/S model and finding priors.")
        batch.design <- design[, 1:n.batch]
        if (!NAs) {
            gamma.hat <- solve(t(batch.design) %*% batch.design) %*%
                t(batch.design) %*% t(as.matrix(s.data))
        } else {
            gamma.hat = apply(s.data, 1, Beta.NA, batch.design)
        }
        delta.hat <- NULL
        for (i in batches) {
            delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, na.rm = T))
        }
        gamma.bar <- apply(gamma.hat, 1, mean)
        t2 <- apply(gamma.hat, 1, var)
        a.prior <- apply(delta.hat, 1, sva:::aprior)
        b.prior <- apply(delta.hat, 1, sva:::bprior)
        if (prior.plots & par.prior) {
            par(mfrow = c(2, 2))
            tmp <- density(gamma.hat[1, ])
            plot(tmp, type = "l", main = "Density Plot")
            xx <- seq(min(tmp$x), max(tmp$x), length = 100)
            lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
            qqnorm(gamma.hat[1, ])
            qqline(gamma.hat[1, ], col = 2)
            tmp <- density(delta.hat[1, ])
            invgam <- 1/rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
            tmp1 <- density(invgam)
            plot(tmp, typ = "l", main = "Density Plot", ylim = c(0, max(tmp$y, tmp1$y)))
            lines(tmp1, col = 2)
            qqplot(delta.hat[1, ], invgam, xlab = "Sample Quantiles",
                   ylab = "Theoretical Quantiles")
            lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
            title("Q-Q Plot")
        }
        gamma.star <- delta.star <- NULL
        if (par.prior) {
            message("Finding parametric adjustments.")
            for (i in 1:n.batch) {
                temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i,], delta.hat[i, ],
                                     gamma.bar[i], t2[i],
                                     a.prior[i], b.prior[i])
                gamma.star <- rbind(gamma.star, temp[1, ])
                delta.star <- rbind(delta.star, temp[2, ])
            }
        } else {
            message("Finding nonparametric adjustments.")
            for (i in 1:n.batch) {
                temp <- sva:::int.prior(as.matrix(s.data[, batches[[i]]]),
                                        gamma.hat[i, ], delta.hat[i, ])
                gamma.star <- rbind(gamma.star, temp[1, ])
                delta.star <- rbind(delta.star, temp[2, ])
            }
        }
        message("Adjusting the Data.")
        bayesdata <- s.data
        j <- 1

    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,] %*% gamma.star))/
            (sqrt(delta.star[j, ]) %*% t(rep(1, n.batches[j])))
        j <- j + 1
    }

        bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, n.array)))) + stand.mean
        return(bayesdata)
    }
}
