
## FIXME: plotting is broken
##' \description{
##'   Modified plotting function to better display labels and print names/colors for 
##'   the constraints. Input is an lpSovleAPI externalPtr object model
##'   Note: only works for simple linear models with two decision variables.
##' }
##' @title Modified plotting function for simple linear models with two desicion
##'     variables
##' @param x lpSolveAPI externalPtr solved model, as passed to `plot`
##' @param y 
##' @param label.width granularity of labels on axes
##' @param ... ignored
plt <- function (x, y, label.width=10, ...)  {
    m <- dim(x)[1]
    n <- dim(x)[2]
    pal <- ifelse(require(RColorBrewer) && !is.null(dimnames(mod)[[1]]),
                  brewer.pal(length(dimnames(mod)[[1]]), name="Dark2"),
                  1:length(dimnames(mod)[[1]]))

    if (n != 2 || any(get.bounds(x)$lower != 0) || 
        any(get.bounds(x)$upper != Inf) || any(get.type(x) != "real")) 
        stop("cannot plot linear program - see help(plot.lpExtPtr)")
    X <- matrix(0, m, n)
    col <- get.column(x, 1)
    if (col$nzrow[1] == 0) {
        col$nzrow <- col$nzrow[-1]
        col$column <- col$column[-1]
    }
    X[col$nzrow, 1] <- col$column
    col <- get.column(x, 2)
    if (col$nzrow[1] == 0) {
        col$nzrow <- col$nzrow[-1]
        col$column <- col$column[-1]
    }
    X[col$nzrow, 2] <- col$column
    X <- rbind(X, diag(2))
    b <- c(get.rhs(x), rep(0, 2))
    pts <- matrix(0, 2, 0)
    for (i in 1:(m + 1)) {
        for (j in (i + 1):(m + 2)) {
            A <- X[c(i, j), ]
            if (abs(det(A)) > 0.000000001) {
                pts <- cbind(pts, solve(A, b[c(i, j)]))
            }
        }
    }
    X <- X[1:m, ]
    b <- b[1:m]
    feasible <- double(0)
    for (j in 1:dim(pts)[2]) if (all(X %*% pts[, j, drop = FALSE] <= 
                                     b) && all(pts[, j] >= 0)) 
                                 feasible <- c(feasible, j)
    pts <- t(pts[, feasible])
    pts <- pts[chull(pts), ]
    limits <- max(pts) * c(-0.1, 1.4)
    box <- max(pts) * c(-0.05, 1.2)
    old.par <- par(mar = c(3, 2, 4, 2) + 0.1, pty = "s")
    on.exit(par(old.par))
    plot(NA, NA, type = "n", axes = FALSE, xlim = limits, ylim = limits, 
         xlab = "", ylab = "", pty = "s", main = "The Feasible Set")
    polygon(pts, col = "gray")
    arrows(0, limits[1], 0, limits[2], lwd = 1.5, length = 0.125)
    arrows(limits[1], 0, limits[2], 0, lwd = 1.5, length = 0.125)
    for (i in seq(0, floor(limits[2]), by=label.width)) {
        segments(i, 0, i, 0.0625)
        text(i, 0, as.character(i), adj = c(0.75, 1.5), srt=45)
        segments(0, i, 0.0625, i)
        text(0, i, as.character(i), adj = 1, srt=45)
    }
    decision.vars <- dimnames(x)[[2]]
    text(limits[2], limits[1], paste(decision.vars[1], " ", sep = ""), 
         adj = 1)
    text(0, limits[2], paste("   ", decision.vars[2], sep = ""), 
         adj = 0)
    for (i in 1:m) {
        if (abs(X[i, 1]) < 0.000000001) 
            segments(box[1], b[i]/X[i, 2], box[2], b[i]/X[i, 
                                                          2])
        else if (abs(X[i, 2]) < 0.000000001) 
            segments(b[i]/X[i, 1], box[1], b[i]/X[i, 1], box[2])
        else {
            y.intercept <- b[i]/X[i, 2]
            if (y.intercept > box[2]) {
                x0 <- (b[i] - X[i, 2] * box[2])/X[i, 1]
                y0 <- box[2]
            }
            else if (y.intercept < box[1]) {
                x0 <- (b[i] - X[i, 2] * box[1])/X[i, 1]
                y0 <- box[1]
            }
            else {
                x0 <- box[1]
                y0 <- (b[i] - X[i, 1] * box[1])/X[i, 2]
            }
            y.intercept <- (b[i] - X[i, 1] * box[2])/X[i, 2]
            if (y.intercept > box[2]) {
                x1 <- (b[i] - X[i, 2] * box[2])/X[i, 1]
                y1 <- box[2]
            }
            else if (y.intercept < box[1]) {
                x1 <- (b[i] - X[i, 2] * box[1])/X[i, 1]
                y1 <- box[1]
            }
            else {
                x1 <- box[2]
                y1 <- (b[i] - X[i, 1] * box[2])/X[i, 2]
            }
            segments(x0, y0, x1, y1, col=pal[i])
            if (!is.null(dimnames(mod)[[1]])) {
                text(x0+label.width, y0, 
                     labels=dimnames(mod)[[1]][i], srt=45, cex=1.1, col=pal[i])
            }
        }
    }
    invisible()
}

## Output LP results and sensitivity report similar to that given 
## by Excel
##' @title Create an LP results / sensitivity report similar to Excel's tables
##' @param mod A solved LP program, an lpExtPtr, created using the lpSolveAPI
##' @param sensitivity return sensitivity results as well
##' @return list of results and sensitivity analyses
##' @author Noah Peart
lp.report <- function(mod, sensitivity=TRUE) {
    stopifnot(require(lpSolveAPI))
    opts <- options(scipen=999)
    on.exit(options(opts))
    m <- dim(mod)[1]           #constraints
    n <- dim(mod)[2]           #decision vars
    
    primals <- get.primal.solution(mod)
    nms <- c("Objective", unlist(dimnames(mod)))
    mi <- 2:(m + 1)            #constr. inds
    ni <- (m + 2):(n + m + 1)  #desc. inds
    
    ## slack
    slack <- abs(round(primals[mi] - get.constr.value(mod), 4)) 
   
    ## Results
    results <- list(
        objective = setNames(primals[1], nms[1]),
        decision = data.frame(
            Name=nms[ni],
            Original.Value=get.bounds(mod)$lower,
            Final.Value=primals[ni]
        ),
        constraints = data.frame(
            Name=nms[mi],
            Value=primals[mi],
            Status=ifelse(slack - 1e-13 < 0, "Binding", "Not Binding"),
            Slack=slack
        )
    )
    
    if (!sensitivity) {
        return(list(results=results, sensitivity=NULL))
    }
    
    ## Sensitivity Analysis
    dual <- get.dual.solution(mod)
    dual.sens <- get.sensitivity.rhs(mod)
    obj.sens <- get.sensitivity.objex(mod)
    coefs <- sapply(seq.int(dim(mod)[2]), function(x) get.column(mod, x)$column[[1]])

    sensitivity <- list(
        decision = data.frame(
            Name=nms[ni],
            Final.Value=primals[ni],
            Reduced.Cost=dual[ni],
            Obj.Coef=coefs,
            Allowed.Increase=coefs - obj.sens$objfrom,
            Allowed.Decrease=obj.sens$objtill - coefs
        ),
        constraints = data.frame(
            Name=nms[mi],
            Final.Value=primals[mi],
            Shadow.Price=dual[mi],
            Constraint.RHS=get.constr.value(mod, "rhs"),
            Allowed.Increase=dual.sens$dualstill[mi-1] - primals[mi],
            Allowed.Decrease=primals[mi] - dual.sens$dualsfrom[mi-1]
        )
    )

    list(results=results, sensitivity=sensitivity)
}

## Print summary of results/sensitivity analysis of LP model
lp.print <- function(results) {
    invisible(lapply(seq_along(results), function(i) {
        cat(names(results[i]))
        if (is.data.frame(results[[i]])) {
            cat("\n---------------")
            print(knitr::kable(format(results[[i]], nsmall=3), 
                               type="html"))
            cat("\n")
        } else {
            cat(": ", results[[i]], "\n\n")
        }
    }))
}

## FIXME: this should be a dispatch from a call to summary
lp.summary <- function(mod, sensitivity=TRUE) {
    stopifnot(require(knitr))
    res <- lp.report(mod, sensitivity)
    
    ## Objective, decision, and constraint results
    lp.print(res$results)
    
    ## Sensitivity Results
    if (sensitivity) {
        cat("\n\t------------------------\n",
            "\t|  Sensitivity Report  |\n",
            "\t------------------------\n")
        lp.print(res$sensitivity)
    }

    invisible(res)
}
