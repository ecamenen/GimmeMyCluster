
################################
#            MAIN
################################

# Global variables settings
TEXT <- TRUE # print values on graph (for optimum partition and heatmap)
NB_ROW_MAX <- 200 # max row to have pdf, otherwise, some plots are in png
MAX_CHAR_LEN <- 25 # maximum length of individual s names

printProgress <- function(verbose = TRUE, val) {
    if (isTRUE(verbose)) {
        cat(paste("\n[", format(Sys.time(), "%X"), "] ", val, "in progress..."), sep = "")
    }
}

getTimeElapsed <- function(start_time) {
    time <- as.numeric(as.difftime(Sys.time() - start_time), units = "secs")
    secs <- time %% 60
    time <- (time - secs) / 60
    mins <- time %% 60
    hours <- time / 60
    time <- paste(mins, "min ", round(secs), "s\n", sep = "")

    if (hours >= 1) {
        time <- paste(round(hours), "h ", time, sep = "")
    }

    cat(paste("\nTime to run the process : ", time, sep = ""))
}

################################
#          Parsing
################################


# rename row and avoid doublets errors
preProcessData <- function(d, header = FALSE, verbose = FALSE) {
    if (ncol(d) == 1) {
        stop(paste("Check for the separator (by default, tabulation)."),
            call. = FALSE
        )
    }

    # Discard line containing the same information on all columns from analysis
    REMOVE_DOUBLETS <- (nrow(d) > NB_ROW_MAX)
    d <- renameRowname(d)
    # remove columns containing characters
    # if(nrow(d) > NB_ROW_MAX) verbose = TRUE
    # get only columns with numeric values
    d <- d[, unlist(sapply(1:ncol(d), function(i) is.numeric(d[, i])))]

    if (isTRUE(REMOVE_DOUBLETS)) {
        printProgress(verbose, "Loading data")
        d <- discardRowCondDoublets(d)
    }

    if (isSymmetric(as.matrix(d)) & !header) {
        colnames(d) <- rownames(d)
    }

    return(d)
}


# avoid doublets in row names
# r: row names vector
renameRownameDoublets <- function(names.row) {
    j <- 1

    for (i in 2:length(names.row)) {
        if (names.row[i] == names.row[i - 1]) {
            j <- j + 1
            names.row[i] <- paste(names.row[i], ".", j, sep = "")
        } else {
            j <- 1
        }
    }
    return(names.row)
}

# rename row and avoid doublets errors
renameRowname <- function(d) {
    # names.row = as.character(d[,1])
    names.row <- rownames(d)
    d <- d[, -1]
    names.row <- renameRownameDoublets(names.row)

    tryCatch(
        {
            substr(names.row, 1, MAX_CHAR_LEN) -> rownames(d)
            return(d)
        },
        warning = function(w) {
            names.row <- renameRownameDoublets(substr(names.row, 1, MAX_CHAR_LEN))
            names.row -> rownames(d)
            return(d)
        },
        error = function(e) {
            return(d)
        }
    )
}

# Discard row from a reaction dataset that have the same conditions in each columns
# x: dataframe
discardRowCondDoublets <- function(x) {
    row_doublets <- list()
    j <- 0

    for (i in 1:nrow(x)) {
        # uniq remove doublets in a vector, so return 1 only if there is only 1
        if ((length(unique(as.integer(x[i, ]))) == 1)) {
            # print(row.names(x[i,]))
            j <- j + 1
            row_doublets[[j]] <- i
        }
    }

    if (length(row_doublets) != 0) {
        removed_reacs <- row.names(x[unlist(row_doublets), ])
        removed_conds <- x[unlist(row_doublets), 1]
        removed <- cbind(removed_reacs, removed_conds)
        colnames(removed) <- c("condition", "")
        writeTsv(removed, verbose = FALSE)
    }

    if (length(row_doublets) > 0) {
        return(x[-unlist(row_doublets), ])
    } else {
        return(x)
    }
}

# rename row and avoid doublets errors
renameRowname <- function(d) {
    names.row <- as.character(d[, 1])
    d <- d[, -1]
    names.row <- renameRownameDoublets(names.row)
    tryCatch(
        {
            substr(names.row, 1, 25) -> rownames(d)
            return(d)
        },
        warning = function(w) {
            names.row <- renameRownameDoublets(substr(names.row, 1, 25))
            names.row -> rownames(d)
            return(d)
        },
        error = function(e) {
            return(d)
        }
    )
}

# Discard row from a reaction dataset that have the same conditions in each columns
# x: dataframe
discardRowCondDoublets <- function(x) {
    row_doublets <- list()
    j <- 0
    for (i in 1:nrow(x)) {
        # uniq remove doublets in a vector, so return 1 only if there is only 1
        if ((length(unique(as.integer(x[i, ]))) == 1)) {
            # print(row.names(x[i,]))
            j <- j + 1
            row_doublets[[j]] <- i
        }
    }
    if (length(row_doublets) != 0) {
        removed_reacs <- row.names(x[unlist(row_doublets), ])
        removed_conds <- x[unlist(row_doublets), 1]
        removed <- cbind(removed_reacs, removed_conds)
        colnames(removed) <- c("condition", "")
        writeTsv(removed, verbose = FALSE)
    }
    if (length(row_doublets) > 0) {
        return(x[-unlist(row_doublets), ])
    } else {
        return(x)
    }
}

# Inputs: x : a matrix
# filename of the saved file
# Prints the matrix, save the matrix
writeTsv <- function(x, f = NULL, cl = FALSE, verbose = TRUE) {
    # print on stdout
    if (isTRUE(verbose)) {
        cat(paste("\n", gsub("_", " ", toupper(x)), ":\n", sep = ""))
    }
    # disabling warning
    options(warn = -1)
    # get variable
    tab <- x
    if (!isTRUE(cl)) {
        output <- as.matrix(rbind(c("", colnames(tab)), cbind(rownames(tab), tab)))
    } else {
        output <- tab
    }
    # discard empty rows
    output <- output[rowSums(is.na(output)) != ncol(output), ]
    # TODOD:
    # output = output[,colSums(is.na(output)) != nrow(output)]
    output[is.na(output)] <- ""
    colnames(output) <- rep("", ncol(output))
    rownames(output) <- rep("", nrow(output))
    if (isTRUE(verbose)) {
        if (!isTRUE(cl)) {
            printed <- round(apply(output[-1, -1], 2, as.numeric), 2)
            rownames(printed) <- rownames(tab)
            colnames(printed) <- colnames(tab)
        } else {
            printed <- output
        }
        print(printed, quote = FALSE)
    }
    if (is.null(f)) {
        f <- paste(deparse(substitute(x)), ".tsv", sep = "")
    }
    write(t(output), f, ncolumns = ncol(output), sep = "\t")
    options(warn = 0)
}

################################
#          Graphic
################################

# Usage: colPers(x), x a number of colours in output
# Gradient of color
colPers <- colorRampPalette(c(
    rgb(0.6, 0.1, 0.5, 1),
    rgb(1, 0, 0, 1),
    rgb(0.9, 0.6, 0, 1),
    rgb(0.1, 0.6, 0.3, 1),
    rgb(0.1, 0.6, 0.5, 1),
    rgb(0, 0, 1, 1)
),
alpha = TRUE
)


setGraphic <- function() {
    setGraphicBasic()
    par(mar = c(5.1, 5.1, 5.1, 2.1))
}

setGraphicBasic <- function() {
    par(
        cex.lab = 1.5,
        font.lab = 3,
        font.axis = 3,
        cex.axis = 0.8,
        cex.main = 2,
        cex = 1,
        lwd = 3
    )
}

plotAxis <- function(side, min, max, interval = 1, lwd = 3) {
    axis(side, seq(min, max, interval), lwd = lwd)
}

#' @export
plotBestClustering <- function(
    sub_title,
    values,
    values_type,
    n_cluster = 2,
    max_cluster = 6,
    interval = 1,
    min_x = 2,
    best = NULL,
    val2 = NULL,
    verbose = FALSE) {

    plotAxis(1, 2, max_cluster)

    if (interval >= 1) {
        axisSeq <- round(values)
    } else {
        axisSeq <- c(0, max(values) + 0.1)
    }

    # case of plotting gap statistics
    if (min_x < 2) {
        best_y <- values[n_cluster]
    } # case of fusion levels
    else if (!is.null(val2)) {
        best_y <- values[n_cluster - 1]
    } else {
        best_y <- max(values)
    }

    # for non-elbow plots
    if (!is.null(val2)) {
        best <- round(max(val2), 2)
    } else if (is.null(best)) {
        best <- round(max(values), 2)
    }

    plotAxis(2, min(axisSeq), max(axisSeq), interval)
    title(
        main = "Optimal number of clusters",
        line = 2,
        cex.main = 2
    )
    mtext(
        text = sub_title,
        font = 3,
        cex = 1.2,
        line = 0.5
    )
    abline(
        v = n_cluster,
        col = "red",
        lty = 2,
        lwd = 2
    )
    points(
        n_cluster,
        best_y,
        pch = 19,
        col = "red",
        cex = 2
    )

    if (!is.null(val2)) {
        t_values <- val2
    } else {
        t_values <- values
    }

    if (isTRUE(TEXT)) {
        text(
            y = values,
            x = min_x:max_cluster,
            labels = round(t_values, 2),
            cex = 1.2,
            pos = 4,
            col = "red"
        )
    }
    if (isTRUE(verbose)) {
        cat(
            "Optimal number of clusters k = ",
            n_cluster,
            "\n",
            "With a",
            values_type,
            " of ",
            best,
            "\n",
            sep = ""
        )
    }
}

# f: filename
savePdf <- function(f) {
    pdf(f)
    setGraphic()
}

################################
#          Statistics
################################

# Get the normalized distance between each points and the center
# Outputs:
# for each column, the mean=0 and the variance is the same
#' @export
scalecenter <- function(x) {
    # output scale function: for each column, mean=0, sd=1
    return(scale(x) * sqrt(nrow(x) / (nrow(x) - 1)))
    # ponderation for sampling index (var use n-1)
    # without this constante, for advanced outputs, total (max_cluster=nrow(data)) will be different from 1
}

# data: dataframe
# method: distance method
#' @export
getDistance <- function(x, method = 1) {
    dists <- c("euclidean", "manhattan", 1, 2, 5, 7)

    if (method < 3) {
        dist(x, method = dists[method])
    } else {
        dist.binary(x, method = as.integer(dists[method]))
    }
}

# x distance object
checkEuclidean <- function(x) {
    if (attributes(x)$method != "euclidean") {
        stop("Distance should be euclidean with this classification method.", call. = FALSE)
    }
}

# x distance object
isSymmetric <- function(x) {
    if (nrow(x) == ncol(x)) {
        isReflexivity <- unique(x[cbind(1:nrow(x), 1:nrow(x))] == 0)

        if (length(isReflexivity) == 1 & isTRUE(isReflexivity)) {
            isCommutativity <- unique(x[lower.tri(x)] == t(x)[lower.tri(x)])

            if (length(isCommutativity) == 1 & isTRUE(isCommutativity)) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

################################
#          Clustering
################################

# Inputs:
# method: number of method of classification
# data: data
# x: distance matrix
# Ouput: Hierarchical classification
getCAH <- function(x, method = 3) {
    if (method > 2) {
        if (method == 8 | method == 9) {
            checkEuclidean(x)
        }

        # cah: classification hierarchic ascending
        cah <- hclust(x, method = getClassifType(method))
        # automaticly ordering by clusters
        return(reorder.hclust(cah, x))
    }
}

# Selects best algo based on cophenetic calculation
# x: data
# d: distance matrix
selectBestCAH <- function(x, dist, verbose = FALSE) {
    temp <- 0
    for (i in 3:9) {
        cah <- getCAH(x, method = i)
        res <- cor(dist, cophenetic(cah))
        if (isTRUE(verbose)) {
            cat(paste(getClassifType(i), ":", round(res, 3), "\n"))
        }
        if (res > temp) {
            temp <- res
            t <- i
        }
    }
    cat(paste("Selected:", getClassifType(t), "\n"))
    return(t)
}

# Inputs:
# x: number of method of classification
getClassifType <- function(x) {
    methods <- c(
        "kmedoids",
        "kmeans",
        "ward.D2",
        "complete",
        "single",
        "average",
        "mcquitty",
        "median",
        "centroid"
    )
    methods[x]
}

# Agglomerative coefficient
# x: clustering object
getCoefAggl <- function(x) {
    coef.hclust(x)
}

# Inputs:
# method: number of method of classification
# data: data
# dist: distance
# n_cluster: number of clustering
# Ouput: Non-hierarchical classification
getCNH <- function(data = NULL, dist = NULL, n_cluster = 2, method = 1, centers = NULL, nrun = 100) {
    if (method == 1) {
        if (!is.null(centers)) {
            return(pam(dist, n_cluster, diss = TRUE, medoids = centers))
        }
        return(pam(dist, n_cluster, diss = TRUE))
    } else if (method == 2) {
        if (!is.null(centers)) {
            n_cluster <- centers
        }
        return(kmeans(data, centers = n_cluster, nstart = nrun))
    }
}

#' @export
getClassif <- function(data = NULL, dist = NULL, method = 3, max_cluster = 6, centers = NULL, nrun = 100) {
    if (method > 2) {
        getCAH(dist, method = method)
    } else {
        list_cnh <- list("method" = getClassifType(method))
        for (k in 2:(max_cluster + 1)) {
            list_cnh[[k]] <- getCNH(data = data, dist = dist, method = method, n_cluster = k, centers = centers[[k - 1]], nrun = nrun)
        }
        return(list_cnh)
    }
}

# Inputs:
# n_cluster: number of clusters
# x: hierarchical classification
# d: data
# Output: partitionning contening n_cluster clusters
#' @export
getClusters <- function(x, n_cluster = 2) {
    if (x$method == "kmedoids") {
        x[[n_cluster]]$clustering
    } else if (x$method == "kmeans") {
        x[[n_cluster]]$cluster
    } else {
        cutree(x, n_cluster)
    }
}

#' @export
getClusterPerPart <- function(x, max_cluster = 6) {
    cl <- list()
    for (k in 2:max_cluster) {
        cl[[k - 1]] <- getClusters(x, n_cluster = k)
    }
    return(cl)
}

# Input:
# x: clusters
colorClusters <- function(x) {
    NB_CLUSTERS <- length(levels(as.factor(x)))
    for (i in 1:NB_CLUSTERS) {
        x[x == i] <- colPers(NB_CLUSTERS)[i]
    }
    return(x)
}

# Inputs:
# cl: clusters
# f : filename
# r: ordered alphabetically
writeClusters <- function(f, sil_k, pca, verbose = FALSE) {
    cluster <- cbind(sil_k[, 1], sil_k[, 3], pca$li[attr(sil_k, "iOrd"), c(1, 2)])
    colnames(cluster) <- c("Cluster", "Silhouette", "Axis1", "Axis2")
    writeTsv(cluster, f, cl = FALSE, verbose = verbose)
}

############################################################
#          Cophenetic (dendrogram distance matrix)
############################################################

# Distance matrix between each leaf of the dendogramm
# Inputs:
# dist : Distance matrix
# x : hierarchical classification
#' @export
plotCohenetic <- function(x, dist, is_png = FALSE, verbose = FALSE) {
    coph_matrix <- cophenetic(x)
    cor_coph <- cor(dist, coph_matrix)
    if (isTRUE(verbose)) {
        cat(
            paste(
                "\nCOPHENETIC:\nExplained variance (%):",
                round(cor_coph^2, 3),
                "\nCorrelation with the data:",
                round(cor_coph, 3),
                "\n"
            )
        )
    }

    # if(is_png) {
    #   #png(paste(opt$output8, ".png", sep=""), DIM_PNG/2, DIM_PNG/2)
    #   par(cex.lab=1.5*2, font.lab=3, font.axis=3, cex.axis=0.8*2, cex.main=2*2, cex=1, lwd=3*2)
    #   par(mar=c(5.1,5.1,5.1,2.1)+7)
    #   lwd=3*2
    #   line.lab = 5
    # }else{
    setGraphic()
    # savePdf(paste(opt$output8, ".pdf", sep=""))
    lwd <- 3
    line.lab <- 3
    # }

    plot(
        dist,
        coph_matrix,
        pch = 19,
        col = alpha("red", 0.2),
        axes = FALSE,
        xlim = c(0, max(dist)),
        xlab = "",
        ylab = "",
        ylim = c(0, max(coph_matrix)),
        asp = 1,
        main = paste("Cophenetic correlation: ", round(cor_coph, 3))
    )
    title(
        xlab = "Distance between metabolites",
        ylab = "Cophenetic distance",
        line = line.lab
    )
    plotAxis(2, 0, max(coph_matrix), lwd = lwd)
    plotAxis(1, 0, max(dist), lwd = lwd)
    abline(0, 1, col = "grey", lty = 2, lwd = lwd)
    # suprLog = dev.off()
}

##############################################
#          Inertia
##############################################

# Relative inter-group inertia for each partitionning
# Inputs:
# max_cluster: maximum number of clusters
# x: dataframe
# cl: list of clusters per partition
#' @export
getRelativeBetweenPerPart <- function(x, cl, max_cluster = 6) {
    x <- as.matrix(x)
    between <- rep(0, max_cluster - 1)
    # total sum of square
    TSS <- sum(scale(x, scale = FALSE)^2)
    for (i in 2:max_cluster) {
        cl_k <- as.factor(cl[[i - 1]])
        # apply(x, 2, mean) : centroids for each column
        # as.vector(table(cl) : size of each clusters
        # t : vector rotation for arithmetic with other row or column vectors
        between[i - 1] <- sum(t((t(getClusterCentroids(x, cl_k)) - apply(x, 2, mean))^2) * as.vector(table(cl_k))) / TSS
    }
    return(100 * between)
}

# x: raw dataset
# cl: list of clusters per partition
# tapply(data[,i], Cla, mean) :
# centroids of each clusters for a column i
# sapply(1:ncol(data), function(i) tapply(data[,i], Cla, mean)) :
# centroids of each clusters for each column
getClusterCentroids <- function(x, cl) {
    sapply(1:ncol(x), function(i) tapply(x[, i], cl, mean))
}

# Difference between each case of a vector
# x: between-inertia vector of integer
#' @export
getBetweenDifferences <- function(x) {
    # apply produce a list, unlist convert in vector
    diff <- unlist(sapply(1:length(x), function(i) x[i] - x[i - 1]))
    return(as.vector(cbind(x[1], t(diff))))
    #-n-1 to remove the last NA value (pairwise comparison)
    # x[1] to get the difference with 1 cluster
}

# x: data
# i_cluster: cluster id
getWithin <- function(x, cl, i_cluster = 2) {
    cluster_size <- length(cl[cl == i_cluster]) # number of individuals in the cluster
    scaled <- scalecenter(x)
    return(cluster_size * sum(getClusterCentroids(scaled, cl = cl)[i_cluster, ]^2) / nrow(x))
}

# cl: list of clusters per partition
#' @export
getRelativeWithinPerCluster <- function(x, cl) {
    max_cluster <- length(cl)
    within <- matrix(NA, max_cluster, max_cluster + 1)
    rownames(within) <- seq(2, max_cluster + 1)
    colnames(within) <- paste("G", seq(1, max_cluster + 1), sep = "")
    for (k in 2:(max_cluster + 1)) {
        cl_k <- cl[[k - 1]]
        for (i in 1:length(table(cl_k))) {
            within[k - 1, i] <- getWithin(x, cl = cl_k, i_cluster = i)
        }
        within[k - 1, ] <- within[k - 1, ] / sum(as.numeric(na.omit(within[k - 1, ])))
    }
    return(within)
}

# Between inertia differences between a partionning and the previous
# x: numerical vector of between differences
plotBetweenDiff <- function(x, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nBETWEEN DIFFERENCES:\n")
    }
    n_cluster <- which.max(x) + 1
    max_cluster <- length(x)
    # savePdf("between_differences.pdf")
    setGraphic()
    plot(
        2:(max_cluster + 1),
        x,
        type = "b",
        ylim = c(round(min(x)) - 1, round(max(x)) + 1),
        xlim = c(2, (max_cluster + 2)),
        xlab = "Nb. of clusters",
        ylab = "Between-cluster variation (%)",
        col = "grey",
        axes = FALSE
    )
    plotBestClustering(
        "Largest between differences method",
        x,
        " variation with the previous partitionning (%)",
        n_cluster,
        verbose = verbose,
        max_cluster = max_cluster
    )
    # suprLog = dev.off()
}

# x: clustering object
#' @export
plotFusionLevels <- function(x, max_cluster = 6, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nFUSION LEVELS:\n")
    }
    fusion <- rev(x$height)
    diff <- unlist(sapply(1:max_cluster, function(i) fusion[i - 1] - fusion[i]))
    fusion <- fusion[1:(max_cluster - 1)]
    n_cluster <- which.max(diff) + 1
    setGraphic()
    # savePdf(paste(opt$output9, ".pdf", sep=""))
    plot(
        2:max_cluster,
        fusion,
        type = "b",
        ylim = c(round(min(fusion)) - 1, round(max(fusion)) + 1),
        xlim = c(2, max_cluster + 1),
        xlab = "Nb. of clusters",
        ylab = "Cophenetic distance",
        col = "grey",
        axes = FALSE
    )
    plotBestClustering(
        "Fusion level method",
        fusion,
        " gain with the previous fusion level",
        n_cluster,
        val2 = diff,
        verbose = verbose,
        max_cluster = max_cluster
    )
    # suprLog = dev.off()
}

# x: vector of between inertia for k partitions
#' @export
plotElbow <- function(x, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nELBOW:\n")
    }
    max_cluster <- length(x) + 1
    within <- c(100, 100 - x)
    ratio <- within[1:(max_cluster - 1)] / within[2:max_cluster]
    best <- round(min(ratio), 2)
    n_cluster <- which.min(ratio)
    setGraphic()
    # savePdf("elbow.pdf")
    plot(
        1:max_cluster,
        within,
        type = "b",
        ylim = c(min(within) - 1, 101),
        xlim = c(1, max_cluster + 1),
        xlab = "Nb. of clusters",
        ylab = "Relative within inertia (%)",
        col = "grey",
        axes = FALSE
    )
    plotBestClustering(
        "Elbow method",
        within,
        " Wk/Wk+1 ratio ",
        n_cluster,
        5,
        1,
        best,
        max_cluster = max_cluster,
        verbose = verbose
    )
    # suprLog = dev.off()
}

################################
#          Silhouette
################################

# Ouput: an ordered silhouette object
#' @export
getSilhouette <- function(dist, cl_k) {
    s <- sortSilhouette(silhouette(cl_k, dist))
    rownames(s) <- labels(dist)[attr(s, "iOrd")]
    return(s)
}

#' @export
getSilhouettePerPart <- function(dist, cl) {
    list_sil <- list()
    for (k in 2:(length(cl) + 1)) {
        list_sil[[k - 1]] <- getSilhouette(dist, cl[[k - 1]])
    }
    return(list_sil)
}

# x: list of silhouettes objects per partition
#' @export
getMeanSilhouettePerPart <- function(x) {
    unlist(sapply(1:length(x), function(i) summary(x[[i]])$avg.width))
}

# Plot the best average silhouette width for all clustering possible
# x: vector of silhouette average width
#' @export
plotSilhouettePerPart <- function(x, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nSILHOUETTE:\n")
    }
    setGraphic()
    # savePdf(opt$output1)
    n_cluster <- which.max(x) + 1
    max_cluster <- length(x)
    plot(
        2:(max_cluster + 1),
        x,
        type = "b",
        xlim = c(2, max_cluster + 2),
        ylim = c(0, max(x) + 0.1),
        col = "grey",
        xlab = "Nb. of clusters",
        ylab = "Average silhouette width",
        axes = FALSE
    )
    plotBestClustering(
        "Silhouette method",
        x,
        "max_cluster average width",
        n_cluster,
        0.1,
        verbose = verbose,
        max_cluster = (max_cluster + 1)
    )
    # suprLog = dev.off()
    # return (n_cluster)
}

# x: a silhouette object
#' @export
plotSilhouette <- function(x, colour = NULL) {
  if (is.null(colour)) {
    colour <- colorClusters(x[, 1])
  }
    # pdf(opt$output2)
    # setGraphicBasic()
    par(mar = c(4, 12, 3, 2))
    plot(
        x,
        max.strlen = MAX_CHAR_LEN,
        main = "",
        sub = "",
        do.clus.stat = TRUE,
        xlab = "Silhouette width",
        cex.names = 0.8,
        col = colour,
        nmax.lab = 100,
        do.n.k = FALSE,
        axes = FALSE
    )
    mtext(
        paste("Average silhouette width:", round(summary(x)$avg.width, 3)),
        font = 2,
        cex = 1.5,
        line = 1
    )
    plotAxis(1, 0, 1, 0.2)
    # suprLog = dev.off()
}

###################################
#          GAP STATISTICS
###################################

# n_bootstrap: nb of NB_BOOTSTRAP
#' @export
getGapPerPart <- function(x, c, max_cluster = 6, n_bootstrap = 500, verbose = FALSE) {
    # FUN mus have only two args in this order and return a list with an object cluster

    if (isTRUE(verbose)) {
        cat("\nGAP STATISTICS:\n")
    }
    if (c$method == "kmeans" & (max_cluster > 10 | nrow(x) >= 100)) {
        plural <- c("few ", "s")
    } else {
        plural <- c("", "")
    }
    if (c$method == "kmeans" |
        nrow(x) >= 100) {
        cat(paste("It could take a ", plural[1], "minute", plural[2], "...\n", sep = ""))
    }

    gapFun <- function(x, k) list(cluster = getClusters(c, n_cluster = k))
    clusGap(x, FUNcluster = gapFun, K.max = max_cluster, verbose = FALSE, B = n_bootstrap)
}

# x: gap statisctics object
getGapBest <- function(x, method = "Tibs2001SEmax") {
    with(x, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method = method))
}

# Plot the gap statistics width for all clustering possible
#' @export
plotGapPerPart <- function(x, max_cluster = 6, verbose = FALSE) {
    setGraphic()
    # savePdf("gap_statistics.pdf")
    n_cluster <- getGapBest(x)
    gap_k <- round(x$Tab, 3)
    best <- gap_k[, "gap"][n_cluster]
    if (n_cluster < max_cluster) {
        best <- paste(best, ">", gap_k[, "gap"][n_cluster + 1], "-", gap_k[, "SE.sim"][n_cluster + 1])
    }
    plot(
        x,
        arrowArgs = list(
            col = "gray",
            length = 1 / 15,
            lwd = 2,
            angle = 90,
            code = 3
        ),
        type = "b",
        xlim = c(1, max_cluster + 1),
        ylim = c(0, max(x$Tab[, "gap"]) + 0.1),
        col = "grey",
        xlab = "Nb. of clusters",
        ylab = expression(Gap[k]),
        main = "",
        axes = FALSE
    )
    plotBestClustering(
        "Gap statistics method",
        x$Tab[, "gap"],
        " gap value",
        n_cluster,
        0.1,
        1,
        best,
        verbose = verbose,
        max_cluster = max_cluster
    )
    # cat(paste("With a corrected index, optimal number of clusters k =",getGapBest(gap,"firstSEmax"), "\n"))
    # suprLog = dev.off()
}

# Plot the gap between the two function: within and random within average
#' @export
plotGapPerPart2 <- function(x, max_cluster = 6) {
    setGraphic()
    # savePdf("log_w_diff.pdf")
    min_y <- round(min(x$Tab[, c(1, 2)]), 1)
    max_y <- round(max(x$Tab[, c(1, 2)]), 1)
    plot(
        0,
        0,
        xlim = c(1, max_cluster),
        ylim = c(min_y - 0.1, max_y + 0.1),
        type = "n",
        xlab = "Nb. of clusters",
        ylab = "log(within-inertia)",
        axes = FALSE
    )
    title(
        main = "Optimal number of clusters",
        line = 2,
        cex.main = 2
    )
    mtext(
        text = "Gap statistics method",
        font = 3,
        cex = 1.2,
        line = 0.5
    )
    n_cluster <- getGapBest(x)
    abline(
        v = n_cluster,
        col = "gray",
        lty = 2,
        lwd = 2
    )
    lines(seq(1:max_cluster), x$Tab[, 1], type = "b", col = "red")
    lines(seq(1:max_cluster), x$Tab[, 2], type = "b", col = "blue")
    plotAxis(1, 1, max_cluster)
    plotAxis(2, min_y, max_y, 0.1)
    legend(
        "topright",
        c("log(W)", "E.log(W)"),
        col = c("red", "blue"),
        lty = 1,
        box.lwd = -1,
        bg = "white"
    )
    # suprLog = dev.off()
}

#' @export
printSummary <- function(between, diff, sil, gap = NULL) {
    # TODO: no max_cluster = nrow(data)
    summary <- cbind(between, diff, 100 - between, sil)
    names <- c(
        "Between-inertia (%)",
        "Between-differences (%)",
        "Within-inertia (%)",
        "Silhouette index"
    )
    if (!is.null(gap)) {
        summary <- cbind(summary, gap$Tab[, "gap"][-1], gap$Tab[, "SE.sim"][-1])
        names <- c(names, "Gap", "Gap SE")
    }
    rownames(summary) <- seq(2, (length(between) + 1))
    colnames(summary) <- names
    return(summary)
}

################################
#          HEATMAP
################################

# Inputs:
# x: vector of size for each clusters
plotRect <- function(x, colors, lwd = 3) {
    # size of each clusters
    temp_size <- 0
    for (i in 1:length(x)) {
        # y begin at the top, so sum(x) must be substracted to y coord.
        # rect(xleft, ybottom, xright, ytop)
        # +0.5 because x, y coord are shifted to 0.5 comparativly to plotcolors functions
        rect(
            temp_size + 0.5,
            sum(x) - temp_size - x[i] + 0.5,
            x[i] + temp_size + 0.5,
            sum(x) - temp_size + 0.5,
            border = colors[i],
            lwd = lwd
        )
        # memorize the size of the cluster (for a bottom-right shift)
        temp_size <- temp_size + x[i]
    }
}

# Outputs:
# lenght of clusters ordered by the clusters order
getOrderedClusterSize <- function(x) {
    nb_cl <- length(levels(as.factor(x)))
    size_cl <- rep(0, nb_cl)
    temp_cl <- rep(0, length(x))
    j <- 0

    for (i in 1:length(x)) {
        if (!x[i] %in% temp_cl) {
            j <- j + 1
        }
        size_cl[j] <- size_cl[j] + 1
        temp_cl[i] <- x[i]
    }
    return(size_cl)
}

# Inputs:
# x: data frame
# dist: a distance object
# s: an organised silhouette object
# c: CAH
# cl: clusters from CAH
#' @export
heatMap <- function(x, dist, s = NULL, c = NULL, cl = NULL, is_png = FALSE, verbose = FALSE) {
    plot.new()
    printProgress(verbose, "Heatmap calculation")
    text <- isTRUE(isTRUE(TEXT) & (nrow(x) < 100))

    if (!is.null(s)) {
        order <- attr(s, "iOrd")
        cl_sizes <- summary(s)$clus.size
        title <- "silhouette\'s scores"
        colors <- colPers(length(levels(as.factor(s[, 1]))))
    } else {
        order <- c$order
        cl_sizes <- getOrderedClusterSize(cl[order])
        title <- "dendrogram"
        colors <- orderColors(c, cl)
    }

    matrix <- as.matrix(dist)
    matrix <- matrix[order, order]
    rownames(matrix) <- rownames(x)[order] -> labels
    # if(tri == TRUE) matrix[!lower.tri(matrix)] = NA
    # image(1:ncol(matrix), 1:ncol(matrix), t(matrix), axes=FALSE, xlab="", ylab="")

    options(warn = -1)
    if (nrow(x) > NB_ROW_MAX) {
        labels <- order
    }
    # png(opt$output4,DIM_PNG, DIM_PNG)
    if (is_png) {
        cex.main <- 5
        cex.legend <- 3
        cex.lab <- 2
        y_top <- 12
        x_lab <- 0.6
        lwd.rect <- 6
    } else {
        # pdf(opt$output4)
        cex.main <- 1.5
        cex.legend <- 0.85
        cex.lab <- 0.7
        y_top <- 8
        x_lab <- 0.5
        lwd.rect <- 3
    }

    par(fig = c(0, 0.9, 0, 1), new = TRUE)
    par(mar = c(1, 8, y_top, 1))
    plotcolors(
        dmat.color(matrix, colors = heat.colors(1000), byrank = FALSE),
        ptype = "image",
        na.color = "red",
        rlabels = FALSE,
        clabels = FALSE,
        border.color = 0
    )
    mtext(
        paste("Distance matrix ordered by", title),
        3,
        line = 6,
        font = 4,
        cex = cex.main
    )
    text(
        -0.5,
        0:(ncol(matrix) - 1) + 1,
        rev(labels),
        xpd = NA,
        adj = 1,
        cex = 0.7
    )
    text(
        0.5:(ncol(matrix) - 0.5),
        ncol(matrix) + 1,
        substr(labels, 0, 20),
        xpd = NA,
        cex = 0.7,
        srt = 65,
        pos = 4
    )
    plotRect(cl_sizes, colors, lwd.rect)
    if (isTRUE(text)) {
        text(expand.grid(1:ncol(matrix), ncol(matrix):1), sprintf("%.2f", matrix), cex = 0.4)
    }

    par(fig = c(0.85, 1, 0.3, 0.8), new = TRUE)
    par(mar = c(5, 0, 4, 0) + 0.1)
    legend_image <- as.raster(matrix(heat.colors(1000), ncol = 1))
    plot(
        c(0, 1),
        c(0, 1),
        type = "n",
        axes = FALSE,
        xlab = "",
        ylab = "",
        main = ""
    )
    rasterImage(legend_image, 0.4, 0, 0.5, 1)
    mtext("   Distance",
        3,
        line = 0.5,
        cex = cex.legend,
        font = 2
    )
    text(
        x = x_lab,
        y = seq(0, 1, l = 3),
        labels = round(seq(max(matrix), 2, l = 3)),
        cex = cex.lab,
        pos = 4
    )

    options(warn = 0)
    if (verbose) {
        cat("done.\n")
    }
    # suprLog = dev.off()
}

################################
#          Dendrogram
################################

# Inputs:
# n_cluster: number of clusters
# x: clustering object
#' @export
plotDendrogram <- function(x, cl, n_cluster = 2, max_cluster = 6) {
    if (length(cl) > NB_ROW_MAX) {
        x$labels <- 1:length(cl)
        cex <- 0.4
    } else {
        cex <- 0.8
    }
    # pdf(opt$output7)
    setGraphicBasic()
    par(mar = c(2, 5, 5, 1))
    plot(
        x,
        hang = -1,
        ylim = c(0, max(x$height)),
        xlim = c(0, length(x$labels)),
        sub = "",
        cex = cex,
        font = 3,
        ylab = "Cophenetic distance",
        main = "Dendrogram",
        axes = FALSE
    )
    plotAxis(2, 0, max(x$height))
    abline(h = rev(x$height)[1:max_cluster], col = "gray", lty = 2, lwd = 1)
    # projection of the clusters
    rect.hclust(x, k = as.numeric(n_cluster), border = orderColors(x, cl))
    # suprLog = dev.off()
}

# Get colors ordered for dendrogram
# x: clustering object
orderColors <- function(x, cl) {
    col_in <- colorClusters(cl)[x$order]
    j <- 1
    col_ordered <- rep(NA, length(table(cl)))
    col_ordered[1] <- col_in[1]
    for (i in 2:length(col_in)) {
        if (col_in[i] != col_in[i - 1]) {
            j <- j + 1
            col_ordered[j] <- col_in[i]
        }
    }
    # vector of color: one by cluster
    return(col_ordered)
}

################################
#            PCA
################################

# nf: number of factorial axis
#' @export
plotPca <- function(x, data, cl, axis1 = 1, axis2 = 2, advanced = FALSE, is_png = FALSE) {
    k <- length(levels(as.factor(cl)))

    if (nrow(data) > NB_ROW_MAX) {
        cpoint <- 0
        cstar <- 0
        cellipse <- 0
        clabel <- 0
        labels <- 1:nrow(data)
    } else {
        cpoint <- 0
        cstar <- 1
        cellipse <- 1
        clabel <- 0
        labels <- rownames(data)
    }

    if (is_png) {
        par(mar = c(0, 0, 18, 0), lwd = 4)
        cex <- 2
        cex.main <- 6
        lwd.line <- 8
        line.main <- 7
    } else {
        # pdf(opt$output3)
        par(mar = c(0, 0, 4.1, 0))
        cex <- 0.8
        cex.main <- 1.5
        lwd.line <- 2
        line.main <- 1
    }

    title <- paste("Cumulated inertia:", round((x$eig[axis1] + x$eig[axis2]) / sum(x$eig), 4) * 100, "%")
    s.class(
        addaxes = FALSE,
        cbind(x$li[, axis1], x$li[, axis2]),
        ylim = c(min(x$li[, axis2]) + min(x$li[, axis2]) / 4, max(x$li[, axis2]) + max(x$li[, axis2]) / 2),
        xlim = c(min(x$li[, axis1]), max(x$li[, axis1])),
        csub = 1.5,
        as.factor(cl),
        grid = FALSE,
        col = colPers(k),
        clabel = clabel,
        cstar = cstar,
        cellipse = cellipse,
        cpoint = cpoint
    )
    mtext(title, font = 2, line = line.main, cex = cex.main)
    abline(h = 0, v = 0, lty = 2, lwd = lwd.line, col = "grey")
    text(
        x = x$li[, axis1],
        y = x$li[, axis2],
        labels = labels,
        col = colorClusters(cl),
        cex = cex
    )
    # colnames(pca_coord) = c("Chemicals", "Axis 1", "Axis 2")

    if (isTRUE(advanced)) {
        par(fig = c(0.8, 1, 0.82, 1), new = TRUE)
        plotInertiaPca(x, data, x$nf)
    }
    # suprLog = dev.off()
}

# nf: number of inertia bar plot corresponding to factorial axis
plotInertiaPca <- function(x, data, nf = 4) {
    if (nrow(data) > NB_ROW_MAX) {
        r_lim <- c(8, 0, 4, 5)
        r_main_cex <- 2.7
        r_main_text <- 2.4
        lwd.hist <- 40
        line.hist <- 2
    } else {
        # r_lim = c(-0.2, 0.3, 1.1, 1.1);
        r_lim <- c(2, 0, 1, 1)
        r_main_cex <- 0.7
        r_main_text <- 0.6
        lwd.hist <- 10
        line.hist <- 0
    }

    inertia <- round(x$eig / sum(x$eig) * 100, 1)
    par(mar = c(r_lim[1], r_lim[2], r_lim[3], r_lim[4]) + 0.1)
    plot(
        inertia,
        type = "h",
        lwd = lwd.hist,
        lend = 1,
        xlim = c(0, nf + 0.2),
        ylim = c(0, max(inertia + 7)),
        col = "grey75",
        font = 2,
        axes = FALSE,
        xlab = "",
        ylab = ""
    )
    title(
        sub = " Inertia (in %)",
        line = line.hist,
        cex.sub = r_main_cex,
        font.sub = 3
    )
    text(1:nf, inertia[1:nf] + 5, inertia[1:nf], cex = r_main_text)
    par(new = TRUE)
    par(mar = c(0, 0, 0, 0))
    plot(0:1, 0:1, axes = FALSE, type = "n")
    rect(0, 0.1, 0.9, 0.9, border = "grey65")
}

#########################################
#            Variables contribution
#########################################

#' @export
getDistPerVariable0 <- function(x, cl) {
    nb_cl <- max(cl)
    i_rows <- sapply(seq(nb_cl), function(i) {
        which(cl == i)
    })
    sapply(seq(nb_cl), function(i) {
        apply(x[i_rows[[i]], , drop = FALSE], 2, function(j) {
            mean(j, na.rm = TRUE)
        })
    }) %>% t() %>% set_rownames(paste0("G", seq(nb_cl)))
}

# For a given partition (cl) and each variables (dataset columns)
# pondered distance between the centroid of each clusters and the global centroid of the cloud
# Inputs:
# x: data
# cl: clusters object
getDistPerVariable <- function(x, cl) {
    # Distance between the centroid of each variables
    # ponderation by the sd of the variable (=total inertia per var)
    x <- scalecenter(x)
    nb_cl <- length(levels(as.factor(cl)))
    ctr <- matrix(0, nrow = nb_cl, ncol = ncol(x))

    for (i in 1:nrow(x)) {
        # get the group number for each row
        cli <- cl[i]
        # in the dataset, for a metabolite row, loop an each metadabolite column
        # values are affected the corresponding cluster row and metabolite column in ctr
        for (j in 1:ncol(x)) {
            if(!is.na(x[i, j])) {
                ctr[cli, j] <- ctr[cli, j] + x[i, j]
            }
        }
    }
    colnames(ctr) <- colnames(x)
    rownames(ctr) <- unique(cl)
    return(ctr)
}

# For a given partition, relative contributions of each metabolites to inertia of each clusters (CTR)
# The total of the clusters for each column corresponds to PDIS
# Inputs:
# n_cluster: number of clusters
# c: hierarchical classification
# x: data
#' @export
getCtrVar <- function(x, cl, n_cluster = 2) {

    # if NA values appear, scale 0/0 could produce NA values, NA could correspond to 0
    nb_cl <- length(levels(as.factor(cl)))
    ncol <- ncol(x)
    ctr <- getDistPerVariable(x, cl)

    rownames(ctr) <- paste("G", seq(1, n_cluster), sep = "")
    colnames(ctr) <- colnames(x)

    for (i in 1:nb_cl) {
        for (j in 1:ncol(x)) {
            ctr[i, j] <- ctr[i, j]^2 / (nrow(x) * length(cl[cl == i]))
        }
    }

    return(ctr)
}

getCtrVar2 <- function(x, cl, n_cluster = 2, scale = TRUE) {
    ctr <- getCtrVar(x, cl, n_cluster)
    nb_cl <- length(levels(as.factor(cl)))

    if (isTRUE(scale)) {
        ctr_part <- getPdis(x, cl, n_cluster)

        for (i in 1:nb_cl) {
            for (j in 1:ncol(x)) {
                if (scale) {
                    ctr[i, j] <- ctr[i, j] / ctr_part[j]
                }
            }
        }
    }

    return(ctr)
}

# Discriminant power (PDIS)
# Relative contributions of the metabolites to inertia of a partitionning (in %)
# Inputs:
# n_cluster: number of clusters
# c: hierarchical classification
# x: data
getPdis <- function(x, cl, n_cluster = 2) {

    # for each metabolite contribution (in column), sum the n_cluster clusters values
    return(apply(getCtrVar(x, cl, n_cluster), 2, sum))
}

# Inputs:
# max_cluster: number max of clusters
# cl: list of clusters
# x: data
# index: pdis or rho2 calculation
#' @export
getPdisPerPartition <- function(x, cl, max_cluster = 6) {
    pdis_per_partition <- matrix(NA, max_cluster - 1, ncol(x))
    rownames(pdis_per_partition) <- seq(2, max_cluster)
    colnames(pdis_per_partition) <- colnames(x)

    for (k in 2:max_cluster) {
        res <- getPdis(x, cl[[k - 1]], k)
        for (i in 1:length(res)) {
            pdis_per_partition[k - 1, i] <- res[i]
        }
    }
    return(pdis_per_partition)
}
