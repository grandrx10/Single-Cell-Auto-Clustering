#' statistical Significance for Hierarchical Clustering (SHC) 
#'
#' Implements the Monte Carlo simulation based significance testing
#' procedure for hierarchical clustering described in Kimes et al. (Biometrics, 2017).
#' Statistical significance is evaluated at each node along the hierarchical
#' tree (dendrogram) starting from the root using a Gaussian null hypothesis test. 
#' A corresponding family-wise error rate (FWER) controlling procedure is
#' provided.
#' 
#' @param x a dataset with n rows and p columns, with observations in rows and
#'        features in columns.
#' @param metric a string specifying the metric to be used in the hierarchical 
#'        clustering procedure. This must be a metric accepted by \code{dist}
#'        or "cor" (to specify 1 - Pearson correlation). (default = "euclidean")
#' @param vecmet a function taking two vectors as input and returning a real-valued
#'        number which specifies how dissimilarities should be computed between rows
#'        of the data matrix (\code{x}). If non-NULL, will take precedence over
#'        \code{metric}. Only one of \code{vecmet} or \code{matmet} can be specified.
#'        (default = NULL)
#' @param matmet a function taking a matrix as input and returning an object of class
#'        \code{dist} representing dissimilarities between rows of the data matrix
#'        (\code{x}). If non-NULL, will take precedence over \code{metric}.
#'        Only one of \code{vecmet} or \code{matmet} can be specified.
#'        (default = NULL)
#' @param linkage a string specifying the linkage to be used in the hierarchical 
#'        clustering procedure. This must be a linkage accepted by 
#'        \code{Rclusterpp.hclust} if \code{rcpp=TRUE}, e.g. "ward",
#'        or \code{stats::hclust} if \code{rcpp=FALSE}, e.g. "ward.D2".
#'        (default = "ward.D2")
#' @param l an integer value specifying the power of the Minkowski distance, if 
#'        used. (default = 2)
#' @param alpha a value between 0 and 1 specifying the desired level of the 
#'        test. If no FWER control is desired, simply set alpha to 1. (default = 1)
#' @param n_sim a numeric value specifying the number of simulations at each node
#'        for Monte Carlo testing. (default = 100) 
#' @param n_min an integer specifying the minimum number of observations needed
#'        to calculate a p-value. (default = 10)
#' @param icovest an integer (1, 2 or 3) specifying the null covariance
#'        estimation method to be used. See \code{\link{null_eigval}} for more
#'        details. (default = 1)
#' @param bkgd_pca a logical value specifying whether to use scaled PCA scores
#'        over raw data to estimate the background noise. When FALSE, raw estimate
#'        is used; when TRUE, minimum of PCA and raw estimates is used.
#'        (default = FALSE)
#' @param rcpp a logical value whether to use the \code{Rclusterpp} package.
#'        Will only be used if the \code{Rclusterpp} package is available.
#'        (default = FALSE)
#' @param ci a string vector specifying the cluster indices to be used for 
#'        testing along the dendrogram. Currently, options include: "2CI", 
#'        "linkage". (default = "2CI")
#' @param null_alg a string vector specifying the clustering algorithms that 
#'        should be used to cluster the simulated null ditributions. Currently, options
#'        include: "2means" and "hclust". While "hclust" typically provides greater power
#'        for rotation invariant combinations of metric and linkage function, 
#'        if a non-rotation invariant metric, e.g. Pearson correlation, is used it is
#'        recommended that the "2means" option is specified. Note, \code{null_alg} and
#'        \code{ci} must be of equal length. (default = "hclust")
#' @param ci_idx an integer between 1 and \code{length(ci)} 
#'        specifiying which CI to use for the FWER stopping rule.
#'        This only has an effect if \code{alpha} is specified to a non-default 
#'        value. (default = 1)
#' @param ci_emp a logical value specifying whether to use the empirical
#'        p-value from the CI based on \code{ci_idx} for the FWER stopping rule.
#'        As with \code{ci_idx}, this only has an effect if \code{alpha} is
#'        specified to a non-default value. (default = FALSE)
#' 
#' @return
#' The function returns a \code{shc} S3-object containing the 
#' resulting p-values. The \code{plot} method can be used to generate a dendrogram
#' with the corresponding p-values placed at each merge. The \code{shc}
#' object has following attributes:
#' \itemize{
#' \item{\code{in_mat}}: {the original data matrix passed to the constructor}
#' \item{\code{in_args}}: {a list of the original parameters passed to the constructor}
#' \item{\code{eigval_dat}}: {a matrix with each row containing the sample eigenvalues for a
#'     subtree tested along the dendrogram}
#' \item{\code{eigval_sim}}: {a matrix with each row containing the estimated eigenvalues used to
#'     simulate null data at a subtree tested along the dendrogram}
#' \item{\code{backvar}}: {a vector containing the estimated background variances used for
#'     computing \code{eigval_sim}}
#' \item{\code{nd_type}}: {a vector of length n-1 taking values in "\code{n_small}",
#'     "\code{no_test}", "\code{sig}", "\code{not_sig}", or "\code{NA}" specifying how each
#'     node along the dendrogram was handled by the iterative testing procedure}
#' \item{\code{ci_dat}}: {a matrix containing the cluster indices for the original
#'     data matrix passed to the constructor}
#' \item{\code{ci_sim}}: {a 3-dimensional array containing the simulated cluster
#'     indices at each subtree tested along the dendrogram}
#' \item{\code{p_emp}}: {a matrix containing the emprical p-values computed at each
#'     subtree tested along the dendrogram}
#' \item{\code{p_norm}}: {a matrix containing the Gaussian approximate p-values
#'     computed at each subtree tested along the dendrogram}
#' \item{\code{idx_hc}}: {a list of tuples containing the indices of clusters joined
#'     at each step of the hierarchical clustering procedure}
#' \item{\code{hc_dat}}: {a \code{hclust} object constructed from the original data matrix and
#'     arguments passed to the constructor}
#' }
#' 
#' @details
#' When possible, the \code{Rclusterpp} should be used for
#' hierarchical clustering by specifying \code{rcpp = TRUE}, except in the case of clustering by Pearson
#' correlation (\code{dist="cor"}), for which we make use of
#' \code{WGCNA::cor} with the usual \code{stats::hclust}.
#'
#' The \code{Rclusterpp} package is no longer available on CRAN and must be
#' installed from GitHub using \code{devtools::install_github("nolanlab/Rclusterpp")}.
#' 
#' For standard minimum variance Ward's linkage clustering, if \code{rcpp} is
#' \code{TRUE}, specify "ward" for \code{linkage}. However, if \code{rcpp} is
#' \code{FALSE}, then "ward.D2" should be specified. See \code{stats::hclust}
#' for details on changes to the function since R >= 3.0.4. 
#' 
#' By default, p-values are computed for all nodes (\code{alpha = 1}).
#' If the procedure should instead terminate when no further nodes
#' meet a desired FWER control threshold, simply specify \code{alpha < 1}.
#' Controlling the FWER  may substantially speed up
#' the procedure by reducing the number of tests considered.
#'
#' Even if the procedure is run using \code{alpha = 1}, the FWER cutoffs may still be
#' computed using \code{fwer_cutoff}. The FWER thresholding can also be visualized by
#' specifying \code{alpha} when calling \code{plot}.
#'
#' The input \code{metric} can either be a character string specifying a metric
#' recognized by \code{dist()} or \code{"cor"} for Pearson correlation.
#' Alternatively, a user-defined dissimilarity function can be specified suing either
#' the \code{matmet} or \code{vecmet} parameter. If specified, \code{matmet} should be a
#' function which takes a matrix as input and returns an object of class \code{dist}
#' from the columns of the matrix. If specified, \code{vecmet} must be a function
#' which computes a real-valued dissimilarity value between two input vectors. This
#' function will be used to compute the dissimilarilty between columns of the data matrix.
#' If either \code{matmet} or \code{vecmet} is specified, \code{metric} will be ignored.
#' If both \code{matmet} and \code{vecmet} are specified, the function exit with an error.
#' Examples using \code{metric}, \code{matmet}, and \code{vecmet} are provided below.
#' 
#' @examples
#' ## using a string input to metric
#' data <- rbind(matrix(rnorm(100, mean = 2), ncol = 2),
#'               matrix(rnorm(100, mean = -1), ncol = 2))
#' shc_metric <- shc(data, metric = "cor", linkage = "average", alpha=0.05)
#' tail(shc_metric$p_norm, 10)
#'
#' ## using a function input to vecmet or any function not optimized for 
#' ## computing dissimilarities for matrices will be incredibly slow
#' ## should be avoided if possible
#' \dontrun{
#' vfun <- function(x, y) { 1 - cor(x, y) }
#' shc_vecmet <- shc(data, vecmet = vfun, linkage = "average", alpha=0.05)
#' tail(shc_vecmet$p_norm, 10)
#' }
#' 
#' ## using a function input to matmet
#' mfun <- function(x) { as.dist(1-cor(x)) }
#' shc_mfun <- shc(data, matmet=mfun, linkage="average", alpha=0.05)
#' tail(shc_mfun$p_norm, 10)
#'
#' @references
#' \itemize{
#'     \item Kimes, P. K., Hayes, D. N., Liu Y., and Marron, J. S. (2016)
#'           Statistical significance for hierarchical clustering.
#'           pre-print available.
#' }
#'
#' @export
#' @seealso \code{\link{plot-shc}} \code{\link{diagnostic}}
#' @import WGCNA methods
#' @importFrom stats as.dendrogram as.dist cutree dist dnorm hclust kmeans rnorm sd pnorm
#' @name shc
#' @aliases shc-constructor
#' @author Patrick Kimes
shc <- function(x, metric = "euclidean", vecmet = NULL, matmet = NULL,
                linkage = "ward.D2", l = 2,
                alpha = 1, icovest = 1, bkgd_pca = FALSE, n_sim = 100,
                n_min = 10, rcpp = FALSE, ci = "2CI", null_alg = "hclust",
                ci_idx = 1, ci_emp = FALSE) {  
  
  n <- nrow(x)
  p <- ncol(x)
  
  if (n < 3) {
    stop("n must be >= 3")
  }
  
  n_ci <- length(ci)
  if (length(null_alg) != n_ci) {
    stop("ci and null_alg must be of same length")
  }
  
  for (ii in 1:n_ci) {
    if (ci[ii] == "linkage" && null_alg[ii] == "2means")
      stop("ci = 'linkage', null_alg = '2means' cannot be specified")
  }
  
  if (ci_idx > n_ci) {
    stop("invalid choice for ci_idx; ci_idx must be < length(ci)")
  }
  
  if (alpha > 1 || alpha < 0) {
    stop("invalid choice for alpha; alpha must be 0 < alpha < 1")
  }
  
  if (!is.matrix(x)) {
    stop("x must be a matrix; use as.matrix if necessary")
  }
  
  if (n_min < 3) {
    stop("n_min must be >= 3")
  }
  
  if (n_min > n) {
    stop("n_min must be <= n")
  }
  
  if (!is.null(vecmet) && !is.null(matmet)) {
    stop("only one of vecmet and matmet can be specified")
  }
  
  if (!is.null(vecmet)) {
    if (!is.function(vecmet)) {
      stop(paste("vecmet must be a function taking two vectors as input",
                 "and returning a real-valued dissimilarity"))
    }
    metric <- NULL
  }
  
  if (!is.null(matmet)) {
    if (!is.function(matmet)) {
      stop(paste("matmet must be a function taking a data matrix as input",
                 "and returning an object of class dist"))
    }
    metric <- NULL
  }
  
  if (rcpp && !requireNamespace("Rclusterpp", quietly = TRUE)) {
    stop("'Rclusterpp' package is not available.\n",
         "Either specify 'rcpp = FALSE' or install 'Rclusterpp' from GitHub using:\n",
         "> devtools::install_github('nolanlab/Rclusterpp')")
  }
  
  ## test vecmet and assign matmet if vecmet specified
  if (!is.null(vecmet)) {
    tryCatch({
      tmp <- vecmet(x[1, ], x[2, ])
    }, warning = function(e) {
      stop(paste0("warning for vecmet specification: ", e))
    }, error = function(e) {
      stop(paste0("error with vecmet specification: ", e))
    })
    matmet <- function(x) {
      as.dist(outer(split(x, row(x)), split(x, row(x)),
                    Vectorize(vecmet)))
    }
  }
  
  ## rclusterpp doesn't recognize 'ward.D2', stop and let user know
  if ((linkage == "ward.D2") & rcpp) {
    stop("Use 'ward' (in place of 'ward.D2') for linkage when rcpp = TRUE.")
  }
  
  ## apply initial clustering
  x_clust <- .initcluster(x, n, p, metric, matmet, linkage, l, 
                          n_ci, ci, rcpp)
  ci_dat <- x_clust$ci_dat
  hc_dat <- x_clust$hc_dat
  idx_hc <- x_clust$idx_hc
  
  ## for plotting purposes, change heights of dendrogram
  if ((linkage == "ward") & rcpp) {
    hc_dat$height <- sqrt(2*hc_dat$height)
  }
  
  ## p-values for all <= (n-1) tests
  p_emp <- matrix(2, nrow=n-1, ncol=n_ci)
  p_norm <- matrix(2, nrow=n-1, ncol=n_ci)
  colnames(p_emp) <- paste(null_alg, ci, sep="_")
  colnames(p_norm) <- paste(null_alg, ci, sep="_")
  
  ## null covariance parameters for all <= (n-1) tests
  eigval_dat <- matrix(-1, nrow=n-1, ncol=p)
  eigval_sim <- matrix(-1, nrow=n-1, ncol=p)
  backvar <- rep(-1, n-1)
  ci_sim <- array(-1, dim=c(n-1, n_sim, n_ci))
  
  ## determine parent nodes for all nodes
  pd_map <- .pd_map(hc_dat, n)
  
  ## compute Meinshausen cutoffs for significance at alpha
  cutoff <- fwer_cutoff(idx_hc, alpha)
  
  ## keep track of each node was tested
  nd_type <- rep("", n-1)
  
  level_index <- rep(0, n)
  # Cluster markers for each cell
  clusters <- rep(NULL, length(x))
  cluster_num <- 0
  
  child_clusters <- rep(NULL, length(x))
  child_borders <- rep(NULL, length(x))
  print(length(child_borders))
  latest_borders <- rep(NA, length(x))
  
  ## move through nodes of dendrogram
  for (k in 1:(n-1)) {
    if (pd_map[k] == 0){
      level <- 0
      level_index[k] <- level + 1
    } else {
      level <- level_index[pd_map[k]]
      level_index[k] <- level + 1
    }
    curr_cutoff <- alpha/(2^(level_index[k]) - 1)
    # curr_cutoff <- cutoff[k]
    
    ## indices for subtree
    idx_sub <- unlist(idx_hc[k, ])
    
    n_sub <- length(idx_sub)
    
    ## only calc p-values for branches w/ more than n_min
    if (n_sub < n_min) {
      nd_type[k] <- "n_small"
      next
    }
    
    ## if parent wasn't significant, skip
    ## - placed after n_min check on purpose
    if ((alpha < 1) && (k > 1) && (nd_type[pd_map[k]] != "sig")) {
      nd_type[k] <- "no_test"
      next
    }
    
    ## estimate null Gaussian
    xk_null <- NULL #null_eigval(x[idx_sub, ], n_sub, p, icovest, bkgd_pca)
    
    ## compute p-values
    m_idx <- NULL#<- colMeans(as.matrix(ci_sim[k, , ]))
    s_idx <- NULL #<- apply(as.matrix(ci_sim[k, , ]), 2, sd)
    
    ## keep everything
    eigval_dat <-NULL #xk_null$eigval_dat
    eigval_sim <-NULL #xk_null$eigval_sim
    backvar <-NULL #xk_null$backvar
    
    hc_isim <- .cluster_shc(x[idx_sub, ], metric, matmet, linkage, l, rcpp)
    split <- cutree(hc_isim, k=2)
    cluster1 <- idx_sub[split == 1]
    cluster2 <- idx_sub[split == 2]
    # print("cluster1 and 2 size:")
    # print(paste(length(cluster1), length(cluster2)))
    
    pval <- 1
    if (length(cluster1) > 1 && length(cluster2) > 1) {
      ## Store or process clusters as needed
      cluster1_obs <- x[cluster1, , drop=FALSE]
      cluster2_obs <- x[cluster2, , drop=FALSE]

      # Find border cells for both clusters
      results_1 <- find_border_cells_simple(cluster1_obs, cluster2_obs)
      results_2 <- find_border_cells_simple(cluster2_obs, cluster1_obs)
      
      # Extract border cells and updated observations
      border_cells_1 <- results_1$border_cells
      
      border_cells_2 <- results_2$border_cells
      
      if (is.null(results_1$error_message) && is.null(results_2$error_message)) {
        # Get the nearest points in updated cluster1 (border points)
        border_points_1 <- cluster1_obs[border_cells_1, , drop = FALSE]
        
        # Exclude border cells from updated cluster1 observations
        non_border_cells_1 <- cluster1_obs[-border_cells_1, , drop = FALSE]
        
        if (nrow(non_border_cells_1) > 0) {
          border_nn_result_1 <- nn2(data = non_border_cells_1, query = border_points_1, k = 1)
          border_to_nonborder_1 <- border_nn_result_1$nn.dists[,1]
        } else {
          print("NO non-border cells for cluster 1")
        }
        
        # -------- Process for updated cluster2 -------------- #
        
        # Get the nearest points in updated cluster2 (border points)
        border_points_2 <- cluster2_obs[border_cells_2, , drop = FALSE]
        
        # Exclude border cells from updated cluster2 observations
        non_border_cells_2 <- cluster2_obs[-border_cells_2, , drop = FALSE]
        # 
        if (nrow(non_border_cells_2) > 0) {
          border_nn_result_2 <- nn2(data = non_border_cells_2, query = border_points_2, k = 1)
          border_to_nonborder_2 <- border_nn_result_2$nn.dists[, 1]
        }else {
          print("NO non-border cells for cluster 2")
        }
        
        # Nearest neighbors from updated cluster1 to updated cluster2 and vice versa
        nn_result_1 <- nn2(data = cluster2_obs, query = border_points_1, k = 1)
        border_to_other_1 <- nn_result_1$nn.dists
        
        nn_result_2 <- nn2(data = cluster1_obs, query = border_points_2, k = 1)
        border_to_other_2 <- nn_result_2$nn.dists
        
        #   
        #   
        #   # Exclude border cells from cluster1 observations
        
        if (exists("border_to_nonborder_1") && exists("border_to_nonborder_2")
            && exists("border_to_other_1") && exists("border_to_other_2")) {
          border_to_other <- c(border_to_other_1, border_to_other_2)
          border_to_self <- c(border_to_nonborder_1, border_to_nonborder_2)
          wilcox_result <- wilcox.test(border_to_other, border_to_self, alternative = "greater", paired=TRUE)
          pval <- wilcox_result$p.value
        }
      }
      
      
    
    } else {
      print("Clusters are not large enough.")
    }
  
    print("--------------------------------")
    p_norm[k, ] <- pval
    p_emp[k, ] <- pval
    p_norm[k, ci == "linkage"] <- pval
    p_emp[k, ci == "linkage"] <- pval
    
    ## update nd_type (node type)
    if (alpha < 1) {
      nd_type[k] <- ifelse(pval < curr_cutoff, #p_norm[k, ci_idx] < cutoff[k],
                             "sig", "not_sig")
      if (nd_type[k] == "not_sig") {
        print("main clustering... not significant.")
        clusters[cluster1] = cluster_num
        # cluster_num = cluster_num + 1
        child_clusters[cluster1] = paste0(cluster_num, "a") 
        clusters[cluster2] = cluster_num
        child_clusters[cluster2] = paste0(cluster_num, "b") 
        
        border_indices_in_x_1 <- cluster1[border_cells_1]
        child_borders[border_indices_in_x_1] = paste0(cluster_num, "a")
        
      
        border_indices_in_x_2 <- cluster2[border_cells_2]
        child_borders[border_indices_in_x_2] = paste0(cluster_num, "b")
        print(length(child_borders))
        
        cluster_num = cluster_num + 1
      } else { # if you were significantly different
        
        border_indices_in_x_1 <- cluster1[border_cells_1]
        latest_borders[cluster1] = NA
        latest_borders[border_indices_in_x_1] = paste0(cluster_num, "a")
        
        border_indices_in_x_2 <- cluster2[border_cells_2]
        latest_borders[cluster2] = NA
        latest_borders[border_indices_in_x_2] = paste0(cluster_num, "b")
      }
    } else {
      nd_type[k] <- "cutoff_skipped"
      
    }
  }

  # cluster1_obs <- x[clust1, , drop=FALSE]
  # cluster2_obs <- x[clust2, , drop=FALSE]
  # 
  # nn_result <- nn2(data = cluster1_obs, query = cluster2_obs, k = 1)
  # nearest_indices <- nn_result$nn.idx
  # idx <- clust1[nearest_indices]
  # clusters[idx] = 50
  
  # clusters[clust1] = 41
  # clusters[clust2] = 42
  # 
  # nn_result <- nn2(data = cluster1_obs, query = cluster2_obs, k = 1)
  # # Get the indices of the nearest neighbors in cluster1 for each point in cluster2
  # nearest_indices_1 <- nn_result$nn.idx
  # border_cells_1 <- unique(nearest_indices_1)
  # 
  # border_idx_1 <- clust1[nearest_indices_1]
  # clusters[border_idx_1] <- 41
  # 
  # non_border_cells <- cluster1_obs[-border_cells_1, , drop = FALSE]
  # border_points_1 <- cluster1_obs[border_cells_1, , drop = FALSE]
  # border_nn_result <- nn2(data = non_border_cells, query = border_points_1, k = 1)
  # idx <- border_nn_result$nn.idx
  # # Find the number of border points
  # num_border_points <- nrow(border_points_1)
  # 
  # # Print the number of border points
  # non_border_cells_1 <- clust1[idx]
  # non_border_cells_1 <- unique(non_border_cells_1)
  # clusters[non_border_cells_1] <- 42
  # # 
  # # # SECOND CLUSTER
  # nn_result <- nn2(data = cluster2_obs, query = cluster1_obs, k = 1)
  # # Get the indices of the nearest neighbors in cluster1 for each point in cluster2
  # nearest_indices_2 <- nn_result$nn.idx
  # border_cells_2 <- unique(nearest_indices_2)
  # 
  # border_idx_2 <- clust2[nearest_indices_2]
  # clusters[border_idx_2] <- 43
  # 
  # non_border_cells <- cluster2_obs[-border_cells_2, , drop = FALSE]
  # border_points_2 <- cluster2_obs[border_cells_2, , drop = FALSE]
  # border_nn_result <- nn2(data = non_border_cells, query = border_points_2, k = 1)
  # idx <- border_nn_result$nn.idx
  # non_border_cells_2 <- clust2[idx]
  # non_border_cells_2 <- unique(non_border_cells_2)
  # clusters[non_border_cells_2] <- 44
  # 
  # for (i in 1:length(clusters)) {
  #   if (is.na(clusters[i]) || clusters[i] > 3) {
  #     clusters[i] = 4
  #   }
  # }
  stop()
  
  
  ## return shc S3 object
  structure(
    list(in_mat = x,
         in_args = list(metric = metric, linkage = linkage, alpha = alpha,
                        l = l, bkgd_pca = bkgd_pca, n_sim = n_sim,
                        n_min = n_min, icovest = icovest, ci = ci,
                        null_alg = null_alg, ci_idx = ci_idx, ci_emp = ci_emp),
         eigval_dat = eigval_dat,
         eigval_sim = eigval_sim,
         backvar = backvar,
         nd_type = nd_type,
         ci_dat = ci_dat,
         ci_sim = ci_sim,
         p_emp = p_emp,
         p_norm = p_norm,
         idx_hc = idx_hc,
         hc_dat = hc_dat,
         clusters = clusters,
         child_clusters = child_clusters,
         child_borders = child_borders,
         latest_borders = latest_borders),
    class = "shc")
}

## #############################################################################
## #############################################################################
## helper functions
find_border_cells <- function(cluster1_obs, cluster2_obs) {
  # Load the RANN package
  library(RANN)
  
  # Make a copy of the original cluster1_obs to preserve the original data
  original_cluster1_obs <- cluster1_obs
  
  # Initialize vector to store indices of border cells
  border_indices <- integer(0)
  
  # Create a vector to keep track of the original indices
  original_indices <- seq_len(nrow(cluster1_obs))
  
  repeat {
    # Check for empty data conditions
    if (nrow(cluster1_obs) < 1 || nrow(cluster2_obs) < 1) {
      print("No data remaining! A cluster has been emptied")
      return(list(
        border_cells = seq_len(nrow(original_cluster1_obs)), error_message = "NOT ENOUGH DATA"
      ))
    }
    
    # Find nearest neighbors in cluster1 for each cell in cluster2
    nn <- nn2(data = cluster1_obs, query = cluster2_obs, k = 1)
    nearest_neighbors <- nn$nn.idx[, 1]
    
    # Calculate the frequency of each neighbor
    neighbor_counts <- table(nearest_neighbors)
    
    # Identify neighbors that are shared by more than 10% of cluster2 cells
    overrepresented <- as.numeric(names(neighbor_counts[neighbor_counts > 0.5 * nrow(cluster2_obs)]))
    
    # If no overrepresented neighbors are found, exit the loop
    if (length(overrepresented) == 0) break
    
    # Add overrepresented neighbors to the border_indices (using original indices)
    border_indices <- c(border_indices, original_indices[overrepresented])
    
    # Remove the overrepresented cells from cluster1
    cluster1_obs <- cluster1_obs[-overrepresented, , drop = FALSE]
    original_indices <- original_indices[-overrepresented]
    print("Overrepresented detection.")
  }
  
  # Find nearest neighbors in the final remaining cluster1_obs for each cell in cluster2_obs
  nn_final <- nn2(data = cluster1_obs, query = cluster2_obs, k = 1)
  nearest_neighbors_final <- nn_final$nn.idx[, 1]
  
  # Add final nearest neighbors to the border_indices (using original indices)
  border_indices <- unique(c(border_indices, original_indices[nearest_neighbors_final]))
  print("Cluster size:")
  print(length(original_indices))
  print("BORDER LENGTH:")
  print(length(border_indices))
  # Use the original_cluster1_obs to extract border cells by their indices
  # border_cells <- original_cluster1_obs[border_indices, , drop = FALSE]
  
  # Return the border cells and any error message
  return(list(border_cells = border_indices, error_message = NULL))
}


find_border_cells_simple <- function(cluster1_obs, cluster2_obs, distance_threshold = 1.2) {
  library(RANN)
  
  # Step 1: Calculate the centroid of cluster1 in the PC space
  cluster1_centroid <- colMeans(cluster1_obs)
  
  # Step 2: Calculate the distances of all cells in cluster1 to the centroid
  cluster1_distances <- sqrt(rowSums((cluster1_obs - cluster1_centroid)^2))
  
  # Step 3: Filter out cells that are farther from the centroid than the threshold
  avg_distance <- mean(cluster1_distances)
  filtered_cluster1_obs <- cluster1_obs[cluster1_distances <= (avg_distance * distance_threshold), , drop = FALSE]
  
  # Save the original indices of the filtered cells
  filtered_indices <- which(cluster1_distances <= (avg_distance * distance_threshold))
  
  # Step 4: Use RANN to find the nearest neighbor of each cell in cluster2 relative to the filtered cluster1
  # Use RANN::nn2 to find nearest neighbors
  nn_result <- nn2(data = filtered_cluster1_obs, query = cluster2_obs, k = 1)
  
  # Get the indices of the nearest neighbors in the filtered cluster1
  nearest_neighbor_indices <- nn_result$nn.idx[, 1]
  
  # Map the nearest neighbor indices back to the original indices of cluster1
  border_cell_indices <- filtered_indices[nearest_neighbor_indices]
  
  return (list(border_cells = border_cell_indices, error_message = NULL))
}






# find_border_cells <- function(cluster1_obs, cluster2_obs) {
#   # Load the RANN package
#   library(RANN)
#   
#   # Check if the input data dimensions match
#   if (ncol(cluster1_obs) != ncol(cluster2_obs)) {
#     print("Dimension Mismatch problem.")
#     return(list(
#       cluster1_obs = NULL, nearest_cells = NULL, error_message = "DIMENSION MISMATCH"
#     ))
#   }
#   
#   # Main loop to find and remove overrepresented neighbors
#   repeat {
#     # Check for empty data conditions
#     if (nrow(cluster1_obs) < 1 || nrow(cluster2_obs) < 1) {
#       print("No data remaining! A cluster has been emptied")
#       return(list(
#         cluster1_obs = NULL, nearest_cells = NULL, error_message = "NOT ENOUGH DATA"
#       ))
#     }
#     
#     # Find nearest neighbors in cluster1 for each cell in cluster2
#     nn <- nn2(data = cluster1_obs, query = cluster2_obs, k = 1)
#     nearest_neighbors <- nn$nn.idx[, 1]
#     
#     # Calculate the frequency of each neighbor
#     neighbor_counts <- table(nearest_neighbors)
#     
#     # Identify neighbors that are shared by more than 10% of cluster2 cells
#     overrepresented <- names(neighbor_counts[neighbor_counts > 0.3 * nrow(cluster2_obs)])
#     
#     nearest_neighbors <- unique(nearest_neighbors)
#     
#     # If no overrepresented neighbors are found, exit the loop
#     if (length(overrepresented) == 0) break
#     
#     # Remove the overrepresented cells from cluster1
#     cluster1_obs <- cluster1_obs[-as.numeric(overrepresented), , drop = FALSE]
#     print("Overrepresented detection.")
#   }
#   
#   # Return the updated cluster1_obs and the matching nearest neighbor indices
#   return(list(cluster1_obs = cluster1_obs, nearest_cells = nearest_neighbors, error_message = NULL))
# }



is_subset <- function(subset_vec, main_vec) {
  return(all(subset_vec %in% main_vec))
}
## identify parent node of each node in dendrogram
.pd_map <- function(hc, n) {
  ## determine parent branch node for all children nodes along dendrogram
  pd_pairs <- rbind(cbind(hc$merge[, 1], 1:(n-1)), 
                    cbind(hc$merge[, 2], 1:(n-1)))
  pd_map <- data.frame(pd_pairs[pd_pairs[, 1] > 0, ])
  names(pd_map) <- c("dtr", "prt")
  pd_map <- pd_map$prt[order(pd_map$dtr)] #the parent of each daughter
  pd_map <- c(pd_map, n) #add final node without a parent
  
  ## flip index, hclust and shc use reversed ordering
  n - rev(pd_map)
}


## determine obs indices at each node of the dendrogram
.idx_hc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idx_hc <- array(list(), c(2*n-1, 2))
  idx_hc[1:n, 1] <- as.list(n:1)
  idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idx_hc
  for (k in 1:(n-1)) {
    idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
    idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idx_hc[(2*n-1):(n+1), ]
}


## calculate sum of squares
.sumsq <- function(x) { norm(sweep(x, 2, colMeans(x), "-"), "F")^2 }


## calculate 2-means cluster index (n x p matrices)
.calc2CI <- function(x1, x2) {
  if (is.matrix(x1) && is.matrix(x2) && ncol(x1) == ncol(x2)) {
    (.sumsq(x1) + .sumsq(x2)) / .sumsq(rbind(x1, x2))
  } else {
    stop(paste("x1, x2 must be matrices with same ncols",
               "for 2CI calculation"))
  }      
}


## parse clustering parameters to produce hclust object
.cluster_shc <- function(x, metric, matmet, linkage, l, rcpp) {
  num_cells <- nrow(x)
  k_value <- min(20, num_cells - 1)  # Default to 20 or smaller than number of cells
  
  # Run SNN construction with adjusted k
  data.SNN <- SNN.Construction(mat = x, k = k_value)
  
  # Proceed with hierarchical clustering
  hc_dat <- HGC.dendrogram(G = data.SNN)
  
  hc_dat
}


## perform hierarchical clustering on the original data and 
## compute the corresponding cluster indices for each merge
.initcluster <- function(x, n, p, metric, matmet, linkage, l, 
                         n_ci, ci, rcpp) {
  
  ## obtain clustering solution
  hc_dat <- .cluster_shc(x, metric, matmet, linkage, l, rcpp)
  
  ## list array of cluster indices at each of the n-1 nodes
  idx_hc <- .idx_hc(hc_dat, n)
  
  ## matrix containing cluster indices
  ci_dat <- matrix(-1, nrow=n-1, ncol=n_ci)
  
  ## calculate cluster index(ices) for merge k
  for (i_ci in 1:n_ci) {
    if (ci[i_ci] == "2CI") {
      for (k in 1:(n-1)) {
        ci_dat[k, i_ci] <- .calc2CI(x[idx_hc[[k, 1]], , drop=FALSE],
                                    x[idx_hc[[k, 2]], , drop=FALSE])
      }
    } else if (ci[i_ci] == "linkage") {
      ## flip index, hclust and shc use revered ordering
      ci_dat[, i_ci] <- rev(hc_dat$height)
    }
  }
  
  list(hc_dat = hc_dat, 
       idx_hc = idx_hc,
       ci_dat = ci_dat)
}


## given null eigenvalues, simulate Gaussian dataset
.simnull <- function(eigval_sim, n, p) {
  simnorm <- matrix(rnorm(n*p, sd=sqrt(eigval_sim)), n, p, byrow=TRUE)
}


## perform hierarchical clustering on a simulated dataset and
## compute the correspond cluster indices for only the final merge
.calcCI_shc <- function(x, p, metric, matmet, linkage, l, 
                        n_ci, ci, null_alg, rcpp) {
  
  ##obtain clustering solution
  hc_isim <- .cluster_shc(x, metric, matmet, linkage, l, rcpp)
  split <- cutree(hc_isim, k=2)
  
  ##row vector containing cluster indices
  ci_isim <- matrix(-1, nrow=1, ncol=n_ci)
  
  for (i_ci in 1:n_ci) {
    if (ci[i_ci] == "2CI") {
      if (null_alg[i_ci] == "hclust") {
        ci_isim[i_ci] <- .calc2CI(x[split==1, , drop=FALSE],
                                  x[split==2, , drop=FALSE])        
      } else if (null_alg[i_ci] == "2means") {
        kmsol <- kmeans(x, centers=2)
        ci_isim[i_ci] <- kmsol$tot.withinss/kmsol$totss
      }
    } else if (ci[i_ci] == "linkage") {
      ci_isim[i_ci] <- hc_isim$height[nrow(x)-1]
    }
  }
  
  ci_isim
}








