#' Historical decomposition
#'
#' Function to compute a historical variance decomposition of a VAR.
#'
#' @param x A \code{bvar} object, obtained from \code{\link{bvar}}.
#' @param type Character scalar. Whether to use median or mean values.
#' @param shock_override Optional numeric vector or matrix of structural shocks to override the default Cholesky identification.
#'                       Must have dimensions (K-1) x M, where K = Nlags * M + 1 and M is the number of variables.
#' @param ... Not used.
#'
#' @return Returns a numerical array (time, variable, shock) with the results
#'         of the historical decomposition.
#'
#' @keywords BVAR historical decomposition hd
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Access a subset of the fred_qd dataset
#' data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]
#' # Transform it to be stationary
#' data <- fred_transform(data, codes = c(5, 5, 1), lag = 4)
#'
#' # Estimate a BVAR using one lag, default settings and very few draws
#' x <- bvar(data, lags = 1, n_draw = 600L, n_burn = 100L, verbose = FALSE)
#'
#' # Compute historical decomposition via Cholesky (default)
#' hist_decomp(x, type = "mean")
#'
#' # Or using a custom shock from sign-restricted identification
#' # ident <- sign_restr(x, sign_rest = restr_mat, M = x$config$n_vars)
#' # shock <- ident$shock[, 1]
#' # hist_decomp(x, type = "mean", shock_override = shock)
#' }
hist_decomp.bvar <- function(x, type = c("mean", "quantile"), shock_override = NULL, ...) {
  
  type <- match.arg(type)
  
  # If no override provided, use default Cholesky-identified shock
  if (is.null(shock_override)) {
    shock <- matrix(0, x[["meta"]][["K"]] - 1, x[["meta"]][["M"]])
    shock[seq(x[["meta"]][["M"]]), ] <- t(chol(vcov(x, type = type)))
  } else {
    # Use user-provided shock vector/matrix
    shock <- shock_override
  }
  
  comp <- companion(x, type = type)
  eps  <- solve(shock[seq(x[["meta"]][["M"]]), ], t(residuals(x, type = type)))
  
  out <- compute_hd(x, shock, comp, eps)
  return(out)
}

#' @rdname hist_decomp.bvar
#' @export
hist_decomp <- function(x, ...) { UseMethod("hist_decomp", x) }

#' @rdname hist_decomp.bvar
#' @export
hist_decomp.default <- function(x, ...) {
  stop("No methods for class ", paste0(class(x), collapse = " / "), " found.")
}

#' Historical variance decompostions draws
#'
#' Computes a historical variance decomposition (HD) using draws of the
#' orthogonalized vcov-matrix (currently retrieved with Cholesky decomposition).
#'
#' @param x Numeric matrix. Contains a draw of impulse responses
#' obtained from \code{\link{compute_irf}}.
#' @param shock Numeric matrix with the shocks.
#' @param comp Numeric matrix with the coefficients in companion form.
#' @param eps Numeric matrix of the inverse shock times residuals.
#' @param lags,N,M,K Integer scalar. Dimensions of the VAR.
#'
#' @return Returns a numeric array of HDs.
#'
#' @noRd
compute_hd <- function(x, shock, comp, eps,
                       lags = x[["meta"]][["lags"]], N = x[["meta"]][["N"]],
                       M = x[["meta"]][["M"]], K = x[["meta"]][["K"]]) {
  
  kind_of_I <- cbind(diag(M), matrix(0, M, (lags - 1) * M))
  hd_decomp <- array(0, dim = c(M, N + 1, M))
  tmp <- array(0, dim = c(K - 1, N + 1, M))
  
  for(j in seq_len(dim(hd_decomp)[3])) { # Over variables
    eps_tmp <- matrix(0, M, N + 1)
    eps_tmp[j, seq(2, N + 1)] <- eps[j, ]
    
    for(i in seq.int(2, dim(hd_decomp)[2])) { # Over observations
      tmp[, i, j] <- shock %*% eps_tmp[, i] + comp %*% tmp[, i - 1, j]
      hd_decomp[, i, j] <- kind_of_I %*% tmp[, i, j]
    }
  }
  
  return(aperm(hd_decomp, c(2, 3, 1)))
}
