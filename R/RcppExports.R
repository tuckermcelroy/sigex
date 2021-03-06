# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

polymul_matt <- function(amat, bmat) {
    .Call(`_sigex_polymul_matt`, amat, bmat)
}

VARMA_auto <- function(param, p, q, maxlag) {
    .Call(`_sigex_VARMA_auto`, param, p, q, maxlag)
}

getEigenValues <- function(M) {
    .Call(`_sigex_getEigenValues`, M)
}

complexExp <- function(x) {
    .Call(`_sigex_complexExp`, x)
}

polymult <- function(a, b) {
    .Call(`_sigex_polymult`, a, b)
}

polymul_mat <- function(amat, bmat) {
    .Call(`_sigex_polymul_mat`, amat, bmat)
}

ar_adjoint <- function(poly_array) {
    .Call(`_sigex_ar_adjoint`, poly_array)
}

auto_VARMA <- function(param, p, q, ps, qs, season, grid, maxlag) {
    .Call(`_sigex_auto_VARMA`, param, p, q, ps, qs, season, grid, maxlag)
}

