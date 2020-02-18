Cm <- 0.219
Rgas <- 8.314
kboltz <- 1.38064852

#' @export
logP <- function(a, f, k) {
  sum(a * f + Cm * k)
}

#' @export
pHAN <- HansenParameter <- function(frag, vol) {
  sum(frag^2)/v
}

#' @export
sHAN <- HansenSolubility <- function(pdisp, pdipole, phbond) {
  sqrt(pdisp^2 + pdipole^2 + phbond^2)
}

#' @export
sHIL <- HildebrandSolubility <- function(ecoh=NULL, vol=NULL, ced=NULL) {
  if(ced) return(sqrt(ced))
  sqrt(ecoh/vol)
}

#' @export
dHmix <- function(pFH, phi, temp) {
  pFH * Rgas * temp * prod(phi)
}

#' @export
xFH <- FloryHugginsInteractionParameter <- function(dHmix=NULL, nsolv=NULL, phi=NULL, temp=NULL, ced=NULL, ced12=NULL, vref=NULL) {
  if(ced & ced12 & vref)
    return((vref * sum(ced*phi) - ced12) / (Rgas * temp))
  dHmix/(kboltz * temp * nsolv * phi)
}

#' @export
A12 <- function(d1, d2, partial=F) {
  a12 <- (d1-d2)^2
  if(!partial) return(sum(a12))
  sum(a12[1] + a12[c(2,3)] * 0.25)
}

#' @export
xFHapprx <- function(vol, A12, temp, beta=0) {
  vol * A12 / (Rgas * temp) + beta
}





