Cm <- 0.219
Rgas <- 8.314
kboltz <- 1.38064852

#' @export
logP <- function(a, f, k) {
  sum(a * f + Cm * k)
}

#' @export
DaviesHLB <- function(H, n) {
  7 + sum(H) - n * 0.475
}

#' @export
GriffinHLB <- function(Mh, M) {
  20 * Mh / M
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

#' Calculate the natural log of the mole fraction of a solute at
#' saturation of an ideal solution as a function of the heat of fusion.
#' @param dHf Molar enthalpy of fusion
#' @param tempS Temperature of the solution
#' @param tempF Melting point of the solid
#' @references Hojjati, H., & Rohani, S. (2006). Measurement and Prediction of Solubility of Paracetamol in Water−Isopropanol Solution. Part 2. Prediction. Organic Process Research & Development, 10(6), 1110–1118. doi:10.1021/op060074g 
#' @export
dHfusion <- function(dHf, tempS, tempF) {
  - dHf / Rgas * (1/tempS - 1/tempF)
}

#' @export
A12 <- function(d1, d2, partial=F) {
  a12 <- (d1-d2)^2
  if(!partial) return(sum(a12))
  sum(a12[1] + a12[c(2,3)] * 0.25)
}

#' @export
xFHapprox <- function(vol, A12, temp, beta=0) {
  vol * A12 / (Rgas * temp) + beta
}

