#' Calculate the minimum amount of surfactant required to form an emulsion
#' @details All densities expressed in grams per cubic centimeter (g/cc)
#' @return Amount expressed in grams of surfactant required to form emulsion
#' @param rd The density of the dispersed phase (g/cc)
#' @param rs The density of the surfactant mixture (g/cc)
#' @param qc The percentage of the continuous phase
#' @param hlbr The required HLB
#' @export
Qs <- function(rd, rs, qc, hlb) {
  6*(rs/rd) / (10 - 0.5*rhlb) + (4*qc/1000)
}

#' Calculate the rate of particle settling in a disperse system using
#' Stoke's Law
#' @details A positive rate reflects a downward creaming, while a negative
#' rate reflects an upward creaming.
#' @param d The particle diameter (cm)
#' @param ri The density of the particle (g/mL)
#' @param re The density of the medium (g/mL)
#' @param nu The viscosity of the medium (g/cm*s)
#' @param G The gravitational constant (980.665 cm/sq.s)
#' @export
StokesRate <- function(d, ri, re, nu, G=980.665) {
  d^2 * (ri - re) * G / (18 * nu)
}

