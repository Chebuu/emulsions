#' Calculate shear stress as a function of yield stress (Pa),
#' shear rate (1/s), flow behavior index, and consistency coefficient.
#' @references http://doi.org/10.1016/j.foodhyd.2017.01.013
#' @export
HerschelBulkley <- function(T0, g, k, n) {
  T0 + k * g^n
}

#' Calculate the Kolmogorov length scale (𝜂) of turbulent flow in a
#' solution as a function of kinematic viscosity and mechanical
#' energy of a stirrer.
#' {\displaystyle \eta =\left({\frac {\nu ^{3}}{\varepsilon }}\right)^{1/4}}
#' @export

KolmogorovLength <- function(v, p, m) {
  (v^3 * p / m)^(1/4)
}

#' Calculate the Kolmogorov time ({\tau _{\eta }}) of turbulent flow in a
#' solution as a function of kinematic viscosity and the mechanical
#' energy of a stirrer.
#' 	{\displaystyle \tau _{\eta }=\left({\frac {\nu }{\varepsilon }}\right)^{1/2}}
#' @export
KolmogorovTime <- function(v, p, m) {
  (v * p / m)^(1/2)
}

#' Calculate the Kolmogorov velocity ({u_{\eta }}) of turbulent flow in a
#' solution as a function of kinematic viscosity and mechanical
#' energy of a stirrer.
#' {\displaystyle u_{\eta }=\left(\nu \varepsilon \right)^{1/4}}
#' @export
KolmogorovVelocity <- function(v, p, m) {
  (v * p / m)^(1/2)
}
