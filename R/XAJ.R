#' @useDynLib XAJ
#' @importFrom Rcpp sourceCpp
NULL


#' Run Xinanjiang (XAJ) model (three sources, lumped).
#'
#' @description An R implementation of three-source Xinanjiang model
#'              by Renjun Zhao, used for daily streamflow simulation.
#' @param PREC Time series of precipitation (daily)
#' @param EVAP Time series of pan evaporation
#'             or potential evapotranspiration (daily), length must
#'             equal to \code{PREC}
#' @param params Parameters (see below)
#' @param UH Unit hydrographs for routing of surface runoff.
#'
#' @details This function is an R implementation of the Xinanjiang (XAJ)
#'          hydrological model. The lumped XAJ model has 13 parameters, including:
#'
#'          KC,   Ratio of potential evap to pan evap
#'
#'          IM,   Fraction of impermeable area
#'
#'          WUM,  Soil moisture capacity of upper layer
#'
#'          WLM,  Soil moisture capacity of lower layer
#'
#'          WDM,  Soil moisture capacity of deep layer
#'
#'          C,    Coefficient of deep evap
#'
#'          B,    Exponent of the soil moisture storage capacity curve
#'
#'          SM,   Areal mean free water capacity of the surface soil layer
#'
#'          EX,   Exponent of the free water capacity curve
#'
#'          KG,   outflow coefficients of the free water storage to interflow
#'
#'          KI,   outflow coefficients of the free water storage to groundwater
#'
#'          CG,   recession constant of the lower interflow storage
#'
#'          CI,   recession constant of groundwater storage
#'
#'          Area of the basin also necessary. The parameter \code{params} must
#'          be a numeric vector looks like:
#'
#'          \code{c(KC, IM, WUM, WLM, WDM, C, B, SM, EX, KG, KI, CG, CI, Area)}
#'
#'          Each parameter in the vector must in this order. Basin area (km^2)
#'          is attached to the end of the vector.
#'
#' @return This function returns a data frame of some common variables of the XAJ
#'         model at each time step, such as evaporation, soil moisture, surface
#'         and underground runoff.
#'
#'         The variables of the table including:
#'
#'         E,   Total evaporation (mm)
#'
#'         EU,  Evaporation (mm) of upper soil layer
#'
#'         EL,  Evaporation (mm) of lower soil layer
#'
#'         ED,  Evaporation (mm) of deep soil layer
#'
#'         W,   Total soil moisture (mm)
#'
#'         WU,  Soil moisture (mm) of upper soil layer
#'
#'         WL,  Soil moisture (mm) of lower soil layer
#'
#'         WD,  Soil moisture (mm) of deep soil layer
#'
#'         R,   Total runoff (mm) of each time step
#'
#'         RS,  Surface runoff (mm) of each time step
#'
#'         RI,  Interflow (mm) of each time step
#'
#'         RG,  Underground runoff (mm) of each time step
#'
#'         Q,   Total runoff (m^3/s) at the outlet of the basin
#'
#'         QS,  Surface runoff (m^3/s) at the outlet of the basin
#'
#'         QI,  Interflow runoff (m^3/s) at the outlet of the basin
#'
#'         QG,  Underground runoff (m^3/s) at the outlet of the basin
#'
#'
#' @references Zhao and Liu, 1995. The Xinanjiang model, Computer Models of
#'             Watershed Hydrology, Water Resources Publication, Highlands
#'             Ranch, CO (1995), pp. 215-232
#' @export
XAJ <- function(PREC, EVAP, params, UH) {
  out <- data.frame(XAJrun(PREC, EVAP, params, UH))
  names(out) <- c("E", "EU", "EL", "ED", "W", "WU", "WL", "WD",
                  "R", "RS", "RI", "RG", "Q", "QS", "QI", "QG")
  out
}