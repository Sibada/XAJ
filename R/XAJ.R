#' XAJ: An R implementation of three-source Xinanjiang model by Renjun Zhao.
#'
#'
#' @name XAJ
#' @aliases XAJ-package
#' @docType package
#' @references Zhao and Liu, 1995. The Xinanjiang model, Computer Models of Watershed Hydrology, Water Resources Publication, Highlands Ranch, CO (1995), pp. 215-232
#' @keywords Xinanjiang hydrological model

#' @useDynLib XAJ
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  op <- options()

  param.names = c('KC', 'IM', 'WUM', 'WLM', 'WDM', 'C', 'B', 'SM', 'EX', 'KI', 'KG', 'CI', 'CG', 'N', 'NK')
  xaj.opts <- list(
    XAJ.param.range = data.frame(
      lower     = c(0.20, 0.00,   5.0,  10.0,  10.0, 0.05, 0.1, 10.0, 0.50, 0.01, 0.01, 0.50, 0.95, 0.1, 1.0),
      upper     = c(1.50, 0.05,   20.,  90.0,  60.0, 0.20, 0.6, 60.0, 2.00, 0.70, 0.70, 0.90, 0.998,5.0, 6.0),
      row.names = param.names
    ),
    XAJ.param.names = param.names
  )

  toset <- !(names(xaj.opts) %in% names(op))
  if(any(toset)) options(xaj.opts[toset])
}

#' Run Xinanjiang (XAJ) model (three sources, lumped).
#'
#' @description An R implementation of three-source Xinanjiang model
#'              by Renjun Zhao, used for daily streamflow simulation.
#' @param PREC Time series of precipitation (daily)
#' @param EVAP Time series of pan evaporation
#'             or potential evapotranspiration (daily), length must
#'             equal to \code{PREC}
#' @param params Parameters (see below)
#' @param pIUH Parameters of instantaneous unit hydrograph (IUH) by Nash
#'             for routing of surface runoff.
#' @param UH Unit hydrograph for routing of surface runoff. Provided
#'           when not use parameters of IUH.
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
#'          KI,   outflow coefficients of the free water storage to interflow
#'
#'          KG,   outflow coefficients of the free water storage to groundwater
#'
#'          CI,   recession constant of the lower interflow storage
#'
#'          CG,   recession constant of groundwater storage
#'
#'          Area of the basin also necessary. The parameter \code{params} must
#'          be a numeric vector looks like:
#'
#'          \code{c(KC, IM, WUM, WLM, WDM, C, B, SM, EX, KI, KG, CI, CG, Area)}
#'
#'          Each parameter in the vector must in this order. Basin area (km^2)
#'          is attached to the end of the vector.
#'
#'          \code{pIUH} provides the two meters, i.e. N and NK to calculate the
#'          instantaneous unit hydrograph (IUH) by Nash for routing of surface
#'          runoff. \code{pIUH} must looks like: \code{c(N, NK)}.
#'
#'          N means number of reservoirs in the instantaneous unit hydrograph,
#'          while NK means common storage coefficient in the instantaneous unit
#'          hydrograph. N and NK might also be the model parameters to be
#'          calibrate.
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
XAJ <- function(PREC, EVAP, params, pIUH = c(3, 6), UH = NULL) {
  # Create instantaneous unit hydrograph (IUH) by parameters of IUH (pIUH)
  # when UH is not provided.
  if(is.null(UH) || !is.numeric(UH))
    UH <- IUH(pIUH[1], pIUH[2], 20)

  # Run XAJ model.
  out <- data.frame(XAJrun(PREC, EVAP, params, UH))
  names(out) <- c("E", "EU", "EL", "ED", "W", "WU", "WL", "WD",
                  "R", "RS", "RI", "RG", "Q", "QS", "QI", "QG")
  out
}

# Create IUH
IUH <- function(N, NK, len) {
  UH <- pgamma(seq(0, 168, length.out = len + 1),
               N, scale = NK)
  UH <- diff(UH)
  UH <- UH/sum(UH)
  UH
}

#' Get upper and lower boundary of the parameters of XAJ model.
#'
#' @description Get upper and lower boundary of the parameters of XAJ model.
#'              Calibration of XAJ model patameters would follow by those ranges.
#'              Meaning of the parameters see `?XAJ`.
#'
#' @export
XAJ_param_range <- function() {
  options('XAJ.param.range')[[1]]
}

#' Set the names of the parameters.
#' @description Set the names of the parameters to output so that to
#'              see what parameter the parameters are.
#' @param params The vector of the XAJ parameters
#'
#' @export
name_params <- function(params) {
  c(params, names = options('XAJ.param.names')[[1]])
}

#par(mar=c(5,5,4,5) + 0.1)
#bar <- barplot(prec[-(1:365)], axes=F, ylim=c(35,0), col="grey",
#               border=NA, xaxs="i", yaxs="i", space=0)
#axis(4)
#par(new = T)
#plot(ds, tmp$Q,ylim=c(0,1.3), xlab="", ylab="Streamflow", type="l", xaxs="i", yaxs="i")
#mtext("Precipitation",side=4,line=3)


