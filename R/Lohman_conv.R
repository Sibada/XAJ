
#' Convolution of the distributed Xinanjiang model output by the routing
#' model of Lohman et al.
#'
#' @description Convolution of the distributed Xinanjiang model output by
#'              the routing model of Lohman et al.
#' @param Qs Runoff time series of each subbasins.
#' @param uh Unit hydrographs of Lohman model of each subbasins.
#'
#' @return The routing result of the XAJ output, that is, the time series
#'         of the streamflow at the outlet of the river.
#'
#' @export

Lohmann_conv <- function(Qs, uh) {
  tmpQ <- Qs %*% t(uh)
  Q <- aux_Lohmann_conv(tmpQ)
  Q
}

#' Conver flow direction ascii grid to river table.
#' @description Conver flow direction ascii grid to rivers.
#'
#' @param dir.file Path of the flow direction ascii grid file.
#' @param grid.loc A two-column array, a column is longitute
#'                 and another column is latitude of the grids
#'                 to be considered in XAJ model.
#' @param arcinfo If the direction code is arcgis style. Default
#'                TRUE.
#'
#' @return A table containing the informations of the river
#'         network. Each row represents a river channel, the first column
#'         is the next river channel to flow in, the second column is the
#'         length [m] of the river channel, while the third column is the
#'         subbasin flow into that river channel.
#' @export
dir2river <- function(dir.file, grid.loc, arcinfo = TRUE) {
  x <- grid.loc[,1]; y <- grid.loc[,2]

  params <- read.table(dir.file, nrows = 6)
  csize <- params[5, 2]
  ndval <- params[6, 2]
  xll <- params[3, 2]
  xll <- ifelse(params[3,1] == 'xllcorner', xll+csize/2, xll)
  yll <- params[4, 2]
  yll <- ifelse(params[4,1] == 'yllcorner', yll+csize/2, yll)

  grid <- t(as.matrix(read.table(dir.file, skip = 6)))
  grid <- grid[,ncol(grid):1]
  grid[grid == ndval] <- NA

  ix = round((x - xll)/csize + 1)
  iy = round((y - yll)/csize + 1)

  if(arcinfo) {
    nextx <- c('1'=1, '2'=1, '4'=0, '8'=-1, '16'=-1, '32'=-1, '64'=0, '128'=1)
    nexty <- c('1'=0, '2'=-1, '4'=-1, '8'=-1, '16'=0, '32'=1, '64'=1, '128'=1)
  } else {
    nextx <- c('3'=1, '4'=1, '5'=0, '6'=-1, '7'=-1, '8'=-1, '1'=0, '2'=1)
    nexty <- c('3'=0, '4'=-1, '5'=-1, '6'=-1, '7'=0, '8'=1, '1'=1, '2'=1)
  }

  # Find the index of the next grid to flow in of each grid
  ng <- nrow(grid.loc)
  nextg <- rep(0, ng)
  distance <- rep(0., ng)
  for(i in 1:ng) {
    idir <- paste(grid[ix[i], iy[i]])
    nx <- ix[i] + nextx[[idir]]
    ny <- iy[i] + nexty[[idir]]


    lat <- y[i]; nlat <- lat + nexty[[idir]]*csize
    lat <- lat*pi/180; nlat <- nlat*pi/180
    idis <- 6371393 * acos(sin(lat)*sin(nlat)+cos(lat)*
                             cos(nlat)*cos(nextx[[idir]]*csize*pi/180))

    distance[i] <- idis

    inextg <- which(ix == nx & iy == ny)
    if(length(inextg) == 0) {
      nextg[i] <- 0
    } else {
      nextg[i] <- inextg[1]
    }
  }

  out <- data.frame(nextr = nextg, dist = distance, sb = 1:ng)
  return(out)
}

#' Create unit hydrograph matrix of the Lohmann routing model.
#'
#' @description Create unit hydrograph matrix of the Lohmann
#'              routing model. Each column is the UH of each
#'              subbasins while each row is the UH of each time
#'              steps.
#'
#' @param river River network table like the output of function
#'              `dir2river()`.
#'
#' @param stn_loc Hydrological station located in which river channel.
#' @param v Flow velocity of the river channel [m].
#' @param diff Water diffusion coefficient of the river channel.
#'
#' @return A matrix of the UH of each gridcells.
Lohmann_uh <- function(river, stn_loc, v = 1.5, diff = 800) {
  ng <- nrow(river)
  if(length(v) < ng) v <- rep(v[1], ng)
  if(length(diff) < ng) diff <- rep(diff[1], ng)

  # Find out the sub-basins controled by the station
  state <- rep(0, max(river[,3]))
  sbs <- river[river[,3] > 0, 3]
  for(i in sbs) {
    og <- which(river[,3] == i)
    if(state[og] != 0)
      next
    while(TRUE) {
      state[i] <- 2
      if(og == stn_loc) {
        state[state == 2] <- 1
        break
      }

      nextg <- river[og, 1]
      if(nextg == 0 | nextg == -1) {
        state[state == 2] <- -1
        break
      }
      if(nextg == 1) {
        state[state == 2] <- 1
        break
      }
      og <- nextg
    }
  }
  basin <- which(state > 0)

  # Calculate uh of each channel
  lcuh <- 48
  h <- sapply(1:ng, function(d) {
    d <- river[i, 2]
    if(d < 2500) { # For channels too short
      ih <- c(1., rep(0, lcuh))
    } else {
      t <- 0:lcuh * 3600
      ih <- d*exp(-(v[i]*t-d)**2/(4*diff[i]*t)) / (2*t*sqrt(pi*diff[i]*t))
      ih[1] <- 0
      if(any(ih > 0)) ih <- ih/sum(ih)
    }
    ih
  })

  # Calculate UHs
  uhlen <- 96
  uh <- matrix(0, ncol=length(state), nrow = uhlen*24)
  for(s in basin) {
    iuh <- rep(1/24, 24)
    chan <- which(river[,3] == s)
    while(TRUE) {
      iuh <- convolve(iuh, rev(h[,chan]), type = 'open')
      chan <- river[chan, 1]
      if(chan == 0)break
    }
    if(any(iuh > 0)) iuh <- iuh/sum(iuh)

    iuhlen <- min(uhlen*24, length(iuh))
    uh[1:iuhlen, s] <- iuh[1:iuhlen]
  }
  ind <- 0:uhlen * 24
  uh <- sapply(1:uhlen, function(i) colSums(uh[(ind[i]+1):ind[i+1], ]))
  uh <- t(uh)

  return(uh)
}
