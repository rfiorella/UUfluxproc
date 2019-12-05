# package_NEON_for_REddyProc.R


extract_NEON_fluxes_REddyProc <- function(site,year=9999,flux.path,met.path,
                                          nee.max,nee.min,lh.max,lh.min) {
  
  # list required packages
  require(rhdf5)
  require(neonUtilities)
  require(xts)
  require(lubridate)
  require(tidyverse)
  
  # stack flux data.
  fluxes <- stackEddy(paste0(flux.path,"/",site,"/"),level="dp04")
  
  fluxes.flat <- fluxes[[site]] # flatten list structure.

  # extract required variables.
  fluxes.reduced <- fluxes.flat %>%
    select(timeBgn,timeEnd, # time variables
           data.fluxCo2.nsae.flux,data.fluxCo2.stor.flux,data.fluxCo2.turb.flux, # CO2 fluxes
           data.fluxH2o.nsae.flux,data.fluxH2o.stor.flux,data.fluxH2o.turb.flux, # H2O fluxes,
           data.fluxTemp.nsae.flux,data.fluxTemp.nsae.flux,data.fluxTemp.turb.flux, # Temp fluxes,
           data.fluxMome.turb.veloFric) # u*
  
  # convert "NEON" time to POSIXct
  fluxes.reduced$timeBgn <- as.POSIXct(fluxes.reduced$timeBgn,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")
  
  # filter out insane values?
  # nee[nee>40] <- nee[nee< -40] <- NA

}
# create xts 
flux.vars <- data.frame(nee,lhf,shf,ustar)
flux.xts <- xts(flux.vars,order.by=flux.time)

# median filter fluxes?
flux.xts <- rollapply(flux.xts,7,median,na.rm=TRUE)
#------------------------------------------------------------
# get Rh, Rg, and Ta data.
Rh.tmp <- loadByProduct("DP1.00098.001",site="NOGP",startdate="2017-01",avg=30,check.size=F)
Rg.tmp <- loadByProduct("DP1.00023.001",site="NOGP",startdate="2017-01",avg=30,check.size=F)
Ta.tmp <- loadByProduct("DP1.00003.001",site="NOGP",startdate="2017-01",avg=30,check.size=F)

# pull out met variables.
Rh <- Rh.tmp$RH_30min %>%
  filter(horizontalPosition == 0) %>%
  select(RHMean,startDateTime) 
Rg <- Rg.tmp$SLRNR_30min %>%
  filter(horizontalPosition == 0) %>%
  select(inSWMean,startDateTime) 
Ta <- Ta.tmp$TAAT_30min %>%
  select(tempTripleMean,startDateTime)

names(Rh) <- c("Rh","time")
names(Rg) <- c("Rg","time")
names(Ta) <- c("Ta","time")

# clean up.
rm(Rh.tmp,Rg.tmp,Ta.tmp)

# change time vars to character, then to posix ct.
Rg$time <- as.POSIXct(Rg$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
Rh$time <- as.POSIXct(Rh$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
Ta$time <- as.POSIXct(Ta$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")

# convert data to xts objects and then merge.
Rg.xts <- xts(Rg$Rg,order.by=Rg$time)
Rh.xts <- xts(Rh$Rh,order.by=Rh$time)
Ta.xts <- xts(Ta$Ta,order.by=Ta$time)

names(Rg.xts) <- "Rg"
names(Rh.xts) <- "Rh"
names(Ta.xts) <- "Ta"

minTime <- min(c(min(index(Rg.xts)),min(index(Rh.xts)),min(index(Ta.xts)),min(index(flux.xts))))
maxTime <- max(c(max(index(Rg.xts)),max(index(Rh.xts)),max(index(Ta.xts)),max(index(flux.xts))))

dummy.ts <- seq.POSIXt(minTime,maxTime,by=1800)
dummy.data <- rep(NA,length(dummy.ts))

dummy.xts <- xts(dummy.data,order.by=dummy.ts)

# create bound xts
all.data <- merge(dummy.xts,Rg.xts,Rh.xts,Ta.xts,flux.xts)

# cut to min/max of flux data.

# convert to MPI required format.
out.data <- coredata(all.data)
time.tmp <- index(all.data)

# pull out time variables as required.
yr <- year(time.tmp)
doy <- yday(time.tmp)
hr <- as.numeric(hour(time.tmp)+minute(time.tmp)/60)

head1 <- c("Year","DoY","Hour","NEE","LH","H","Rg","Tair","Tsoil","rH","VPD","Ustar")
head2 <- c("--","--","--","umol-2-s","Wm-2","Wm-2","Wm-2","degC","degC","%","hPa","ms-1")
data.out <- data.frame(yr,doy,hr,out.data[,"nee"],
                       out.data[,"lhf"],out.data[,"shf"],out.data[,"Rg"],
                       out.data[,"Ta"],rep(NA,nrow(out.data)),
                       out.data[,"Rh"],rep(NA,nrow(out.data)),out.data[,"ustar"])

data.out <- do.call(rbind,list(head1,head2,data.out))

write.table(data.out,"~/Desktop/NOGP_forREddyProc.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)



