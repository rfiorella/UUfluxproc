# package_NEON_for_REddyProc.R
# rich fiorella 191127

# remove existing data?
rm(list=ls())

# load packages
library(rhdf5)
library(neonUtilities)
library(REddyProc)
library(xts)
library(lubridate)
library(tidyverse)

# need to gather and preprocess data.
fluxes <- stackEddy("~/Dropbox/NEON/DP4_00200_001/NOGP/",level="dp04")

nee <- fluxes$NOGP$data.fluxCo2.nsae.flux
lhf <- fluxes$NOGP$data.fluxH2o.nsae.flux
shf <- fluxes$NOGP$data.fluxTemp.nsae.flux
ustar <- fluxes$NOGP$data.fluxMome.turb.veloFric
flux.time <- fluxes$NOGP$timeBgn

flux.time <- as.POSIXct(flux.time,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")

rm(fluxes)

# filter out insane values.
nee[nee>40] <- nee[nee< -40] <- NA

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
all.data <- all.data["2018-01-01/2019-12-31"]

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

#------------------------------------------------------------
#------------------------------------------------------------
# now go through REddyProc

rm(list=ls())

EddyData <- fLoadTXTIntoDataframe("NOGP_forREddyProc.txt",Dir="~/Desktop")
EddyData <- filterLongRuns(EddyData,"NEE")
# convert VPD to numeric.
EddyData$VPD <- as.numeric(EddyData$VPD)

#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- fConvertTimeToPosix(
  EddyData, 'YDH',Year = 'Year',Day = 'DoY', Hour = 'Hour') %>% 
  filterLongRuns("NEE")
#+++ Initalize R5 reference class sEddyProc for post-processing of eddy data
#+++ with the variables needed for post-processing later
EProc <- sEddyProc$new(
  'US-Wref', EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))

EProc$sPlotFingerprintY('NEE', Year = 2019)

EProc$sEstimateUstarScenarios(
  nSample = 100L, probs = c(0.05, 0.5, 0.95))

EProc$sGetUstarScenarios()

EProc$sMDSGapFillUStarScens('NEE')

#EProc$sPlotFingerprintY('NEE_U50_f', Year = 2018)

# partition (nighttime method.)
EProc$sSetLocationInfo(LatDeg = 46.77, LongDeg = -100.92, TimeZoneHour = 0)
EProc$sMDSGapFill('Tair',FillAll = FALSE, minNWarnRunLength = NA)
EProc$sMDSGapFill('VPD', FillAll = FALSE, minNWarnRunLength = NA)

EProc$sMRFluxPartitionUStarScens()

# plot GPP and Reco
yr.to.print <- 2018

pdf("NOGP_2018.pdf",width=12,height=8)
par(mfrow=c(1,3))
par(mar=c(1,1,0,0),oma=c(1,1,1,1))
EProc$sPlotFingerprintY('NEE_U50_f', Year = yr.to.print)
EProc$sPlotFingerprintY('GPP_U50_f', Year = yr.to.print)
EProc$sPlotFingerprintY('Reco_U50', Year = yr.to.print)
dev.off()

yr2.to.print <- 2019

pdf("NOGP_2019.pdf",width=12,height=8)
par(mfrow=c(1,3))
par(mar=c(1,1,0,0),oma=c(1,1,1,1))
EProc$sPlotFingerprintY('NEE_U50_f', Year = yr2.to.print)
EProc$sPlotFingerprintY('GPP_U50_f', Year = yr2.to.print)
EProc$sPlotFingerprintY('Reco_U50', Year = yr2.to.print)
dev.off()

# pdf("testlegends.pdf",width=12,height=3)
# par(mfrow=c(1,3))
# par(mar=c(0.1,0.1,0.1,0.1),oma=c(0.1,0.1,0.1,0.1))
# EProc$sPlotFingerprintY('NEE_U50_f', Year = yr.to.print, onlyLegend = TRUE)
# EProc$sPlotFingerprintY('GPP_U50_f', Year = yr.to.print, onlyLegend = TRUE)
# EProc$sPlotFingerprintY('Reco_U50', Year = yr.to.print, onlyLegend = TRUE)
# dev.off()

