# plot_monthly_PMinv_vars.R
# rich fiorella 190505

# plots monthly data required to perform the PM inversion:
# LH, air T, wind speed, radiation, soil heat flux, VPD, air density
rm(list=ls())

library(rhdf5)
library(tidyverse)
library(xts)
library(gridExtra)

#-------------------------------------------------
do.filtering <- TRUE
nee.max.threshold <- 50
nee.min.threshold <- -50
lh.max.threshold <- 500
lh.min.threshold <- -200

# get a list of available files.
flist <- list.files("~/Desktop/NEON/",pattern='.h5',full.names=TRUE,recursive=TRUE)

# open csv of NEON site metadata.
neon_sites <- read_csv("../../metadata/190115-field-sites.csv")

#------------------------------------------------------------------
# extract list of stations that script will run through.
sites.tmp <- strsplit(flist,split=".",fixed=TRUE)   # split file names to extract part of file name w/ site ID
codes <- unique(sapply(sites.tmp,'[[',3))           # extract and unlist unique site IDs
names <- sapply(codes,
                function(x){neon_sites$`Site Name`[neon_sites$`Site ID` == x]}) # match site names to list of unique site IDs.
domain.number <- sapply(codes,
                        function(x){neon_sites$`Domain Number`[neon_sites$`Site ID` == x]}) # get domain numbers
domain.name <- sapply(codes,
                      function(x){neon_sites$`Domain Name`[neon_sites$`Site ID` == x]}) # get domain numbers


# clean up
rm(neon_sites)
rm(sites.tmp)

# create a vector of variables to *ALWAYS* keep through the loop
vars.to.keep <- c("vars.to.keep","flist","codes","names","domain.number","domain.name")
#------------------------------------------------------------------
# loop through sites, extract CO2 and H2O calibration data.
for (i in 1:length(codes)) {
  
  # print site code for tracking
  print(codes[i])
  
  neon.path = "~/Desktop/NEON/"
  
  # get all the files containing the relevant data.
  flux.files <- list.files(paste0(neon.path,codes[i]),pattern='.h5',full.names=TRUE,recursive=TRUE) # EC files
  rh.files <- list.files(paste0(neon.path,codes[i]),pattern='DP1.00098.001.000',full.names=TRUE,recursive=TRUE) # RH files
  bp.files <- list.files(paste0(neon.path,codes[i]),pattern='DP1.00004.001.000',full.names=TRUE,recursive=TRUE) # pressure
  rnt.files <- list.files(paste0(neon.path,codes[i]),pattern='DP1.00023.001.000',full.names=TRUE,recursive=TRUE) # radiation
  shf.files <- list.files(paste0(neon.path,codes[i]),pattern='DP1.00040.001.0',full.names=TRUE,recursive=TRUE) # soil heat flux
  airT.files <- list.files(paste0(neon.path,codes[i]),pattern='DP1.00003.001.000',full.names=TRUE,recursive=TRUE) # triple asp. air T
  wspd.files <- list.files(paste0(neon.path,codes[i]),pattern='DP1.00001.001.000',full.names=TRUE,recursive=TRUE) # wind speed.

  # determine how many levels are present at this site.
  #--------------------------------------------------------
  attrs <- h5readAttributes(flux.files[1],codes[i])
  
  # extract heights from attrs list.
  #--------------------------------------------------------
  heights <- as.numeric(attrs$DistZaxsLvlMeasTow)
  rm(attrs)

  # get list of available months for this site
  #--------------------------------------------------------
  extract_sitemons <- function(x) {
    # break files names into list, pull out the part that is the yr/mn combo.
    
    pieces <- strsplit(x,split="/",fixed=TRUE)
    yrmn <- sapply(pieces,'[[',7) # NOTE: THIS LINE CHANGES IF DIRECTORY STRUCTURE CHANGES!!!
    
    return(yrmn)
  }
  
  # return list of months for each data file.
  flux.yrmn <- extract_sitemons(flux.files)
  rh.yrmn <- extract_sitemons(rh.files)
  bp.yrmn <- extract_sitemons(bp.files)
  rnt.yrmn <- extract_sitemons(rnt.files)
  shf.yrmn <- extract_sitemons(shf.files)
  airT.yrmn <- extract_sitemons(airT.files)
  wspd.yrmn <- extract_sitemons(wspd.files)
  
  # open PDF
  pdf(paste0("~/Dropbox/NEON_Plots/190401/monthly/PMinv_",codes[i],".pdf"))

  # extract carbon isotope data.
  for (j in 1:length(flux.files)) {
    
    # print j for tracking
    print(j)

    #----------------------------------------
    # Get, arrange all monthly data.
    #----------------------------------------
    #========================================
    # flux files.
    tmp.lh.data <- h5read(flux.files[j],paste0('/',codes[i],'/dp04/data/fluxH2o'))

    # once again, wishing on the HDF5 monkey paw.
    h5closeAll()
    
    # get LH timeseries.
    lh <- tmp.lh.data[['nsae']] # units are W/m2
    rm(tmp.lh.data)
    
    # contains NaN - replace w/ NA.
    lh$flux[is.na(lh$flux)] <- NA
    
    # filter, if do.filtering is TRUE.
    if (do.filtering) {
      # set to NA outside of thresholds.
      lh$flux[lh$flux>lh.max.threshold] <- NA
      lh$flux[lh$flux<lh.min.threshold] <- NA
    }
    
    #---------------------------------------------------
    # convert to xts, sort by time, plot.
    #---------------------------------------------------
    # first, the times.
    t1 <- as.POSIXct(lh$timeBgn,format="%Y-%m-%dT%H:%M:%S.%OSZ",
                     tz="UTC",origin="1970-01-01")
    t2 <- as.POSIXct(lh$timeEnd,format="%Y-%m-%dT%H:%M:%S.%OSZ",
                     tz="UTC",origin="1970-01-01")
    
    # weird R tirefire here that I need to convert to numeric before
    # taking mean time.
    t1 <- as.numeric(unlist(t1))
    t2 <- as.numeric(unlist(t2))
    
    tdf <- data.frame(t1,t2)
    
    # take mean
    t.mean <- rowMeans(tdf)
    
    # convert back to times.
    t.mean <- as.POSIXct(t.mean,origin="1970-01-01",tz="UTC")

    # convert to xts
    lh.plot <- xts(lh$flux,order.by=t.mean)
    
    # set xts object names.
    names(lh.plot)  <- "LH"  
    
    #---------------------------------------------------------
    # remaining tower met variables need a helper function that extracts
    # data from: 1) highest position on tower, and 2) only from files 
    # in that month.
    
    select_metfile <- function(x,x.yearmonth,flux.yearmonth) {
      
      # check to see which files in x are in the specified yearmonth
      # yearmonth should be a string in YYYY-MM format. flux indicates
      # it comes from EC flux file list, x.yearmonth corresponds to 
      # var x.
      
      in.month <- x[x.yearmonth == flux.yearmonth]
      
      # if in.month returns only 1 file = we're done.
      # if in.month returns null = there are no data for this var form that month.
      # if in.month returns >1 file, need to select among different heights on the
      # tower (most likely) *SOIL HEAT FLUX is an exception, this is not measured
      # on the tower!!!
      
      if (length(in.month) == 0) {
        
        return(NULL) # must be a better way?
        
      } else if (length(in.month) == 1) {
        
        # we're done, just return in.month
        return(in.month)
    
      } else {
    
        # if we're interested in value highest up on tower, then 
        # take the highest true value.
        return(tail(in.month,n=1))
        
      }
      
    }
    
    #=========================================================
    # net radiation
    
    # find which files correspond to yrmn of flux file.
    rnt.files.in.month <- select_metfile(rnt.files,rnt.yrmn,flux.yrmn[j])

    # open file.
    options(readr.num_columns = 0)
    
    if (!is.null(rnt.files.in.month)) {
      rnt.data <- read_csv(rnt.files.in.month)
      
      # cut down to only the variables needed to calculate Rnet
      rnt.data <- rnt.data %>%
        select(startDateTime,endDateTime,inSWMean,outSWMean,inLWMean,outLWMean) %>%
        drop_na()  %>%# get rid of rows w/ missing vals.
        mutate(Rnet = (inSWMean - outSWMean) + (inLWMean - outLWMean))
      
      # create xts object.
      rnt.data <- rnt.data %>%
        mutate(meanTime = mean(c(startDateTime,endDateTime)))
      
      # create xts object
      rnt.plot <- xts(rnt.data$Rnet,order.by=rnt.data$startDateTime) # NOTE! FLUX VAR USES MEAN TIME!!! OFFSET BY 15 min.
      names(rnt.plot) <- "Rnet"
    } else {
      rnt.plot <- xts(order.by=t.mean)
      rnt.plot <- merge(rnt.plot,Rnet=NA)
    }
    #=========================================================
    # soil heat flux.
    
    # find which files correspond to yrmn of flux file.
    shf.files.in.month <- select_metfile(shf.files,shf.yrmn,flux.yrmn[j])
    
    # open file.
    if (!is.null(shf.files.in.month)) {
      shf.data <- read_csv(shf.files.in.month,
                           col_type=cols(SHFMean=col_double(),
                                         SHFMinimum=col_double(),
                                         SHFMaximum=col_double(),
                                         SHFVariance=col_double(),
                                         SHFNumPts=col_double(),
                                         SHFExpUncert=col_double()))
      
      # cut down to only the variables needed to calculate Rnet
      shf.data <- shf.data %>%
        select(startDateTime,endDateTime,SHFMean) %>%
        drop_na() # get rid of rows w/ missing vals.
      
      # create xts object
      shf.plot <- xts(shf.data$SHFMean,order.by=shf.data$startDateTime) # NOTE! FLUX VAR USES MEAN TIME!!! OFFSET BY 15 min.
      names(shf.plot) <- "SHF"
    } else {
      shf.plot <- xts(order.by=t.mean)
      shf.plot <- merge(shf.plot,SHF=NA)
    }
    
    #=========================================================
    # air temperature
    
    # find which files correspond to yrmn of flux file.
    airT.files.in.month <- select_metfile(airT.files,airT.yrmn,flux.yrmn[j])
    
    # open file.
    if (!is.null(airT.files.in.month)) {
      airT.data <- read_csv(airT.files.in.month)
      
      # cut down to only the variables needed to calculate Rnet
      airT.data <- airT.data %>%
        select(startDateTime,endDateTime,tempTripleMean) %>%
        drop_na() # get rid of rows w/ missing vals.
      
      # create xts object
      airT.plot <- xts(airT.data$tempTripleMean,order.by=airT.data$startDateTime) # NOTE! FLUX VAR USES MEAN TIME!!! OFFSET BY 15 min.
      names(airT.plot) <- "AirT"
    } else {
      airT.plot <- xts(order.by=t.mean)
      airT.plot <- merge(airT.plot,AirT=NA)
    }
    
    
    #=========================================================
    # relative humidity - use along w/ Air T to get VPD.
    
    # find which files correspond to yrmn of flux file.
    rh.files.in.month <- select_metfile(rh.files,rh.yrmn,flux.yrmn[j])
    
    # open file.
    if (!is.null(airT.files.in.month) & !is.null(rh.files.in.month)) {
      rh.data <- read_csv(rh.files.in.month)
      
      # cut down to only the variables needed to calculate Rnet
      rh.data <- rh.data %>%
        select(startDateTime,endDateTime,RHMean) %>%
        drop_na() # get rid of rows w/ missing vals.
      
      # use airT to get saturation vapor pressure.
      vpd.data <- full_join(rh.data,airT.data,by="startDateTime")
      
      vpd.data <- vpd.data %>%
        mutate(TK = tempTripleMean + 273.15) %>% # convert to Kelvin
        mutate(e.sat = ifelse(TK>273.15,
                              exp(53.67957 - 6743.769/TK - 4.8451*log(TK)), # vapor over liquid saturation
                              exp(23.33086 - 6111.72784/TK + 0.15215*log(TK)))) %>% # vapor over ice saturation
        mutate(e = e.sat*(RHMean/100)) %>%
        mutate(VPD = (e.sat - e)/10) # dividing by 10 here converts hPa to kPa.
      
      # create xts object
      vpd.plot <- xts(vpd.data$VPD,order.by=vpd.data$startDateTime) # NOTE! FLUX VAR USES MEAN TIME!!! OFFSET BY 15 min.
      names(vpd.plot) <- "VPD"
    } else {
      vpd.plot <- xts(order.by=t.mean)
      vpd.plot <- merge(vpd.plot,VPD=NA)
    }
    
    #--------------------------------------------------------------------
    # Barometric pressure.
    
    # find which files correspond to yrmn of flux file.
    bp.files.in.month <- select_metfile(bp.files,bp.yrmn,flux.yrmn[j])
    
    # open file.
    if (!is.null(bp.files.in.month) & !is.null(rh.files.in.month) & !is.null(airT.files.in.month)) {
      bp.data <- read_csv(bp.files.in.month)
      
      # cut down to only the variables needed to calculate Rnet
      bp.data <- bp.data %>%
        select(startDateTime,endDateTime,staPresMean) %>%
        mutate(staPresMean = 10*staPresMean) %>% # BaroPressure given in kPa - converting to more typical hPa (at least for atm folks...)
        drop_na() # get rid of rows w/ missing vals.
      
      # convert pressure to air density.
      dens.data <- full_join(vpd.data,bp.data,by="startDateTime")
      
      dens.data <- dens.data %>%
        mutate(density = 100*(staPresMean - e)/(287.058*TK)) %>% # density = p (in Pa)/RT
        filter(density > 0.5 & density < 2) # some crazy values...why? 
      
      # create xts object
      dens.plot <- xts(dens.data$density,order.by=dens.data$startDateTime) # NOTE! FLUX VAR USES MEAN TIME!!! OFFSET BY 15 min.
      names(dens.plot) <- "Rho"
    } else {
      dens.plot <- xts(order.by=t.mean)
      dens.plot <- merge(dens.plot,Rho=NA)
    }
    
    #--------------------------------------------------------------------
    # wind speed.
    
    # find which files correspond to yrmn of flux file.
    wspd.files.in.month <- select_metfile(wspd.files,wspd.yrmn,flux.yrmn[j])

    # open file.
    if (!is.null(wspd.files.in.month)) {
      wspd.data <- read_csv(wspd.files.in.month)
      
      # cut down to only the variables needed to calculate Rnet
      wspd.data <- wspd.data %>%
        select(startDateTime,endDateTime,windSpeedMean) %>%
        drop_na() # get rid of rows w/ missing vals.
      
      # create xts object
      wspd.plot <- xts(wspd.data$windSpeedMean,order.by=wspd.data$startDateTime) # NOTE! FLUX VAR USES MEAN TIME!!! OFFSET BY 15 min.
      names(wspd.plot) <- "WSpd"
    } else {
      wspd.plot <- xts(order.by=t.mean)
      wspd.plot <- merge(wspd.plot,WSpd=NA)
    }
    
    #=====================================================================
    # Make plots
    
    # set common minimum and maximum time.
    min.time <- min(c(index(lh.plot),
                      index(rnt.plot),
                      index(shf.plot),
                      index(airT.plot),
                      index(vpd.plot),
                      index(dens.plot),
                      index(wspd.plot)))
    
    max.time <- max(c(index(lh.plot),
                      index(rnt.plot),
                      index(shf.plot),
                      index(airT.plot),
                      index(vpd.plot),
                      index(dens.plot),
                      index(wspd.plot)))
    
    # LH plot.
    p1 <- ggplot(data=lh.plot,aes(x=index(lh.plot),y=LH)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("",limits=c(min.time,max.time),expand = c(0,0)) +
      scale_y_continuous("LH")# + # W/m2
      #ggtitle(paste(flux.yrmn[j],codes[i],names[i],domain.number[i],domain.name[i]))
    
    # make plot of radiation.
    p2 <- ggplot(data=rnt.plot,aes(x=index(rnt.plot),y=Rnet)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("",limits=c(min.time,max.time),expand = c(0,0)) +
      scale_y_continuous("Rnet") # W/m2
    
    # make plot of soil heat flux
    p3 <- ggplot(data=shf.plot,aes(x=index(shf.plot),y=SHF)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("",limits=c(min.time,max.time),expand = c(0,0)) +
      scale_y_continuous("Soil HF") # W/m2
    
    # air temperature plot
    p4 <- ggplot(data=airT.plot,aes(x=index(airT.plot),y=AirT)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("",limits=c(min.time,max.time),expand = c(0,0)) +
      scale_y_continuous("T") # °C 
  
    # vpd plot
    p5 <- ggplot(data=vpd.plot,aes(x=index(vpd.plot),y=VPD)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("",limits=c(min.time,max.time),expand = c(0,0)) +
      scale_y_continuous("VPD") # kPa 
  
    # barometric pressure
    p6 <- ggplot(data=dens.plot,aes(x=index(dens.plot),y=Rho)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("",limits=c(min.time,max.time),expand = c(0,0)) +
      scale_y_continuous("Rho") # kg/m3
      
    # wind speed
    p7 <- ggplot(data=wspd.plot,aes(x=index(wspd.plot),y=WSpd)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("date",limits=c(min.time,max.time),expand = c(0,0)) +
      scale_y_continuous("Wind") # m/s
    
    grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=7,top=paste(flux.yrmn[j],codes[i],names[i],domain.number[i],domain.name[i]))
  } # j
  
  dev.off() # close pdf for that site and...

  
} # close i - iterate to next site...
    
    
    