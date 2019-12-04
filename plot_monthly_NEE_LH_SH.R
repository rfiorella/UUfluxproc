# plot_monthly_NEE_and_LH.R
# rich fiorella 190429
# updated 190919 to include sensible heat
# and add a 7pt median filter.

# plots monthly NE, LH, SH data from NEON h5 files.

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
sh.max.threshold <- 1400
sh.min.threshold <- -500

# get a list of available files.
flist <- list.files("/uufs/chpc.utah.edu/common/home/bowen-group2/NEON/raw/DP4_00200_001",pattern='.h5',full.names=TRUE,recursive=TRUE)

# open csv of NEON site metadata.
neon_sites <- read_csv("../../metadata/190115-field-sites.csv")

#------------------------------------------------------------------
# extract list of stations that script will run through.
sites.tmp <- strsplit(flist,split=".",fixed=TRUE)   # split file names to extract part of file name w/ site ID
codes <- unique(sapply(sites.tmp,'[[',5))           # extract and unlist unique site IDs
names <- sapply(codes,
                function(x){neon_sites$`Site Name`[neon_sites$`Site ID` == x]}) # match site names to list of unique site IDs.
domain.number <- sapply(codes,
                        function(x){neon_sites$`Domain Number`[neon_sites$`Site ID` == x]}) # get domain numbers
domain.name <- sapply(codes,
                      function(x){neon_sites$`Domain Name`[neon_sites$`Site ID` == x]}) # get domain numbers

# create a vector of variables to *ALWAYS* keep through the loop
vars.to.keep <- c("vars.to.keep","flist","neon_sites","sites.tmp","codes","names","domain.number","domain.name")
#------------------------------------------------------------------
# loop through sites, extract CO2 and H2O calibration data.
for (i in 1:length(codes)) {
  
  # open PDF
  pdf(paste0("/uufs/chpc.utah.edu/common/home/bowen-group2/NEON/plots/1910xx/monthly/",codes[i],"_NEE_LH.pdf"))
  
  # print site code for tracking
  print(codes[i])

  # return list of files associated only w/ this site
  site.files <- list.files(paste0("/uufs/chpc.utah.edu/common/home/bowen-group2/NEON/raw/DP4_00200_001/",codes[i]),pattern='.h5',full.names=TRUE,recursive=TRUE)

  # loop through site files, make sure heights don't change
  #--------------------------------------------------------
  attrs <- h5readAttributes(site.files[1],codes[i])
  
  # extract heights from attrs list.
  heights <- as.numeric(attrs$DistZaxsLvlMeasTow)
  
  # get site month/year for plot title.
  slist.tmp <- strsplit(site.files,split="/",fixed=TRUE)
  yrmn <- sapply(slist.tmp,'[[',11)
  rm(slist.tmp)
 
  print(site.files)
  print(yrmn) 
  # extract carbon isotope data.
  for (j in 1:length(site.files)) {
    
    # print j for tracking
    print(j)

    # get data.
    tmp.nee.data <- h5read(site.files[j],paste0('/',codes[i],'/dp04/data/fluxCo2'))
    tmp.lh.data <- h5read(site.files[j],paste0('/',codes[i],'/dp04/data/fluxH2o'))
    tmp.sh.data <- h5read(site.files[j],paste0('/',codes[i],'/dp04/data/fluxTemp'))
    
    # get NEE timeseries.
    nee <- tmp.nee.data[['nsae']] # nee$flux is the flux term in units umol CO2 m-2 s-1
    lh <- tmp.lh.data[['nsae']]
    sh <- tmp.sh.data[['nsae']]
    
    # contains NaN - replace w/ NA.
    nee$flux[is.na(nee$flux)] <- NA
    lh$flux[is.na(lh$flux)] <- NA
    sh$flux[is.na(sh$flux)] <- NA
    
   if (do.filtering == TRUE) {
      # set to NA outside of thresholds.
      nee$flux[nee$flux>nee.max.threshold] <- NA
      nee$flux[nee$flux<nee.min.threshold] <- NA
	
      lh$flux[lh$flux>lh.max.threshold] <- NA
      lh$flux[lh$flux<lh.min.threshold] <- NA
	
      sh$flux[sh$flux>sh.max.threshold] <- NA
      sh$flux[sh$flux<sh.min.threshold] <- NA
   }

    #---------------------------------------------------
    # convert to xts, sort by time, plot.
    #---------------------------------------------------
    # first, the times.
    t1 <- as.POSIXct(nee$timeBgn,format="%Y-%m-%dT%H:%M:%S.%OSZ",
                     tz="UTC",origin="1970-01-01")
    t2 <- as.POSIXct(nee$timeEnd,format="%Y-%m-%dT%H:%M:%S.%OSZ",
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
    
    nee.plot <- xts(nee$flux,order.by=t.mean)
    lh.plot <- xts(lh$flux,order.by=t.mean)
    sh.plot <- xts(sh$flux,order.by=t.mean)
    
    # filter out large differences, if do.filtering is TRUE.
    if (do.filtering) {
       # apply median filter.
       nee.plot <- rollapply(nee.plot,7,median,na.rm=TRUE,fill=NA)
       lh.plot <- rollapply(lh.plot,7,median,na.rm=TRUE,fill=NA)
       sh.plot <- rollapply(sh.plot,7,median,na.rm=TRUE,fill=NA)
    }
    #---------------------------------------------------------
    # set xts object names.
    names(nee.plot) <- "NEE"
    names(lh.plot)  <- "LH"
    names(sh.plot)  <- "SH"
    
    p1 <- ggplot(data=nee.plot,aes(x=index(nee.plot),y=NEE)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("date") +
      scale_y_continuous("NEE (umol CO2 m-2 s-1)") +
      ggtitle(paste(yrmn[j],codes[i],names[i],domain.number[i],domain.name[i]))
    
    #print(p1)
    
    p2 <- ggplot(data=lh.plot,aes(x=index(lh.plot),y=LH)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("date") +
      scale_y_continuous("LH (W m-2)") +
      ggtitle(paste(yrmn[j],codes[i],names[i],domain.number[i],domain.name[i]))
    
    p3 <- ggplot(data=sh.plot,aes(x=index(sh.plot),y=SH)) +
      geom_line() + 
      theme_bw() +
      scale_x_datetime("date") +
      scale_y_continuous("SH (W m-2)") +
      ggtitle(paste(yrmn[j],codes[i],names[i],domain.number[i],domain.name[i]))
    
    grid.arrange(p1,p2,p3,nrow=3)
  } # j
  
  dev.off() # close pdf for that site and...
  
  # once again, wishing on the HDF5 monkey paw.
  h5closeAll()
  
} # close i - iterate to next site...
    
    
    
