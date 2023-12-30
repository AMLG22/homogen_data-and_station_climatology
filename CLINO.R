setwd('~/CLINO/')
#- CLINO.- Climatological Standard Normals calculation and CSV output files.
#Retired from the Spanish State Meteorological Agency (AEMET).
#Member of WMO Expert Team on Data Requirements for Climate services (ET-DRC).
#-- Version history. (Many thanks to Belinda Mariño, AEMET colleague, for
#   testing the different versions and finding some bugs in the code.):
#- This version fixes a small bug in the message of not enough data.
#- 1.5 completes the fix of version 1.4
#- 1.4 fixes a critical bug which misidentified the number of years with data
#  in supplementary variables, not writing them into the CSV files.
#- 1.3 forces reading station codes as class character to avoid
#  reading codes as e.g. 1234E as numeric in scientific format.
#- 1.2 extracts extreme values from original observations, to avoid
#  the report of extreme values different from real observations.
#- 1.1 corrected errors when data were provided for one station only,
#  a station had less than 24 years of data or the homogenized file had
#  supplementary stations, plus bad report of extreme daily values.
#  It also added the calculation of annual precipitation quintiles
#  and the date_format parameter to allow more flexibility.)
CLINO <- function(period='1981-2021',date_format='%Y-%m-%d') {
#period: period from which to compute the Climatological Standard Normals
#  (1991-2020 by default).
#-----------------
#This function needs four configuraton files (in CSV format, with header):
# CLINO_stations.csv : stations metadata
# CLINO_variables.csv : variables, file names and related parameters
# CLINO_parameters.csv : list of CLINO parameters (predefined and custom)
# CLINO_calculations.csv : calculation types and their CLINO codes
#-----------------  
#Daily homogenized data are provided in the files listed in CLINO_variables.csv
#They can be the *.rda results of Climatol daily homogenizations or CSV files
#  with stations daily data in columns, station codes in the header and a
#  first column with dates in YYYY-MM-DD format.
#----------------- Example:
#An example can be run by uncompressing CLINO_tool.tgz into our R working
#directory and running:   source('CLINO.R'); CLINO()
#----------------- Initial notice:
  cat('\nCLINO.- Climatological Standard Normals calculation and CSV output files.\n')
  cat('Version 1.5.1 (Jose A. Guijarro, 2023-02-16, under license GPL >= 3.\n')
  cat('Directions of use are available at:\n\n')
  cat('   http://www.climatol.eu/CLINO/CLINO-en.pdf (English)\n')
  cat('   http://www.climatol.eu/CLINO/CLINO-fr.pdf (French)\n')
  cat('   http://www.climatol.eu/CLINO/CLINO-es.pdf (Spanish)\n')
#----------------- Read input files:
  cat('\nReading input files...\n')
  cClas <- c("character","character","character","numeric","numeric",
    "numeric","character","character") #read stations metadata:
  stm <- read.csv('CLINO_stations.csv',as.is=TRUE,colClasses=cClas)
  nst <- nrow(stm) #no. of stations
  vrb <- read.csv('CLINO_variables.csv',as.is=TRUE) #variables, files and par.
  nv <- nrow(vrb) #no. of variables
  pmt <- read.csv('CLINO_parameters.csv',as.is=TRUE) #parameters list
  clc <- read.csv('CLINO_calculations.csv',as.is=TRUE) #calculations list
  clcn <- tolower(clc[,1]) #calculation names
#----------------- Calculate the CLINO values:
  cat('\nCalculation of the Climate Normals:\n\n')
  yrs <- unlist(strsplit(period,'-')) #first and last years of the CLINO period
  ny <- as.integer(yrs[2])-as.integer(yrs[1])+1 #no. of years
  dv <- seq(as.Date(sprintf('%s-01-01',yrs[1])),
    as.Date(sprintf('%s-12-31',yrs[2])),by='1 day') #target dates vector
  mv <- strftime(dv,'%m') #months vector
  yv <- strftime(dv,'%Y') #years vector
  dl <- split(dv,mv) #list of dates split by months
  cln <- array(NA,c(150,12,nst)) #initialization of CLINO values
  mnm <- array(NA,c(150,12,nst)) #initialization of monthly minimum values
  mnd <- array(NA,c(150,12,nst)) #initialization of monthly minimum dates
  mxm <- array(NA,c(150,12,nst)) #initialization of monthly maximum values
  mxd <- array(NA,c(150,12,nst)) #initialization of monthly maximum dates
  xnd <- array(NA,c(150,12,nst)) #initialization of daily extreme dates
  noy <- array(NA,c(150,12,nst)) #initialization of no. of years with data
  aqt <- matrix(NA,6,nst) #initialization of annual precipitation quintiles
  for(kv in 1:nv) { #for every variable:
    cat(sprintf('Processing %s data from file %s\n',vrb[kv,1],vrb[kv,2]))
    ddf <- vrb[kv,2] #daily data file
    z <- nchar(ddf); ext <- substring(ddf,z-2) #file extension
    if(ext=='rda') {
      load(ddf) #load climatol homogenization results (dat, dah, x, ...)
      stid <- est.c[1:nei,4] #station IDs
      dat <- as.matrix(dat[x%in%dv,]) #raw data of target dates
    } else if(ext=='csv') {
      dah <- read.csv(ddf,check.names=FALSE,as.is=TRUE)
      x <- as.Date(dah[,1],format=date_format) #vector of dates in the file
      if(ncol(dah)==2) cname=names(dah)[2] #avoid errors if only 1 station
      dah <- dah[,2:ncol(dah)] #get rid of the dates column
      stid <- names(dah); if(is.null(stid)) stid <- cname #station IDs
      dah <- as.matrix(dah) #convert to matrix
    }
    else error(sprintf('Extension of %s is neither "rda" nor "csv"',ddf))
    dah <- as.matrix(dah[x%in%dv,]) #keep only data of target dates
    #if raw data are unknown, use homogenized data instead:
    if(ext=='csv') dat <- as.matrix(dah)
    prm <- unlist(strsplit(vrb[kv,3],'-')) #parameters to calculate
    for(k in 1:nst) { #for every station:
      ke <- which(stid==stm[k,1])
      if(length(ke)==0) {
        cat(sprintf('  No %s data for station %s\n',vrb[kv,1],stm[k,1]))
        next
      }
      noyv <- noyr(dat[,ke],mv,ny) #no. of years with data of this variable
      for(i in 1:length(prm)) { #for every parameter:
        qp <- which(pmt[,1]==prm[i]) #row no. in the parameter list
        if(length(qp)==0) {
          cat(sprintf('  Parameter %s not found in CLINO_parameters.csv\n',
            prm[i]))
          next
        }
        #calculate monthly data:
        if(is.na(pmt[qp,5])) {
          if(pmt[qp,4]=='max'|pmt[qp,4]=='min') { #extremes from observed data:
            options(warn=-1) #suppress warnings due to months without data
            # which produce infinite values (which will be deleted below)
            md <- aggregate(dat[,ke],list(mv,yv),pmt[qp,4],na.rm=TRUE)
            options(warn=0) #restore normal warnings
        } else md <- aggregate(dah[,ke],list(mv,yv),pmt[qp,4])
        } else md <- aggregate(dah[,ke],list(mv,yv),pmt[qp,4],pmt[qp,5])
        md <- t(matrix(md[,3],12,ny))
        md[abs(md)==Inf] <- NA #delete spurious infinite values
        #calculate CLINO values:
        cln[qp,,k] <- round(apply(md,2,mean),1)
        mnm[qp,,k] <- round(apply(md,2,min,na.rm=TRUE),1) #MinMon
        mnd[qp,,k] <- as.integer(yrs[1])-1+apply(md,2,which.min) #years
        mxm[qp,,k] <- round(apply(md,2,max,na.rm=TRUE),1) #MaxMon
        mxd[qp,,k] <- as.integer(yrs[1])-1+apply(md,2,which.max) #years
        noy[qp,,k] <- noyv #no. of years with data
        #extreme dates (from original observations, rather than homogenized):
        if(pmt[qp,4]=='max'|pmt[qp,4]=='min') {
          z <- tapply(dat[,ke],mv,ifelse(pmt[qp,4]=='max',which.max,which.min))
          for(j in 1:12) xnd[qp,j,k] <- dl[[j]][z[j]]
        }
        if(qp==1) { #calculate monthly and annual precipitation quintiles:
          z <- round(apply(md,2,quantile,probs=c(0,.2,.4,.6,.8,1)),1)
          aqt[,k] <- round(quantile(apply(md,1,sum,na.rm=TRUE),
            probs=c(0,.2,.4,.6,.8,1)),1)
          for(zk in 1:6) {
            zq <- match(sprintf('Q%d',zk-1),pmt[,2])
            cln[zq,,k] <- z[zk,]
            noy[zq,,k] <- noy[1,,k] #copy precipitation no. of years
          }
        }
      }
    }
  }
#----------------- write the CLINO values into CSV files:
  cat('\nWriting the Climate Normals into CSV files:\n\n')
  for(k in 1:nst) { #for every station:
    #check if there are at least 80% (24 years) of observations:
    anoy <- apply(noy[,,k],1,min) #minimum annual no. of years
    if(max(anoy[1:8],na.rm=TRUE)<24) {
      cat(sprintf('Station %s %s has less than 24 years of data (skipped)\n',stm[k,2],stm[k,7]))
      next
    }
    stname <- iconv(stm[k,7], '', 'ascii', '')
    stname <- gsub('[^[:alnum:]]','',stname)
    fname <- sprintf('%s_%s.csv',stname,stm[k,2])
    cat('Writing file',fname,'\n')
    Fs <- file(fname,'w') #open output file
    #------- write file header:
    write(sprintf('World Meteorological Organization Climate Normals for %s,,,,,,,,,,,,,,,,',period),Fs)
    write('Single Station Data Sheet For All Climatological Surface Parameters,,,,,,,,,,,,,,,,',Fs)
    write(',,,,,,,,,,,,,,,,',Fs)
    write('Station Header Record,,,,,,,,,,,,,,,,',Fs)
    write(',,,,,,,,,,,,,,,,',Fs)
    write(sprintf('Country_Name,%s,,,,,,,,,,,,,,,',stm[k,8]),Fs)
    write(sprintf('Station_Name,%s,,,,,,,,,,,,,,,',stm[k,7]),Fs)
    write(',,,,,,,,,,,,,,,,',Fs)
    write('WMO_Number,Latitude,Longitude,Station_Height,,,,,,,,,,,,,,',Fs)
    lat <- deg2ddmmss(stm[k,4],'lat')
    lon <- deg2ddmmss(stm[k,5],'lon')
    write(sprintf('%s,%s,%s,%f,,,,,,,,,,,,,',stm[k,2],lat,lon,stm[k,6]),Fs)
    write(',,,,,,,,,,,,,,,,',Fs)
    write('WMO Integrated Global Observing System (WIGOS) Station Identifier (if available)',Fs)
    write(sprintf('%s,,,,,,,,,,,,,,,,',stm[k,3]),Fs)
    write(',,,,,,,,,,,,,,,,',Fs)
    write(',,,,,,,,,,,,,,,,',Fs)
    #------- write values of calculated parameters:
    write('Principal Climatological Surface Parameters,,,,,,,,,,,,,,,,',Fs)
    for(kp in 1:nrow(pmt)) { #for every parameter:
      qp <- pmt[kp,1] #parameter code
      if(qp>99) break #end of predefined parameters
      if(kp==9) {
        write(',,,,,,,,,,,,,,,,',Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        write('Secondary and Other Climatological Surface Parameters (add as needed),,,,,,,,,,,,,,,,',Fs)
      }
      if(qp==11) { #precipitation quintiles go together:
        i <- match(11,pmt[,1])
        Kp <- match(110:115,pmt[,1])
        write(',,,,,,,,,,,,,,,,',Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        write('Parameter_Code,Parameter_Name,Units,,,,,,,,,,,,,,',Fs)
        write(sprintf('%d,%s,%s,,,,,,,,,,,,,,',pmt[i,1],pmt[i,2],pmt[i,3]),Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        write('WMO_Number,Parameter_Code,Calculation_Name,Calculation_Code,January,February,March,April,May,June,July,August,September,October,November,December,Annual',Fs)
        for(iq in Kp) {
          z <- sprintf('%s,%d,%s,%s',stm[k,2],pmt[i,1],clc[iq-75,1],
            clc[iq-75,2])
          for(j in 1:12) z <- paste(z,sprintf('%.1f',cln[iq,j,k]),sep=',')
          kaqt <- match(pmt[iq,2],c('Q0','Q1','Q2','Q3','Q4','Q5'))
          z <- paste(z,',',sprintf('%.1f',aqt[kaqt,k]),sep='')
          write(paste(z,',',sep=''),Fs)
        }
        z <- sprintf('%s,11,NOY,98',stm[k,2])
        write(paste(z,paste(noy[1,,k],collapse=','),',',sep=','),Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        next
      }
      if(is.na(cln[kp,1,k])) { #check for multiple instances of the parameter:
        z <- qp*10+0:9 #possible codes of the parameter instances
        Kp <- match(z,pmt[,1])
        if(sum(!is.na(Kp))==0) Kp <- kp #no instances found; keep kp
        else Kp <- Kp[!is.na(Kp)] #list of found instances
      } else Kp <- kp
      for(i in Kp) {
        if(is.na(anoy[i])) next #no data for this parameter
        if(anoy[i]<24) next #less than 24 years of observations
        write(',,,,,,,,,,,,,,,,',Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        write('Parameter_Code,Parameter_Name,Units,,,,,,,,,,,,,,',Fs)
        write(sprintf('%d,%s,%s,,,,,,,,,,,,,,',pmt[i,1],pmt[i,2],pmt[i,3]),Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        write('WMO_Number,Parameter_Code,Calculation_Name,Calculation_Code,January,February,March,April,May,June,July,August,September,October,November,December,Annual',Fs)
        qc <- match(pmt[i,4],clcn) #calculation line
        if(is.na(qc)) qc <- match(pmt[i,3],clcn)
        z1 <- sprintf('%s,%d',stm[k,2],pmt[i,1]) #WMO no. & parameter code
        z <- sprintf('%s,%s,%s',z1,clc[qc,1],clc[qc,2])
        if(pmt[i,6]=='max') zv <- mxm[i,,k]
        else if(pmt[i,6]=='min') zv <- mnm[i,,k]
        else zv <- cln[i,,k]
        for(j in 1:12) z <- paste(z,sprintf('%.1f',zv[j]),sep=',')
        if(is.na(pmt[i,6])) z <- paste(z,',',sep='')
        else z <- paste(z,sprintf('%.1f',round(eval(call(pmt[i,6],zv))
          ,1)),sep=',')
        write(z,Fs)
        if(!is.na(xnd[i,j,k])) { #dates of extreme values:
          if(pmt[i,4]=='max') {
            z <- sprintf('%s,MaxDate,15',z1)
            jx <- which.max(cln[i,,k]) #month of maximum value
          } else {
            z <- sprintf('%s,MinDate,16',z1)
            jx <- which.min(cln[i,,k]) #month of minimum value
          }
          for(j in 1:12) z <- paste(z,sprintf('%s',
            strftime(as.Date(xnd[i,j,k],origin='1970-01-01'),'%Y/%d')),sep=',')
          paste(z,sprintf('%s',strftime(as.Date(xnd[i,jx,k],
            origin='1970-01-01'),'%Y/%m/%d')),sep=',')
          write(z,Fs)
        } else if(!is.na(mnm[i,j,k])) { #maximum and minimum monthly values:
          z <- sprintf('%s,MinMon,17',z1)
          for(j in 1:12) z <- paste(z,sprintf('%.1f',mnm[i,j,k]),sep=',')
          write(paste(z,',',sep=''),Fs)
          z <- sprintf('%s,%d,DMinMon,18',stm[k,2],pmt[i,1])
          for(j in 1:12) z <- paste(z,sprintf('%d',mnd[i,j,k]),sep=',')
          write(paste(z,',',sep=''),Fs)
          z <- sprintf('%s,%d,MaxMon,19',stm[k,2],pmt[i,1])
          for(j in 1:12) z <- paste(z,sprintf('%.1f',mxm[i,j,k]),sep=',')
          write(paste(z,',',sep=''),Fs)
          z <- sprintf('%s,%d,DMaxMon,20',stm[k,2],pmt[i,1])
          for(j in 1:12) z <- paste(z,sprintf('%d',mxd[i,j,k]),sep=',')
          write(paste(z,',',sep=''),Fs)
        }
        z <- sprintf('%s,%d,NOY,98',stm[k,2],pmt[i,1])
        write(paste(z,paste(noy[i,,k],collapse=','),',',sep=','),Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
        write(',,,,,,,,,,,,,,,,',Fs)
      }
    }
    close(Fs)
  }
  cat('\n')
}

#===== Auxiliary functions: ===============================

#- deg2ddmmss.- Convert degrees with decimals to degrees, minutes and seconds.
deg2ddmmss <- function(deg,coord) {
#deg: degrees with decimals
#coord: either 'lat' or 'lon', according with the deg coordinate
  if(deg<0) { 
    deg <- -deg
    if(coord=='lat') L <- 'S' else L <- 'W'
  }
  else { if(coord=='lat') L <- 'N' else L <- 'E' }
  dd <- floor(deg); deg <- (deg - dd) * 60
  mm <- floor(deg); deg <- (deg - mm) * 60
  ss <- round(deg)
  if(ss==60) {
    ss <- 0; mm <- mm+1
    if(mm==60) { mm <- 0; deg <- deg+1 }
  }
  if(coord=='lat') return(sprintf('%02d|%02d|%02d|%s',dd,mm,ss,L))
  else return(sprintf('%03d|%02d|%02d|%s',dd,mm,ss,L))
}

#- greq.- Count no. of values greater or equal to a threshold
greq <- function(dat,thr) return(sum(dat >= thr, na.rm=TRUE))

#- grth.- Count no. of values greater than a threshold
grth <- function(dat,thr) return(sum(dat > thr, na.rm=TRUE))

#- loeq.- Count no. of values lower or equal to a threshold
loeq <- function(dat,thr) return(sum(dat <= thr, na.rm=TRUE))

#- loth.- Count no. of values lower than a threshold
loth <- function(dat,thr) return(sum(dat < thr, na.rm=TRUE))

#- noyr.- Rounded no. of years of original data
noyr <- function(d,mv,ny) {
  nas <- aggregate(d,list(mv),is.na) #missing data (by month)
  nm <- unlist(lapply(nas$x,sum)) #monthly no. of missing data
  nd <- unlist(lapply(nas$x,length)) #monthly no. of days in the period
  return(round(ny*(1-nm/nd)))
}
