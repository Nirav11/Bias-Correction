rm(list = ls())

library(CRHMr)
library(Reanalysis)
library(raster)
library(ncdf4)
library(stringr)


setwd('D:/Langtang/ERA5_Langtang/')

# Creating SSRD Obs files for different CSV
ssrd_path <- './SSR'
ssrd_list = list.files(ssrd_path,pattern="nc",full.names=T)


# print files in short wave list
for  (i in 1:length(ssrd_list)){
  ssrd <- ERAgetNearestTimeseries(
    ssrd_list[i],
    'ssrd',
    85.56169, #For Kyanjing Station
    28.21081,
    projection = "+proj=utm +zone=45 +ellps=WGS84",
    timezone = "Etc/GMT+1",
    quiet = TRUE,
    logfile = "test_SSRD")
  writeObsFile(ssrd,file.path(ssrd_path,paste0('Kyanjing_',str_sub(ssrd_list[i], start = 7, end = 22),'.obs')))
    
  
}

ssrd_list = list.files(ssrd_path,pattern="obs",full.names=T)
final_ssrd <- createObsDataframe('1980-12-31','1981-01-01', timestep = 1, 
                                 variables = c('Qsi'), 
                                 reps = 1, timezone = 'Etc/GMT+1')

for  (i in 1:length(ssrd_list)){
  ssrd <- readObsFile(ssrd_list[i],timezone = 'Etc/GMT+1')
  names(ssrd)[c(2)] <- c('Qsi')
  ssrd[c(2)] <- ssrd[c(2)] / 3600
  final_ssrd <- appendObs(final_ssrd,ssrd)
}

#-----------------------------------------------

# Creating tsr Obs files for different CSV
tsr_path <- './TSR'
tsr_list = list.files(tsr_path,pattern="nc",full.names=T)


# print files in short wave list
for  (i in 1:length(tsr_list)){
  tsr <- ERAgetNearestTimeseries(
    tsr_list[i],
    'strd',
    85.56169, #For Kyanjing Station
    28.21081,
    projection = "+proj=utm +zone=45 +ellps=WGS84",
    timezone = "Etc/GMT+1",
    quiet = TRUE,
    logfile = "test_tsr")
  writeObsFile(tsr,file.path(tsr_path,paste0('Kyanjing_',str_sub(tsr_list[i], start = 7, end = 22),'.obs')))
  
}

tsr_list = list.files(tsr_path,pattern="obs",full.names=T)
final_tsr <- createObsDataframe('1980-12-31','1981-01-01', timestep = 1, 
                                 variables = c('Qli'), 
                                 reps = 1, timezone = 'Etc/GMT+1')

for  (i in 1:length(tsr_list)){
  tsr <- readObsFile(tsr_list[i],timezone = 'Etc/GMT+1')
  names(tsr)[c(2)] <- c('Qli')
  tsr[c(2)] <- tsr[c(2)] / 3600
  final_tsr <- appendObs(final_tsr,tsr)
}

#------------------------------------------------

# Creating temp Obs files for different CSV
temp_path <- './Temp'
temp_list = list.files(temp_path,pattern=".nc",full.names=T)


# print files in short wave list
for  (i in 1:length(temp_list)){
  temp <- ERAgetNearestTimeseries(
    temp_list[i],
    't2m',
    85.56169, #For Kyanjing Station
    28.21081,
    projection = "+proj=utm +zone=45 +ellps=WGS84",
    timezone = "Etc/GMT+1",
    quiet = TRUE,
    logfile = "test_temp")
  writeObsFile(temp,file.path(temp_path,paste0('Kyanjing_',str_sub(temp_list[i], start = 8, end = 23),'.obs')))
  
}

temp_list = list.files(temp_path,pattern="obs",full.names=T)
final_temp <- createObsDataframe('1980-12-31','1981-01-01', timestep = 1, 
                                variables = c('t'), 
                                reps = 1, timezone = 'Etc/GMT+1')

for  (i in 1:length(temp_list)){
  temp <- readObsFile(temp_list[i],timezone = 'Etc/GMT+1')
  names(temp)[c(2)] <- c('t')
  final_temp <- appendObs(final_temp,temp)
}

names(final_temp)[c(2)] <- c('t.1')
#------------------------------------------------------
# print files for Dew Point Temperature list

dew_path <- './Dew'
dew_list = list.files(dew_path,pattern=".nc",full.names=T)
# print files in Dew list
for  (i in 1:length(dew_list)){
  dew <- ERAgetNearestTimeseries(
    dew_list[i],
    'd2m',
    85.56169, #For Kyanjing Station
    28.21081,
    projection = "+proj=utm +zone=45 +ellps=WGS84",
    timezone = "Etc/GMT+1",
    quiet = TRUE,
    logfile = "test_dew")
  writeObsFile(dew,file.path(dew_path,paste0('Kyanjing_',str_sub(dew_list[i], start = 10, end = 25),'.obs')))
  
}

dew_list = list.files(dew_path,pattern="obs",full.names=T)
final_dew <- createObsDataframe('1980-12-31','1981-01-01', timestep = 1, 
                                 variables = c('d2m'), 
                                 reps = 1, timezone = 'Etc/GMT+1')

for  (i in 1:length(dew_list)){
  dew <- readObsFile(dew_list[i],timezone = 'Etc/GMT+1')
  names(dew)[c(2)] <- c('d2m')
  final_dew <- appendObs(final_dew,dew)
}


#------------------------------------------------------

uwind_path <- './UWind'
uwind_list = list.files(uwind_path,pattern=".nc",full.names=T)

# print files in uwind list
for  (i in 1:length(uwind_list)){
  uwind <- ERAgetNearestTimeseries(
    uwind_list[i],
    'u10',
    85.56169, #For Kyanjing Station
    28.21081,
    projection = "+proj=utm +zone=45 +ellps=WGS84",
    timezone = "Etc/GMT+1",
    quiet = TRUE,
    logfile = "test_uwind")
  writeObsFile(uwind,file.path(uwind_path,paste0('Kyanjing_',str_sub(uwind_list[i], start = 9, end = 26),'.obs')))
  
}

uwind_list = list.files(uwind_path,pattern="obs",full.names=T)
final_uwind <- createObsDataframe('1980-12-31','1981-01-01', timestep = 1, 
                                variables = c('u10'), 
                                reps = 1, timezone = 'Etc/GMT+1')

for  (i in 1:length(uwind_list)){
  uwind <- readObsFile(uwind_list[i],timezone = 'Etc/GMT+1')
  names(uwind)[c(2)] <- c('u10')
  final_uwind <- appendObs(final_uwind,uwind)
}


#------------------------------------------------------

vwind_path <- './VWind'
vwind_list = list.files(vwind_path,pattern=".nc",full.names=T)

# print files in vwind list
for  (i in 1:length(vwind_list)){
  vwind <- ERAgetNearestTimeseries(
    vwind_list[i],
    'v10',
    85.56169, #For Kyanjing Station
    28.21081,
    projection = "+proj=utm +zone=45 +ellps=WGS84",
    timezone = "Etc/GMT+1",
    quiet = TRUE,
    logfile = "test_vwind")
  writeObsFile(vwind,file.path(vwind_path,paste0('Kyanjing_',str_sub(vwind_list[i], start = 9, end = 26),'.obs')))
  
}

vwind_list = list.files(vwind_path,pattern="obs",full.names=T)
final_vwind <- createObsDataframe('1980-12-31','1981-01-01', timestep = 1, 
                                  variables = c('v10'), 
                                  reps = 1, timezone = 'Etc/GMT+1')

for  (i in 1:length(vwind_list)){
  vwind <- readObsFile(vwind_list[i],timezone = 'Etc/GMT+1')
  names(vwind)[c(2)] <- c('v10')
  final_vwind <- appendObs(final_vwind,vwind)
}

final_wind <- final_uwind
final_wind[2] <- sqrt((final_uwind[2]*final_uwind[2])+(final_vwind[2]*final_vwind[2]))
names(final_wind)[c(2)] <- c('u.1')

# ----------------------------------------------------------

# For Precipitation

precip_path <- './Precip/'
precip_list = list.files(precip_path,pattern=".nc",full.names=T)

# print files in precip list
for  (i in 1:length(precip_list)){
  precip <- ERAgetNearestTimeseries(
    precip_list[i],
    'tp',
    85.56169, #For Kyanjing Station
    28.21081,
    projection = "+proj=utm +zone=45 +ellps=WGS84",
    timezone = "Etc/GMT+1",
    quiet = TRUE,
    logfile = "test_precip")
  writeObsFile(precip,file.path(precip_path,paste0('Kyanjing_',str_sub(precip_list[i], start = 10, end = 26),'.obs')))
  
}

precip_list = list.files(precip_path,pattern="obs",full.names=T)
final_precip <- createObsDataframe('1980-12-31','1981-01-01', timestep = 1, 
                                  variables = c('tp'), 
                                  reps = 1, timezone = 'Etc/GMT+1')

for  (i in 1:length(precip_list)){
  precip <- readObsFile(precip_list[i],timezone = 'Etc/GMT+1')
  names(precip)[c(2)] <- c('tp')
  final_precip <- appendObs(final_precip,precip)
}

final_precip$tp.1 <- final_precip$tp.1*1000
names(final_precip)[c(2)] <- c('p.1')


deaccum <- ERA5deaccum(final_precip, colnum = 1)
deaccum$tp.1[deaccum$tp.1 < 0] <- 0

# Temporary variable for storing temp and dew point temp in o Celcius

dew_temp <- final_dew
dew_temp[2] <- dew_temp[2] - 273.15

temp_temp <- final_temp
temp_temp[2] <- temp_temp[2] - 273.15
names(temp_temp)[c(2)] <- c('t.1')

RH_Era <- dew_temp
RH_Era[2] <- 0
RH_Era[2] <- 100 * (exp(17.625 * dew_temp[2]/(243.04+dew_temp[2]))/exp(17.625 * temp_temp[2]/(243.04+temp_temp[2])))
names(RH_Era)[c(2)] <- c('rh.1')

dd <- 0
dd <- merge(temp_temp,RH_Era)
dd <- merge(dd, final_ssrd)
dd <- merge(dd, final_tsr)
dd <- merge(dd, final_wind)
dd <- merge(dd, deaccum)

writeObsFile(dd, file.path('../For Bias 11 Dec 2022/Forcing Input/ERA5/','Kyanjing_ERA5_1981_2020.obs'), comment = "Unbiased ERA5 for temp, rh, Qsi, Qli, u and p")

ERA5 <- readObsFile('../For Bias 11 Dec 2022/Forcing Input/ERA5/Kyanjing_ERA5_1981_2020.obs', timezone = 'Etc/GMT+1')

ea <- changeRHtoEa(ERA5)
writeObsFile(ea,file.path('../For Bias 11 Dec 2022/Forcing Input/ERA5/','Kyanjing_ERA5_1981_2020_ea.obs'))



obs.file <- read.csv('../For Bias 11 Dec 2022/Kyanjing_2012_2019.csv')
obs.file$datetime <- as.POSIXct(paste(obs.file$DATE, obs.file$TIME), format="%m/%d/%Y %H:%M", tz = "UTC")
obs <- subset(obs.file, select = c(datetime, TAIR, RH,KINC,LINC,WSPD,PVOL))
names(obs)[c(2:7)] <- c('t.1','rh.1','Qsi.1','Qli','u.1','p.1')
writeObsFile(obs,'../For Bias 11 Dec 2022/Forcing Input/Observed/Kyanjing_AWS_2012_2019.obs')


# ------------------------------------------------------------
merra_path <- '../For Bias 11 Dec 2022/Forcing Input/MERRA2/'

merra_list = list.files(merra_path,pattern=".csv",full.names=T)
# print files in merra list
df <- 0

merra_first <- read.csv(merra_list[1], header=TRUE)[1]     # gene names
rf    <- do.call(cbind,lapply(merra_list,function(fn)read.csv(fn, header=TRUE)[2]))
rf    <- cbind(merra_first,rf)
rf$datetime <- as.POSIXct(rf$datetime,format = "%m/%d/%Y %H:%M", tz = "Etc/GMT+1")
names(rf)[c(2:8)] <- c('Qli.1','p.1','Qsi.1','t.1', 'd.1','ut.1','vt.1')
rf$u.1 <- sqrt((rf$ut.1*rf$ut.1)+(rf$vt.1*rf$vt.1)) 
rf$rh.1 <- 100 * (exp(17.625 * rf$d.1/(243.04+rf$d.1))/exp(17.625 * rf$t.1/(243.04+rf$t.1)))
rf <- subset(rf, select = c(datetime, t.1, rh.1, Qsi.1, Qli.1, u.1, p.1))

writeObsFile(rf, file.path(merra_path,'Merra2trhQsiQliup.obs'))


#100 × {exp[17.625 × Dp/(243.04 + Dp)]/exp[17.625 × T/(243.04 + T)]}
# nc <- nc_open('../ERA5_Langtang/Precip/lg_prcp_1981_1990.nc')
# time <- ncvar_get(nc, "time")
# dim_time <- length(time)
# precip <- ncvar_get(nc, "tp", start = c(1, 1, 1), count = c(1, 1, dim_time))
# nc_close(nc)
# 
# start <- 100
# end <- 124
# df <- data.frame(time[start:end], precip[start:end])
# names(df) <- c("hour", "total_precip")
# plot(df)
# 
# 
# df <- data.frame(time, precip)
# names(df) <- c("datetime", "precip")
# df$datetime <- as.POSIXct(df$datetime*3600, origin= "1900-01-01 00:00", tz = "GMT")
# deaccum <- ERA5deaccum(df, colnum = 1)
# deaccum$total_precip <- df$precip
# 
# p <- plot(deaccum$datetime[start:end], deaccum$precip[start:end], type = "p", col = "red")
# 
# 
