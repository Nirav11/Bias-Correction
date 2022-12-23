# 21 Dec 2022: qmap doQmapRQUANT tricub wet.day=FALSE qstep=0.1,nboot=10
# bias-correct reanalysis data with observation at Kyanjing Station
# MERRA2, ERA5, WFDEI

rm(list=ls())
library(CRHMr)
library(ggplot2)
library(hydroGOF) #for rmse
library(sirad) #modeval
library(qmap)
library(raster)
library(ncdf4)
library(Reanalysis)


#define working directory
setwd('D:/Langtang/For Bias 11 Dec 2022/')
outPath <- file.path('./Output')

source('./paul_Wang-Bovik.R') #paul's function to get WBI


#function for stat
#----------------------------------------------------------------
dostat <- function(test, para, what, observed, simulated){
  stat <- modeval(simulated, observed, minlength = 30) # if below 30 datasets, produce NAs
  wb <- Wang_Bovik(simulated, observed)
  stat0 <- data.frame(test, para, what,
                      N=stat$N, pearson=stat$pearson, MBE=stat$MBE, 
                      MAE=stat$MAE, RMSE=stat$RMSE,
                      r2= stat$R2, slope=stat$slope, intercept = stat$intercept, 
                      WBI=wb[1], Mxy=wb[2], Vxy=wb[3], Rxy=wb[4])
  return(stat0)
}

d_stat <- data.frame() #(stat)

#function to make legend transparent
#----------------------------------------------------------------
transparent_legend =  theme(legend.background = element_rect(fill = "transparent"),
                            legend.key = element_rect(fill = "transparent", color = "transparent"))
#----------------------------------------------------------------

#Path for input forcing files
watch <- file.path('./Forcing Input/Watch')
era5 <- file.path('./Forcing Input/ERA5')
merra2 <- file.path('./Forcing Input/MERRA2')


#first read all the datasets 
###############################################
# A. Reanalysis - MERRA2
MERRA2 <- readObsFile(file.path(merra2, 'Merra2trhQsiQliup.obs'),
                    timezone = 'Etc/GMT+1')

# MERRA2 <- subset(MERRA2, select = c(datetime, t.1, rh.1, u.1, Qsi.1, Qli.1, p.1))
MERRA2 <- changeRHtoEa(MERRA2)

# B. Reanalysis - ERA5
ERA5 <- readObsFile(file.path(era5,'Kyanjing_ERA5_1981_2020.obs'),
                   timezone = 'Etc/GMT+1')
names(ERA5)[7] <- 'p.1'
ERA5 <- changeRHtoEa(ERA5)

#ERA5 <- subset(ERA5, select = c(datetime, t.1, ea.1, u.1, Qsi.1, Qli.1, p.1))

# C. Reanalysis - WFDEI
WFDEI <- readObsFile(file.path(watch,'Langtang_WFDEI_teauQsiQlip_1979_2018.obs'),
                     timezone = 'Etc/GMT+1')

# D. Obs data
peyto <- readObsFile('Kyanjing_AWS_2012_2019.obs', 
                   timezone = 'Etc/GMT+1')
peyto <- subset(peyto, select = c(datetime, t.1, rh.1, Qsi.1, Qli.1, u.1, p.1))
peyto <- changeRHtoEa(peyto)
peyto <- peyto[complete.cases(peyto),]


###############################################

biasfunc <- function(insitu = pyeto, re, test){
  # re <- subset(re, datetime >= db1 & datetime <= db4)
  df <- merge(insitu, re, by='datetime', suffixes = c('.obs', '.re'))
  df$month <- as.numeric(format(df$datetime, format = '%m'))
  re$month <- as.numeric(format(re$datetime, format = '%m'))

  bias.corrected <- createObsDataframe('1981-01-01', '2020-12-31', timestep = 1,
                                  variables = c('t', 'ea', 'Qsi', 'Qli', 'u', 'p'),
                                  reps = 1, timezone = 'Etc/GMT+1')
  bias.corrected$month <- as.numeric(format(bias.corrected$datetime, format = '%m'))

  for(j in 1:6){ #j <- 1 #for t, ea, u, Qsi, Qli, p
  for(i in 1:12){ #i <- 2 #for 12 months
    obs<-df[,j+1][which(df$month==i)]
    print(obs)
    mod<-df[,j+7][which(df$month==i)]
    print(mod)
    qm.fit <- fitQmapRQUANT(obs,mod,qstep=0.1,nboot=10,wet.day=FALSE)
    # print(qm.fit)
    bias.corrected[,j+1][bias.corrected$month==i] <- doQmapRQUANT(re[,j+1][re$month==i],
                                                         qm.fit,type="tricub")
    }
  }
  return(bias.corrected)
}

#apply the biasfunc and get bias corrected variables
MERRA2.bias <- biasfunc(insitu = peyto, re=MERRA2, test = 'MERRA2')
ERA5.bias <- biasfunc(insitu = peyto, re=ERA5, test='ERA5')
WFDEI.bias <- biasfunc(insitu = peyto, re=WFDEI, test = 'WFDEI')


MERRA2.bias.rh <- changeEatoRH(MERRA2.bias)
writeObsFile(MERRA2.bias.rh[-8], obsfile = file.path(outPath, 'MERRA2_tRHuQsiQlip_BiasCorrect2PeytoBow_Jan1979_Dec2019.obs'),
             comment = 'Bias corrected to Peyto Main for teauQsiQli, and to BowSummit for p')

ERA5.bias.rh <- changeEatoRH(ERA5.bias)
writeObsFile(ERA5.bias.rh[-8], obsfile = file.path(outPath, 'ERA5_tRHQsiQliu_BiasCorrectKyanjing_1981_2020.obs'),
             comment = 'Bias corrected to Kyanjing for teaQsiQliu')

WFDEI.bias.rh <- changeEatoRH(WFDEI.bias)
writeObsFile(WFDEI.bias.rh[-8], obsfile = file.path(outPath, 'WFDEI_tRHuQsiQlip_BiasCorrect2PeytoBow_Jan1979_Dec2018.obs'),
             comment = 'Bias corrected to Peyto Main for teauQsiQli, and to BowSummit for p')


#----------------------------------------------------------------


# DO STAT
#######################

# first, merge obs data and bias corrected data
statfunc <- function(insitu, bias.corrected, test){
  obs.re <- merge(insitu, bias.corrected, by='datetime', suffixes = c('.obs', '.re'))
  what='bias_correction'
  for(i in 1:6){ # i = 6;
    observed <- obs.re[i+1]
    print(observed)
    simulated <- obs.re[i+7]
    para <- stringr::str_split(names(observed), pattern = '[.]')[[1]][1]
    names(observed)[1] <- 'var'
    names(simulated)[1] <- 'var'
    observed <- observed$var
    simulated <- simulated$var
    stat0 <- dostat(test, para, what, observed, simulated)
    d_stat <- rbind(d_stat, stat0)
  }
  return(d_stat)
}

#---------------------------------------------------------------

MERRA2.stat <- statfunc(insitu=peyto, bias.corrected=MERRA2.bias, test='MERRA2')
ERA5.stat <- statfunc(insitu=peyto, bias.corrected=ERA5.bias, test='ERA5')
WFDEI.stat <- statfunc(insitu=peyto, bias.corrected=WFDEI.bias, test='WFDEI')

statALL <- rbind(MERRA2.stat, ERA5.stat, WFDEI.stat)
statALL[5:15] <- round(statALL[5:15],2)
write.csv(statALL, file.path(outPath, 'Bias_correct_stat.csv'), row.names = F)

#---------------------------------------------------------------

OBSMERRA2 <- merge(peyto, MERRA2.bias, by='datetime', suffixes = c('.obs', '.MERRA2'))
ERA5WFDEI <- merge(ERA5.bias, WFDEI.bias, by='datetime', suffixes = c('.ERA5', '.WFDEI'))
ALL <- merge(OBSMERRA2, ERA5WFDEI, by='datetime')


# get ppt
######
ALL.ppt <- subset(ALL, select = c(datetime, p.1.obs, p.1.MERRA2, p.1.ERA5, p.1.WFDEI))
# ALL.ppt.CFSR <- subset(ALL_CFSR, select = c(datetime, p.1.obs, p.1.CFSR))
# monthly aggregation to get stat
ALL.ppt.m <- aggDataframe(ALL.ppt, columns = c(1:4), period = 'monthly', funs = c('sum'))

# Cumulative precip for plotting
ALL.ppt[2:5] <- cumsum(ALL.ppt[2:5])

######

#plot cumulative
###########
g <- ggplot(ALL.ppt, aes(x=datetime)) +
  geom_line(aes(y=p.1.obs, col='Observation'))+
  geom_line(aes(y=p.1.MERRA2, col='MERRA2'))+
  geom_line(aes(y=p.1.ERA5, col='ERA5'))+
  geom_line(aes(y=p.1.WFDEI, col='WFDEI')) + theme_bw() +
  ylab('Cumulative precipitation [mm]') + xlab('Year') +
  scale_color_manual(name='', values = c('Observation'='red', 'MERRA2'='purple',
                                         'ERA5'='blue', 'WFDEI'='green'),
                     guide='legend')+transparent_legend+
  theme(legend.position=c(.25,.80))
ggsave(g, filename = file.path(outPath, 'cumulative_precip.png'),
       width = 6, height = 5)

#########

##get stat from monthly aggregated values for precip
#######################
p_stat <- NULL  
para='precip.monthly'
what="Bias.correct"
  for(i in 3:5){ # i = 3; 
    observed <- ALL.ppt.m[2]
    simulated <- ALL.ppt.m[i]
    test <- stringr::str_split(names(simulated), pattern = '[.]')[[1]][3]  
    names(observed)[1] <- 'var'
    names(simulated)[1] <- 'var'
    observed <- observed$var
    simulated <- simulated$var
    stat0 <- dostat(test, para,what, observed, simulated)
    p_stat <- rbind(p_stat, stat0)
  }
  p_stat[5:15] <- round(p_stat[5:15],2)
  write.csv(p_stat, file.path(outPath, 'Bias_correct_p_stat.csv'), row.names = F)
#######################
  
  
#======================================================================================
######  GGPLOT          Function to Plot xy and label  MBE RMSE R2   ##############
#======================================================================================
ggxy <- function(data, obs, sim, xlab, ylab, title,xpos=3/5, ypos=1/5, outplot){
  #need to determine limits for plot
  xymax <- max(max(obs), max(sim))
  xymin <- min(min(obs), min(sim))
  int <- abs(xymax-xymin)/25
  #stat
  stat <- modeval(sim, obs, minlength = 30) # if below 30 datasets, produce NAs
  wb <- Wang_Bovik(sim, obs)
  statALL <- data.frame(N=stat$N, MBE=stat$MBE, MAE=stat$MAE, RMSE=stat$RMSE, 
                      r2= stat$R2, slope=stat$slope, intercept=stat$intercept,
                      WBI=wb[1], Mxy=wb[2], Vxy=wb[3], Rxy=wb[4])

  textlab2 <- paste("WBI = ", round(statALL$Wang.Bovik.Index,2))
  textlab3 <- paste("RMSE = ", round(statALL$RMSE,2))
  textlab4 <- paste("MBE = ", round(statALL$MBE,2))
  textlab5 <- paste('N = ', statALL$N)
  plot <- ggplot(data,aes(obs,sim), environment=environment())+
    xlab(xlab)+ylab(ylab)+ggtitle(title)+
    geom_point(shape=1)+ # Use hollow circles
    geom_abline()+ #(intercept=0,slope=1)
    theme_bw(10)+
    geom_smooth(method=lm, colour="blue")+
    coord_equal(ratio = 1) +
    xlim(c(xymin, xymax))+ylim(c(xymin, xymax))+
    # annotate("text", x = xpos*xymax, y = xymin+ypos*(xymax+abs(xymin))-0*int, label = textlab1, color="darkblue", size = 4.5, parse=FALSE)+
    annotate("text", x = xpos*xymax, y = xymin+ypos*(xymax+abs(xymin))-1*int, label = textlab2, color="darkblue", size = 3.1, parse=FALSE)+
    annotate("text", x = xpos*xymax, y = xymin+ypos*(xymax+abs(xymin))-2*int, label = textlab3, color="darkblue", size = 3.1, parse=FALSE)+
    annotate("text", x = xpos*xymax, y = xymin+ypos*(xymax+abs(xymin))-3*int, label = textlab4, color="darkblue", size = 3.1, parse=FALSE)+
    annotate("text", x = xpos*xymax, y = xymin+ypos*(xymax+abs(xymin))-4*int, label = textlab5, color="darkblue", size = 3.1, parse=FALSE)
  # return(plot)
  ggsave(plot, file= outplot, width = 5, height = 5)
}
#######################################################################################
#define obseveration and simulation data, before plot

  # xlab <- 'test'
  # ylab <- 'test'
  # title <- 'test'
  # data <- OBSMERRA2
  # obs <- data$t.1.obs
  # sim <- data$t.1.MERRA2
  

#temp 
##########
# HA. TempERA5ture (y=1/3 for annotate) 
# HA1. WFDEI
outplot=file.path(outPath, "Temp_WFDEI.png")
ggxy(data=ALL, obs=ALL$t.1.obs, sim=ALL$t.1.WFDEI,
     xlab=expression(paste('Observation (',degree, 'C)')), # [�C]
     ylab=expression(paste('WFDEI (',degree, 'C)')),
     title = "Hourly Temperature",
     xpos = 3/5, ypos = 1/5, outplot)

# HA2. ERA5
outplot=file.path(outPath, "Temp_ERA5.png")
ggxy(data=ALL, obs=ALL$t.1.obs, sim=ALL$t.1.ERA5, 
     xlab=expression(paste('Observation (',degree, 'C)')), # [�C]
     ylab=expression(paste('ERA5 (',degree, 'C)')),
     title = "Hourly Temperature ", 
     xpos = 3/5, ypos = 1/5, outplot)

# # HA3. MERRA2
outplot=file.path(outPath, "Temp_MERRA2.png")
ggxy(data=ALL, obs=ALL$t.1.obs, sim=ALL$t.1.MERRA2,
     xlab=expression(paste('Observation (',degree, 'C)')), # [�C]
     ylab=expression(paste('MERRA2 (',degree, 'C)')),
     title = "Hourly Temperature",
     xpos = 3/5, ypos = 1/5, outplot)


##########

#vap
##########
# HB. Vapour Pressure  
# HB1. WFDEI
outplot=file.path(outPath, "Vap_WFDEI.png")
ggxy(data=ALL, obs=ALL$ea.1.obs, sim=ALL$ea.1.WFDEI, 
     xlab='Observation (kPa)', 
     ylab='WFDEI (kPa)',
     title = "Hourly Vapour Pressure", 
     xpos = 4/5, ypos = 1/5, outplot)

# HB2. ERA5
outplot=file.path(outPath, "Vap_ERA5.png")
ggxy(data=ALL, obs=ALL$ea.1.obs, sim=ALL$ea.1.ERA5, 
     xlab='Observation (kPa)', 
     ylab='ERA5 (kPa)',
     title = "Hourly Vapour Pressure", 
     xpos = 4/5, ypos = 1/5, outplot)

# HB3. MERRA2
outplot=file.path(outPath, "Vap_MERRA2.png")
ggxy(data=ALL, obs=ALL$ea.1.obs, sim=ALL$ea.1.MERRA2, 
     xlab='Observation (kPa)', 
     ylab='MERRA2 (kPa)',
     title = "Hourly Vapour Pressure", 
     xpos = 4/5, ypos = 1/5, outplot)

# HB3. MERRA2
outplot=file.path(outPath, "Vap_MERRA2.png")
ggxy(data=ALL, obs=ALL$ea.1.obs, sim=ALL$ea.1.MERRA2, 
     xlab='Observation (kPa)', 
     ylab='MERRA2 (kPa)',
     title = "Hourly Vapour Pressure", 
     xpos = 4/5, ypos = 1/5, outplot)
##########

#wind
##########
# HC. Wind Speed (y=1/5, colour = 'darkblue' in annotate)
# HC1. WFDEI
outplot=file.path(outPath, "wind_WFDEI_qmap.png")
ggxy(data=ALL, obs=ALL$u.1.obs, sim=ALL$u.1.WFDEI, 
     xlab= expression(paste('Observation (m',s^-1,')',sep = '')),
     ylab= expression(paste('WFDEI (m',s^-1,')',sep = '')),
     title = "Hourly Wind Speed", 
     xpos = 2/15, ypos = 14/15, outplot)

# HC2. ERA5
outplot=file.path(outPath, "wind_ERA5.png")
ggxy(data=ALL, obs=ALL$u.1.obs, sim=ALL$u.1.ERA5, 
     xlab= expression(paste('Observation (m',s^-1,')',sep = '')),
     ylab= expression(paste('ERA5 (m',s^-1,')',sep = '')),
     title = "Hourly Wind Speed", 
     xpos = 2/15, ypos = 14/15, outplot)

# HC3. MERRA2
outplot=file.path(outPath, "wind_MERRA2.png")
ggxy(data=ALL, obs=ALL$u.1.obs, sim=ALL$u.1.MERRA2, 
     xlab= expression(paste('Observation (m',s^-1,')',sep = '')),
     ylab= expression(paste('MERRA2 (m',s^-1,')',sep = '')),
     title = "Hourly Wind Speed", 
     xpos = 2/15, ypos = 14/15, outplot)

##########

#Inc.SW
##########
# HD. Inc. SW (y=1/5, colour = 'darkblue' in annotate)
# HD1. WFDEI
outplot=file.path(outPath, "Qsi_WFDEI.png")
ggxy(data=ALL, obs=ALL$Qsi.1.obs, sim=ALL$Qsi.1.WFDEI, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('WFDEI (W',m^-2,')',sep = '')),
     title = "Hourly Incoming Shortwave Radiation", 
     xpos = 6.3/7, ypos = 1.3/7, outplot)

# HD2. ERA5
outplot=file.path(outPath, "Qsi_ERA5.png")
ggxy(data=ALL, obs=ALL$Qsi.1.obs, sim=ALL$Qsi.1.ERA5, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('ERA5 (W',m^-2,')',sep = '')),
     title = "Hourly Incoming Shortwave Radiation", 
     xpos = 6.3/7, ypos = 1.3/7, outplot)

# HD3. MERRA2
outplot=file.path(outPath, "Qsi_MERRA2.png")
ggxy(data=ALL, obs=ALL$Qsi.1.obs, sim=ALL$Qsi.1.MERRA2, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('MERRA2 (W',m^-2,')',sep = '')),
     title = "Hourly Incoming Shortwave Radiation", 
     xpos = 6.3/7, ypos = 1.3/7, outplot)
##########

#Inc. LW
##########
# HE. Inc. LW (y=1/5, colour = 'darkblue' in annotate)
# HE1. WFDEI
outplot=file.path(outPath, "Qli_WFDEI.png")
ggxy(data=ALL, obs=ALL$Qli.1.obs, sim=ALL$Qli.1.WFDEI, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('WFDEI (W',m^-2,')',sep = '')),
     title = "Hourly Incoming Longwave Radiation", 
     xpos = 6.3/7, ypos = 1.3/7, outplot)

# HE2. ERA5
outplot=file.path(outPath, "Qli_ERA5.png")
ggxy(data=ALL, obs=ALL$Qli.1.obs, sim=ALL$Qli.1.ERA5, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('ERA5 (W',m^-2,')',sep = '')),
     title = "Hourly Incoming Longwave Radiation", 
     xpos = 6.3/7, ypos = 1.3/7, outplot)

# HE3. MERRA2
outplot=file.path(outPath, "Qli_MERRA2.png")
ggxy(data=ALL, obs=ALL$Qli.1.obs, sim=ALL$Qli.1.MERRA2, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('MERRA2 (W',m^-2,')',sep = '')),
     title = "Hourly Incoming Longwave Radiation", 
     xpos = 6.3/7, ypos = 1.3/7, outplot)
##########

plotDailyPath <- './Output/DailyOutput/'
  
  
  
#######################             DAILY 
#######################                     
ALL.da <- aggDataframe(ALL, columns = c(1:27), period = 'daily', funs = 'mean')

#DA. Temp (y=1/3 for annotate) 
##########
# DA1. WFDEI
outplot=file.path(plotDailyPath, "Temp_WFDEI_daily.png")
ggxy(data=ALL.da, obs=ALL.da$t.1.obs.mean, sim=ALL.da$t.1.WFDEI.mean, 
     xlab=expression(paste('Observation (',degree, 'C)')), # [�C]
     ylab=expression(paste('WFDEI (',degree, 'C)')),
     title = "Daily Temperature", 
     xpos = 4/5, ypos = 1/5, outplot)

# DA2. ERA5
outplot=file.path(plotDailyPath, "Temp_ERA5_daily.png")
ggxy(data=ALL.da, obs=ALL.da$t.1.obs.mean, sim=ALL.da$t.1.ERA5.mean, 
     xlab=expression(paste('Observation (',degree, 'C)')), # [�C]
     ylab=expression(paste('ERA5 (',degree, 'C)')),
     title = "Daily TempERA5ture", 
     xpos = 4/5, ypos = 1/5, outplot)

# DA3. MERRA2
outplot=file.path(plotDailyPath, "Temp_MERRA2_daily.png")
ggxy(data=ALL.da, obs=ALL.da$t.1.obs.mean, sim=ALL.da$t.1.MERRA2.mean, 
     xlab=expression(paste('Observation (',degree, 'C)')), # [�C]
     ylab=expression(paste('MERRA2 (',degree, 'C)')),
     title = "Daily TempERA5ture", 
     xpos = 4/5, ypos = 1/5, outplot)
##########

#DB. Vap Pressure
#########
# DB1. WFDEI
outplot=file.path(plotDailyPath, "ea_WFDEI_daily.png")
ggxy(data=ALL.da, obs=ALL.da$ea.1.obs.mean, sim=ALL.da$ea.1.WFDEI.mean, 
     xlab= 'Observation (kPa)',
     ylab= 'WFDEI (kPa)',
     title = "Daily Vapour Pressure", 
     xpos = 6/7, ypos = 1/6, outplot)
# DB2. ERA5
outplot=file.path(plotDailyPath, "ea_ERA5_daily.png")
ggxy(data=ALL.da, obs=ALL.da$ea.1.obs.mean, sim=ALL.da$ea.1.ERA5.mean, 
     xlab= 'Observation (kPa)',
     ylab= 'ERA5 (kPa)',
     title = "Daily Vapour Pressure", 
     xpos = 6/7, ypos = 1/6, outplot)
# DB3. MERRA2
outplot=file.path(plotDailyPath, "ea_MERRA2_daily.png")
ggxy(data=ALL.da, obs=ALL.da$ea.1.obs.mean, sim=ALL.da$ea.1.MERRA2.mean, 
     xlab= 'Observation (kPa)',
     ylab= 'MERRA2 (kPa)',
     title = "Daily Vapour Pressure", 
     xpos = 6/7, ypos = 1/6, outplot)
##########

#DB. wind
##########
# C. Wind Speed (y=1/5, colour = 'darkblue' in annotate)
# C1. WFDEI
outplot=file.path(plotDailyPath, "u_WFDEI_daily.png")
ggxy(data=ALL.da, obs=ALL.da$u.1.obs.mean, sim=ALL.da$u.1.WFDEI.mean, 
     xlab= expression(paste('Observation (m',s^-1,')',sep = '')),
     ylab= expression(paste('WFDEI (m',s^-1,')',sep = '')),
     title = "Daily Wind Speed", 
     xpos = 6/7, ypos = 1/6, outplot)
# C2. ERA5
outplot=file.path(plotDailyPath, "u_ERA5_daily.png")
ggxy(data=ALL.da, obs=ALL.da$u.1.obs.mean, sim=ALL.da$u.1.ERA5.mean, 
     xlab= expression(paste('Observation (m',s^-1,')',sep = '')),
     ylab= expression(paste('ERA5 (m',s^-1,')',sep = '')),
     title = "Daily Wind Speed", 
     xpos = 6/7, ypos = 1/6, outplot)
# C3. MERRA2
outplot=file.path(plotDailyPath, "u_MERRA2_daily.png")
ggxy(data=ALL.da, obs=ALL.da$u.1.obs.mean, sim=ALL.da$u.1.MERRA2.mean, 
     xlab= expression(paste('Observation (m',s^-1,')',sep = '')),
     ylab= expression(paste('MERRA2 (m',s^-1,')',sep = '')),
     title = "Daily Wind Speed", 
     xpos = 6/7, ypos = 1/6, outplot)
##########

#DB. inc. sw
##########
# D. Inc. SW (y=1/5, colour = 'darkblue' in annotate)
# D1. WFDEI
outplot=file.path(plotDailyPath, "qsi_WFDEI_daily.png")
ggxy(data=ALL.da, obs=ALL.da$Qsi.1.obs.mean, sim=ALL.da$Qsi.1.WFDEI.mean, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('WFDEI (W',m^-2,')',sep = '')),
     title = "Daily Incoming Shortwave Radiation", 
     xpos = 6/7, ypos = 1/6, outplot)
# D2. ERA5
outplot=file.path(plotDailyPath, "qsi_ERA5_daily.png")
ggxy(data=ALL.da, obs=ALL.da$Qsi.1.obs.mean, sim=ALL.da$Qsi.1.ERA5.mean, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('ERA5 (W',m^-2,')',sep = '')),
     title = "Daily Incoming Shortwave Radiation", 
     xpos = 6/7, ypos = 1/6, outplot)

# D3. MERRA2
outplot=file.path(plotDailyPath, "qsi_MERRA2_daily.png")
ggxy(data=ALL.da, obs=ALL.da$Qsi.1.obs.mean, sim=ALL.da$Qsi.1.MERRA2.mean, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('MERRA2 (W',m^-2,')',sep = '')),
     title = "Daily Incoming Shortwave Radiation", 
     xpos = 6/7, ypos = 1/6, outplot)
##########

#DB. inc. LW
##########
# E. Inc. LW (y=1/5, colour = 'darkblue' in annotate)
# E1. WFDEI
outplot=file.path(plotDailyPath, "qli_WFDEI_daily.png")
ggxy(data=ALL.da, obs=ALL.da$Qli.1.obs.mean, sim=ALL.da$Qli.1.WFDEI.mean, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('WFDEI (W',m^-2,')',sep = '')),
     title = "Daily Incoming Longwave Radiation", 
     xpos = 6/7, ypos = 1/7, outplot)

# E2. ERA5
outplot=file.path(plotDailyPath, "qli_ERA5_daily.png")
ggxy(data=ALL.da, obs=ALL.da$Qli.1.obs.mean, sim=ALL.da$Qli.1.ERA5.mean, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('ERA5 (W',m^-2,')',sep = '')),
     title = "Daily Incoming Longwave Radiation", 
     xpos = 6/7, ypos = 1/7, outplot)

# E3. MERRA2
outplot=file.path(plotDailyPath, "qli_MERRA2_daily.png")
ggxy(data=ALL.da, obs=ALL.da$Qli.1.obs.mean, sim=ALL.da$Qli.1.MERRA2.mean, 
     xlab= expression(paste('Observation (W',m^-2,')',sep = '')),
     ylab= expression(paste('MERRA2 (W',m^-2,')',sep = '')),
     title = "Daily Incoming Longwave Radiation", 
     xpos = 6/7, ypos = 1/7, outplot)
##########

#plot precip daily
ALL.da <- aggDataframe(ALL, columns = c(1:27), period = 'daily', funs = 'sum')
ALL.da <- aggDataframe(ALL, columns = c(1:27), period = 'daily', funs = 'sum')


#DF. precip
#########
# DF1. WFDEI
outplot=file.path(plotDailyPath, "p_WFDEI_daily.png")
ggxy(data=ALL.da, obs=ALL.da$p.1.obs.sum, sim=ALL.da$p.1.WFDEI.sum, 
     xlab= 'Observation (mm)',
     ylab= 'WFDEI (mm)',
     title = "Daily Precipitation", 
     xpos = 6/7, ypos = 1/6, outplot)

# DB2. ERA5
outplot=file.path(plotDailyPath, "p_ERA5_daily.png")
ggxy(data=ALL.da, obs=ALL.da$p.1.obs.sum, sim=ALL.da$p.1.ERA5.sum, 
     xlab= 'Observation (mm)',
     ylab= 'ERA5 (mm)',
     title = "Daily Precipitation", 
     xpos = 6/7, ypos = 1/6, outplot)

# DB3. MERRA2
outplot=file.path(plotDailyPath, "p_MERRA2_daily.png")
ggxy(data=ALL.da, obs=ALL.da$p.1.obs.sum, sim=ALL.da$p.1.MERRA2.sum, 
     xlab= 'Observation (mm)',
     ylab= 'MERRA2 (mm)',
     title = "Daily Precipitation", 
     xpos = 6/7, ypos = 1/6, outplot)
##########


