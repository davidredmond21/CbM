library(moments)
library(ggplot2)
library(robustbase)
library(reshape2)
library(tibble)
library(dplyr)
library(entropy)
library(readr)
library(magrittr)
library(spectral)
library(seewave)
library(viridis)
library(lubridate)
library(stringr)
library(e1071)
library(readr)
library(ggthemes)
# Source all the libraries and directory path names
# Summary files first to generate TRENDS
library(rhdf5)
dataframe_dir = "C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/dataframes"
# Load data.

hdf_account="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/raw_adc_data/merck_MCSVX5074_adc_ffts.hdf"

hdf_account="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/dataframes/merck_summary.hdf"

hdf_db<-function(hdf_account){
  h5ls(hdf_account)
  df   = h5read(file=hdf_account, name="main")  #27 labels
  axis0   = h5read(file=hdf_account, name="main/axis0")  #27 labels
  axis1   = h5read(file=hdf_account, name="main/axis1")
  blocki  = h5read(file=hdf_account, name="main/block0_items") ## 21 labels  
  block0  = h5read(file=hdf_account, name="main/block0_values") ## 21 columns 
  #
  dataset<- t(block0)
  colnames(dataset)<-axis0
  return(dataset)
}
hdf_account="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/raw_adc_data/merck_MCSVE6B8S_adc_samps.hdf"
hdf_account="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/raw_adc_data/merck_MCSVX5074_adc_samps.hdf"

db_adc<- hdf_db(hdf_account)

# pick a sensor 
sensor_list <- db_adc%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
timestamp   <- db_adc%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
# split the sensor_list into 4 
tim_axis <- seq(1,nrow(db_adc))/48000
#db_sub<- cbind(tim_axis[0:8000],db_adc[0:8000,c(2705:2707)])
db_sub<- cbind(tim_axis[0:8000],db_adc[0:8000,c(700:704)])%>%as.tibble()
glimpse(db_sub)
str(db_sub)
db_sub%>%as.tibble()%>%melt(id.vars=c("V1"))%>%
  ggplot(aes(V1,value,col=variable))+geom_line()+
  #facet_wrap(~variable,scale="free")+
  ggtitle(paste(sensor_list,"raw data plot before, during and after significant event"))


### reading in the SUMMARY data : Initially to look at the RPM 
#db_sum <- read_csv(summ_account)
db_sum<-hdf_db(hdf_account)
db_sum$timestamp%<>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
db_sum$sensorid%<>%as.factor()
db_sum$rmsvel  %<>%as.numeric()
db_sum$filt_rms  %<>%as.numeric()
db_sum$rpm  %<>%as.numeric()
sensor_list <- db_sum$sensorid%>%as.factor%>%levels
list_a <- sensor_list
rpm_est<-db_sum%>%filter(machine_on==1)%>%select(rpm,sensorid)%>%aggregate(.~sensorid,data=.,FUN="getmode")





db_fft      <- read_csv(fft_account[i])
timestamp   <- db_fft%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
sensor_id   <- db_fft%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
rpm_id      <- rpm_est[which(rpm_est$sensorid==sensor_id),]%$%rpm
harmonics    <- seq(1:20)*rpm_id/60  

db_fft_plt<- cbind(db_fft[,1],db_fft[,c(2706:2708)])
colnames(db_fft_plt) <- c("freq", "before","during","after")
db_fft_plt%>%as.tibble()%>%melt(id.vars=c("freq"))%>%filter(freq<1200)%>%
  ggplot(aes(freq,value,col=variable))+geom_line()+
  ggtitle(paste(sensor_id,"spectrum  plot before, during and after significant event"))+
  geom_vline(xintercept = c(harmonics),col="red",linetype="dashed")



db_fft_plt%>%as.tibble()%>%melt(id.vars=c("freq"))%>%filter(freq<400)%>%
  ggplot(aes(freq,value,col=variable))+geom_line()+
  ggtitle(paste(sensor_id,"spectrum  plot before, during and after significant event"))+
  geom_vline(xintercept = c(harmonics[1:8]),col="red",linetype="dashed")





db_sum%>%
  filter(machine_on==1)%>%
  select(timestamp,stdDev,skewness,sensorid, rms,ptp,min,max,kurtosis,crest)%>%
  filter(sensorid==sensor_id)%>%na.omit()%>%
  filter(timestamp>"2018-10-30")%>%
  melt(id.vars=c("timestamp"))%>%
  ggplot(aes(timestamp,value,col=variable))+geom_point()+
  facet_wrap(~variable,scale="free")





#points(db_fft$frequencies[1:1000],db_fft[1:1000,(n_fft-1)],type="l",lty=3,col="green")
#points(db_fft$frequencies[1:1000],db_fft[1:1000,(n_fft-2)],type="l",lty=3,col="green")
rpm_est <- which.max(db_fft[1:100,n_fft])-1
rpm_lines = seq(1,12,1)*(rpm_est/1)
abline(v=rpm_lines,col="forestgreen",lty=3)
for( j in 1:12){
  text((trunc(rpm_lines[j])-10),7,paste0(j,"X rpm"),col="red",cex=0.75)
}
# takes out the harmonic energy at the RPM and +2, -1
db_fft[rpm_lines[1],n_fft]  <-0
db_fft[rpm_lines[1]+1,n_fft]<-0
db_fft[rpm_lines[1]-1,n_fft]<-0
db_fft[rpm_lines[1]+2,n_fft]<-0
spur_line <- which.max(db_fft[1:1000,n_fft])-1 
abline(v=spur_line,col="purple",lty=3,lwd=3)
text((spur_line+5),5,paste0(spur_line,"Hz next highest harmonic"),col="purple",cex=0.8)
## adding reference plots


############### ~ Monday 
raw_dir="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/raw_adc_data_csv/"


rms_d<-function(xx){
  yy<-diff(xx)
  sqrt(sum(yy^2)/length(yy))}
kurt_d<-function(xx){
  yy<-diff(xx)
  kurtosis(yy)}
#function for RMS difference
rms_x = function(xx){( sqrt(sum(xx^2)/length(xx) ) )}
ent_f <-function(xx){ entropy(discretize(xx,numBins=length(xx))) }


for (i in 1:6) {
  db_samps <- read_csv(adc_account[i])
  #glimpse(db_samps)
  timestamp   <- db_samps%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  sensor_id <- db_samps%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  time_axis <- seq(1,nrow(db_samps))/48000
  db_ssamps   <-  cbind(time_axis[0:12000],db_samps[0:12000,])
  #db_1<- db_samps[,2707]
  s_ncol <- ncol(db_samps)
  p_col  <- s_ncol-300
  db_sub      <- db_samps[,p_col:s_ncol]
  timestamp   <- db_sub%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  rms_dif <-apply(db_sub,2,rms_d)%>%as.vector()
  std_vec <-apply(db_sub,2,sd)%>%as.vector()
  krt_vec <-apply(db_sub,2,kurtosis)%>%as.vector()
  df<- cbind(timestamp,rms_dif,std_vec,krt_vec)%>%as.tibble()
  df$timestamp%<>%as.POSIXct(tz='UTC',origin="1970-01-01")
  df%>% melt(id.vars=c("timestamp"))%>%
    ggplot(aes(timestamp,value,col=variable))+geom_line()+
    scale_x_datetime()+
    ggtitle(paste(sensor_id,"trending features"))
    ggsave(file=paste(sensor_id,"tending_features.pdf"),plot=last_plot(),dpi=900,path=raw_dir)#,width=1080,height=182,units=c("mm"))
}
#  facet_wrap(~variable)
# When did it happen?
timestamp[which.max(std_vec) ]
timestamp[which.max(std_vec[0:481])]
# Prior to that when was the last peak 

library(tidyverse)
library(forecast)
install.packages(c("rgenoud","parallel","furrr","tsibble","brotools"))
library(rgenoud)
library(parallel)
library(lubridate)
library(furrr)
library(tsibble)
library(brotools)
ihs <- function(x){
  log(x + sqrt(x**2 + 1))
}

to_tibble <- function(forecast_object){
  point_estimate <- forecast_object$mean %>%
    as_tsibble() %>%
    rename(point_estimate = value,
           date = index)
  
  upper <- forecast_object$upper %>%
    as_tsibble() %>%
    spread(key, value) %>%
    rename(date = index,
           upper80 = `80%`,
           upper95 = `95%`)
  
  lower <- forecast_object$lower %>%
    as_tsibble() %>%
    spread(key, value) %>%
    rename(date = index,
           lower80 = `80%`,
           lower95 = `95%`)
  reduce(list(point_estimate, upper, lower), full_join)
}
order_list <- list("p" = seq(0, 3),
                   "d" = seq(0, 2),
                   "q" = seq(0, 3)) %>% cross() %>% map(lift(c))
list("p" = seq(0, 3),
     "d" = seq(0, 2),
     "q" = seq(0, 3)) %>%
  cross() %>%
  head(3)


