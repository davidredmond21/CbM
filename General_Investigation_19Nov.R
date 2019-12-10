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
# one-liner functions

rms_d<-function(xx){
  yy<-diff(xx)
  sqrt(sum(yy^2)/length(yy))}
kurt_d<-function(xx){
  yy<-diff(xx)
  kurtosis(yy)}
rms_x = function(xx){( sqrt(sum(xx^2)/length(xx) ) )}
ent_f <-function(xx){ entropy(discretize(xx,numBins=length(xx))) }
hdf_db<-function(hdf_account){
  #h5ls(hdf_account)
  #axis1   = h5read(file=hdf_account, name="main/axis1")
  #blocki  = h5read(file=hdf_account, name="main/block0_items") ## 21 labels  
  axis0   = h5read(file=hdf_account, name="main/axis0")  #27 labels
  block0  = h5read(file=hdf_account, name="main/block0_values") ## 21 columns 
  #
  dataset<- t(block0)
  colnames(dataset)<-axis0
  return(dataset)
}
#
#
dataframe_dir = "C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data"
raw_csv_dir="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/raw_adc_data_csv/"
raw_hdf_dir="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/raw_adc_data/"

# Load data.
#### RAW plots
account= "merck"
#account= "dairygold"
hdf_files <- list.files(file.path(raw_hdf_dir),pattern="*hdf", full.names=TRUE)
hdf_samp  <- hdf_files[grep(hdf_files, pattern = paste0(account,"_........._adc_samps.hdf"))]

lsb_scale = 1/32.768001556396484
# loop to generate some TIME plots
for (i in 1:length(hdf_samp)){
  db_adc  <- hdf_db(hdf_samp[i])
  db_adc2  <- hdf_db(hdf_samp[i-1])
  
  sensor_id <- db_adc%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  sensor_id2 <- db_adc2%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  
  timestamp   <- db_adc%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  tim_axis <- seq(1,nrow(db_adc))/48000
  db_adc   <- apply(db_adc,2, function(x) {x*lsb_scale })
  db_adc2   <- apply(db_adc2,2, function(x) {x*lsb_scale })
  
  # take the last 4 smaples to use in the TIMESERIES plots
  n_col = ncol(db_adc)-5
  p_col = n_col - 3
  special_timestamps <- which(timestamp >"2018-10-11 00:00:00 UTC") [1:18]
  ref_timestamps <- which(timestamp >"2018-10-11 00:00:00 UTC")[1:4]
  fault_timestamps <- which(timestamp >"2018-11-11 00:00:00 UTC")[1:18]
  
  
  db_sub<- cbind(tim_axis[0:3000],db_adc[0:3000,fault_timestamps[c(1,18)]])%>%as.tibble()
  colnames(db_sub)<-c("msec","Oct-10","Nov-13")
  
  db_sub%>%as.tibble()%>%melt(id.vars=c("msec"))%>%
  ggplot(aes(msec,value,col=variable))+geom_line()+
  theme_bw()+
    ylab("Accel in g")+
    xlab("msec")+
  facet_wrap(~variable, scales="fixed")+
    theme(legend.position = "none")+
    ggtitle(paste(sensor_id,"raw data plot for samples", timestamp[fault_timestamps[1]],"and", timestamp[fault_timestamps[18]] ))

    ggsave(file=paste(sensor_id,"raw_data.pdf"),plot=last_plot(),dpi=900,path=raw_hdf_dir,width=25,height=22,units=c("cm"))
###
# loop to generate some Feature trending
#change the number of samples to the last 300
  n_col = ncol(db_adc)
  n_col2 = ncol(db_adc2)
  p_col  <-  max(c(n_col - 236),0)
  p_col2  <-  max(c(n_col2 - 236),0)
  
  db_sub  <- db_adc[,p_col:n_col]
  db_sub2 <- db_adc2[,p_col2:n_col2]
  timestamp   <- db_sub%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  
  df<- cbind(timestamp,rms_dif,std_vec,krt_vec,rms_vec)%>%as.tibble()
  df$timestamp%<>%as.POSIXct(tz='UTC',origin="1970-01-01")
  df2<- cbind(timestamp,std_vec,std_vec2)
  colnames(df2) <-c("timestamp",sensor_id,sensor_id2)
  df2%<>%as.tibble
  df2$timestamp%<>%as.POSIXct(tz='UTC',origin="1970-01-01")
  rms_dif <-apply(db_sub,2,rms_d)%>%as.vector()
  rms_vec <-apply(db_sub,2,rms)%>%as.vector()
  std_vec <-apply(db_sub,2,sd)%>%as.vector()
  std_vec2 <-apply(db_sub2,2,sd)%>%as.vector()
  krt_vec <-apply(db_sub,2,kurtosis)%>%as.vector()
  df2%>%filter(timestamp>"2018-11-05")%>%melt(id.vars=c("timestamp"))%>%
    ggplot(aes(timestamp,value,col=variable))+geom_point()+
    theme_stata()+
    ylab("Std Dev Values")+
    xlab("Sample Dates")+
    #facet_wrap(~variable)+
    #theme(legend.position = "none")+
    scale_x_datetime()  +
    ggtitle(paste(sensor_id,sensor_id2,"trending features"))
  
  
    ggsave(file=paste(sensor_id,"tending_features.pdf"),plot=last_plot(),dpi=900,path=raw_hdf_dir,width=25,height=22,units=c("cm"))
    }
# 
hdf_ffts  <- hdf_files[grep(hdf_files, pattern = paste0(account,"_........._adc_ffts.hdf"))]

for (i in 1:length(hdf_ffts)){
  db_fft      <- hdf_db(hdf_ffts[i])
  db_fft2      <- hdf_db(hdf_ffts[i-1])
  
  sensor_id   <- db_fft%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  sensor_id2   <- db_fft2%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  timestamp   <- db_fft%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  
  timestamp2   <- db_fft2%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  n_col = ncol(db_fft)
  n_col2 = ncol(db_fft2)
  normal=  737-9 
  abnorm= 866-9
#  abnorm = 2413
#  abnorm2 = 903
 
  db_fft_plt<- db_fft[,c(1,ref_timestamps[2],  fault_timestamps[18])]
  db_fft_plt%<>% as.tibble()
  
  colnames(db_fft_plt)<-c("frequencies","Oct-10","Nov-13")
  harmonics <- seq(1:12) * db_fft_plt$`Oct-10`[0:80]%>%which.max()
  
 # no_fault  <-  db_fft_plt$Normal[0:580]%>%which.max()
  #fault     <-  db_fft_plt$Abnormal[0:580]%>%which.max()
  
  db_fft_plt%>%as.tibble()%>%melt(id.vars=c("frequencies"))%>%filter(frequencies<2600)%>%
    ggplot(aes(frequencies,value,col=variable))+geom_line()+
    #theme_stata()+
    theme_bw()+
    ylab("Amp Accel ")+
    xlab("Frequency ")+
    #theme(legend.position = "none")+
    facet_wrap(~variable, scales="fixed")+
    geom_vline(xintercept=harmonics,col="green", lty="dotted")+
    #geom_vline(xintercept=fault,col="blue", lty="dashed")+
    #geom_vline(xintercept=no_fault,col="purple", lty="dashed")+
    ggtitle(paste(sensor_id," Pump B P2102 for samples ",timestamp[ref_timestamps[1]], "and", timestamp[fault_timestamps[16]] ))

  
  ############ Low Frequency
  db_fft_plt%>%as.tibble()%>%melt(id.vars=c("frequencies"))%>%filter(frequencies<2600)%>%
    ggplot(aes(frequencies,value,col=variable))+geom_line()+
    #theme_stata()+
    theme_bw()+
    ylab("Amp Accel ")+
    xlab("Frequency ")+
    #theme(legend.position = "none")+
    facet_wrap(~variable, scales="fixed")+
    geom_vline(xintercept=harmonics[1]*7*4.02,col="green", lty="dotted")+
    #geom_vline(xintercept=fault,col="blue", lty="dashed")+
    #geom_vline(xintercept=no_fault,col="purple", lty="dashed")+
    ggtitle(paste(sensor_id," Pump B P2102 for samples ",timestamp[ref_timestamps[1]], "and", timestamp[fault_timestamps[18]] ))
  
  
  
  
  db_fft_plt2%>%as.tibble()%>%melt(id.vars=c("frequencies"))%>%filter(frequencies<3000)%>%
    ggplot(aes(frequencies,value,col=variable))+geom_line()+
    theme_stata()+
    ylab("Amp Accel ")+
    xlab("Frequency ")+
    #theme(legend.position = "none")+
    facet_wrap(~variable)+
    ggtitle(paste(sensor_id2," Pump A P2102  ",timestamp2[normal2], "and ", timestamp2[abnorm2] ))
  
  
  ggsave(file=paste(sensor_id,"FFT_features.pdf"),plot=last_plot(),dpi=900,path=raw_hdf_dir,width=25,height=22,units=c("cm"))

  ## Last Weeks
  db_fft_lst<- db_fft[,c(1,n_col)]
  db_fft_lst%>%as.tibble%>%melt(id.vars=c("frequencies"))%>%filter(frequencies<2000)%>%
    ggplot(aes(frequencies,value,col=variable))+geom_line()+
    theme_stata()+
    ylab("Amp Accel ")+
    xlab("Frequency ")+
    theme(legend.position = "none")+
    ggtitle(paste(sensor_id,"Most Recent Spectral features ",timestamp[n_col]))
  ggsave(file=paste(sensor_id,"FFT_latest.pdf"),plot=last_plot(),dpi=900,path=raw_hdf_dir,width=25,height=22,units=c("cm"))
  
  ## High Freq Spectrum
  db_fft_plt%>%as.tibble()%>%melt(id.vars=c("frequencies"))%>%filter(frequencies<3000)%>%filter(frequencies>800)%>%
    ggplot(aes(frequencies,value,col=variable))+geom_line()+
    theme_stata()+
    ylab("Amp Accel ")+
    xlab("Frequency ")+
    theme(legend.position = "none")+
    facet_wrap(~variable)+
    ggtitle(paste(sensor_id,"High Freq Spectral features ",timestamp[n_col]))
  ggsave(file=paste(sensor_id,"FFT_Hi_freq_features.pdf"),plot=last_plot(),dpi=900,path=raw_hdf_dir,width=25,height=22,units=c("cm"))
  }




