#urpose of this script is
#.libPaths(Sys.getenv('RLIBPATHS'))
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
#
options(mc.cores = parallel::detectCores())
#
# Functions
#
abs_max = function(xx){ max(abs(xx))}   ## F1  ( need to auto-scale it)
ave_abs = function(xx){ sum(abs(xx))/length(xx)}  ##F2 (need to auto-scale)
pk_pk = function(xx){ max((xx)-min(xx)) }  ## F3 (auto-scaled)
## F4 variance....
## F5 is Std dev
## F6 is skewness
## f7 is kurtosis
## f8 is rms
rms_x = function(xx){( sqrt(sum(xx^2)/length(xx) ) )}
## F9
crest = function(xx){ ( max(abs(xx)/(sqrt(sum(xx^2))/length(xx)))) }
## F10
clearence = function(xx) { max(abs(xx))/rms_x(xx) }
## F11
impulse_f = function(xx) { max(abs(xx))/ave_abs(xx) }
## F12
shape_f = function(xx)  { rms_x(xx)/ave_abs(xx) }
# Filtering to 12.5kHz on a vector
ent_f <-function(xx){entropy(discretize(xx,numBins=length(xx))) }
mah   =  function(xx,mah_offset=0) {mahalanobis(xx,colMeans(xx),cov(xx))-mah_offset }
###### - reserving these arrays variables
abs_vec <- 0
sd_vec  <- abs_vec
skw_vec <- abs_vec
krt_vec <- abs_vec
rms_vec <- abs_vec
#
lmod<-function(x,bw=5000){scale(x,center=T,scale=FALSE)%>%fft%>%Mod}
# Set up basics Bandwidth of interest upto 5kHz
bw=5000
# ignore DC
dc_offset=2
# Given RPM is 
rpm=2870/60%>%trunc()  
#
#
# Simple display of the Frequency content for a specific sample ####
# Prior Knowledge is the RPM speed for the Vlines
#
##### Select which ACCOUNT 
account= "boortmalt"
account= "dairygold"
account= "medite"
account= "merck"
#
#########################
root_dir<-paste("/home/david/Downloads/MediteData")
#root_dir<-paste("C:/Users/dredmond/Desktop/DataSet1/",account, sep="")
tot_dir =list.dirs(root_dir,full.names=TRUE,  recursive=T)
use_dir_raw<-tot_dir[grep(tot_dir, pattern = "raw")]
use_rpm_raw<-tot_dir[grep(tot_dir, pattern = "raw/time_freq_plots")]
########
infiles= 0
j=1
i=1
h=1
db.mod <- as.matrix(NA,nrow=5000,ncol=ncol)
#
#
#### RAW plots
for (h in 1:length(use_dir_raw))
   {
  infiles =list.files(file.path(use_dir_raw[h]),pattern="*.csv", full.names=TRUE)
  #rpm_file<-list.files(file.path(use_rpm_raw[h]),pattern="*.csv", full.names=TRUE)
  #rpm_est<- read_delim(rpm_file,delim=",",col_types = "?d", col_names = F,skip=1)
  #rpm_est$X1 <- as.POSIXct(rpm_est$X1,tz='UTC',origin="1970-01-01")
  db1 <- read_csv(infiles[112],col_names = FALSE)
  db2 <- read_csv(infiles[113],col_names = FALSE)
  sensor_id<-str_extract(infiles[113],"MCSV.....")}
#  
  ncol=length(infiles)
  #
  
  db.mod <- matrix(NA,nrow=1000,ncol=length(infiles))
    for (i in 1:length(infiles))  {
      db<-read_csv(infiles[i],col_names = "X" )
      sd_vec[i]  <- apply(db,2,sd)
      skw_vec[i]  <- apply(db,2,skewness)
      krt_vec[i]  <- apply(db,2,kurtosis)
      rms_vec[i]  <- apply(db,2,rms_x)
      db.mod[,i] <- apply(db,2,lmod)[1:1000]%>%log()
      i=i+1
    }

  ####### nwaterfall
  idx=seq_along(db.mod[,1])
  jdx=seq_along(db.mod[1,])

 
  library(rgl)
  open3d()
  bg3d("white")
  surface3d(idx,jdx,db.mod)
  plot3d(idx,jdx,db.mod)
  axes3d()
  for (i in 1:nx) lines3d(x[i], y, z[i, ], col = 'white', lwd = 2)
  rgl.postscript("persptrial_060514.eps","eps")
  
  
  
  
  
### ON/OFF machine
  mc_state <- kmeans(sd_vec,2)$cluster
  ggplot(rpm_est,aes(X1,X2,col=mc_state))+geom_line()
#
#sensor_list<-str_extract(infiles,"MCSV....")%>%unique()%>%as.array
#
  db_fft1<-db1%>%scale(.,center=T,scale=FALSE)%>%fft%>%Mod%>%log()
  db_fft2<-db2%>%scale(.,center=T,scale=FALSE)%>%fft%>%Mod%>%log()
  freq= seq(from =1, to = length(db_fft1))
  db_fft<-cbind(db_fft1,db_fft2,freq)%>%as.tibble
  colnames(db_fft)<-c("X1","X2","Freq")
#
# a generic FFT plot of 2 samples
#
  rpm=2870/60%>%trunc()  
  RPM_lines= seq(from =1, to =80)*rpm
  dc_offset=2.00
  lim_freq = 3000
  db_fft%>%slice(.,dc_offset:lim_freq) %>%
  ggplot(aes(Freq,X1))+
  geom_line() +
  geom_line(aes(Freq,X2),col="green")+
  xlab("Frequency Hz")+
  ylab("Amplitude dB")+
  scale_x_log10()+
  scale_y_log10()+
  geom_vline(xintercept = RPM_lines,col="blue", alpha=0.24)+
  scale_fill_discrete(name="Frequency\nSpectrum",
                      breaks=c("Sample #1", "Successive sample"),
                      labels=c("Frequency", "4 hours later"))+
  annotate("text", x= (rpm+25), y= 10e-2, label= "RPM 24.3Hz and harmonics",col="blue")+
  ggtitle(paste(sensor_id,"High frequency Harmonics for" ,subtitle = "2 successive samples"))
  #
  ggsave(file=paste(sensor_id,"freq_plot.jpg"),plot=last_plot() ,dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
#
#db1 <- read_csv(infiles[121],col_names = FALSE)
#db2 <- read_csv(infiles[122],col_names = FALSE)
#
  adc_sample<- seq(from =1, to = nrow(db1))
  db_raw<- cbind(db1,db2,adc_sample)
  colnames(db_raw)<-c("X1","X2","Sample")
  start_point = 1
  stop_point  = trunc(48000/rpm)*3  # shows 3 cycles
  RPM_period  = seq(from =1, to =3)*stop_point/6
  db_raw%>%slice(.,start_point:stop_point) %>%
  ggplot(aes(y=X1,x=Sample,alpha=0.5))+
  geom_point()+
  geom_jitter(width = 0.0,height = 0.015)+
  xlim(start_point,stop_point)+
  xlab("Vibration measurements")+
  ylab("Acceleration in g")+
  geom_hline(yintercept = c(max(db1$X1),min(db1$X1)),alpha=0.4, col="darkgreen")+
  geom_vline(xintercept = RPM_period,col="blue", alpha=0.24)+
  annotate ("text",x=1100,y=min(db1$X1)+0.04,label="Lower bound", color="darkgreen")+
  annotate ("text",x=1100, y= max(db1$X1)+0.04,label="Upper bound", color="darkgreen")+
  annotate("text", x= (stop_point/4+900), y= min(db1$X1)+0.04, label= "RPM period",col="blue")+
  theme(legend.position = "none")+
  ggtitle(paste(sensor_id,"RAW data measurements",subtitle = " sub-sample repeated every 4 hours"))
  #
  ggsave(file=paste(sensor_id,"Vibration_Data.jpg"),plot=last_plot(),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
#
########################### Lets look at FAULT frequencies
#
lmod<-function(x,bw=5000){scale(x,center=T,scale=FALSE)%>%fft%>%Mod}
# Set up basics Bandwidth of interest upto 5kHz
bw=5000
# ignore DC
dc_offset=2
# Given RPM is 
rpm=2870/60%>%trunc()  
#
### building a sample 
db<-  read_delim(infiles[2],delim=",",col_types = "d")
db.lmod <- apply(db,2,lmod)%>%log%>%as.tibble()
#colnames(db.lmod) <-c("X_axis","Y_axis","Z_axis")
colnames(db.lmod)  <-c("X-axis")
ref_freq<- seq(from=1,to=nrow(db.lmod))
RPM_lines= seq(from =2, to =bw/rpm, by=12)*rpm%>%trunc
#
### Algorithm
#
# Start with a GIVEN rpm- value
# Build the fault frequencies table
# Fault freq, and their integer harmonics for XYZ
# With 1 Hz resolution fomr Zybo - calculate the Fault Energy  
# Calculate the Energy RATIO from 2 consequative measurement samples
#
# From the name plate 6202 ----> 
fault_mult=c(3.58,5.42,0.40,2.34)
fault_name=c("BPFO","BPFI","FTF","BSF")
fault_freq=rpm*fault_mult
harm_ff = as.matrix(t(rpm*fault_mult))
colnames(harm_ff)<-fault_name
reps=40
XX<-matrix(rep(harm_ff,reps),nrow=ncol(harm_ff),ncol=reps)
full_ff<- sweep(t(XX),1,1:reps,`*`)%>%trunc
colnames(full_ff)<-fault_name
full_ff1<- sweep(t(XX),1,1:reps,`*`)%>%trunc%>%+1
colnames(full_ff1)<-fault_name
full_ff2<- sweep(t(XX),1,1:reps,`*`)%>%trunc%>%-1
colnames(full_ff2)<-fault_name
full_ff3<- sweep(t(XX),1,1:reps,`*`)%>%trunc%>%+2
colnames(full_ff3)<-fault_name
full_ff<-rbind(full_ff,full_ff1,full_ff2)%>%as.tibble() 
#full_ff<-rbind(full_ff,full_ff1,full_ff2,full_ff3)%>%as.tibble() 
#  full_ff is the list of freq where we may have Fault harmonics +/- 1Hz 
# this acts as the index to db.lmod.
# For each axis sum the fault freq, for each fault type
# initialize the dataframes
db_bpfo=db.lmod[1:length(infiles),]%>%cbind(.,"idx"=1)
db_bpfi=db.lmod[1:length(infiles),]%>%cbind(.,"idx"=1)
db_ftf =db.lmod[1:length(infiles),]%>%cbind(.,"idx"=1)
db_bsf =db.lmod[1:length(infiles),]%>%cbind(.,"idx"=1)
j=1
k <- 3:7 # these are the important fault freq harmonics 
for (i in infiles){
  db = read_delim(i,delim=",",col_types = "d")
  db.lmod <- apply(db,2,lmod)%>%as.tibble() # generate the FFT,Log,MOD
  colnames(db.lmod) <-c("X-axis")
  db_bpfo[j,]<-db.lmod[full_ff$BPFO[k],]%>%apply(.,2,rms)%>%log%>%t()%>%cbind(.,j)%>%as.tibble()
  db_bpfi[j,]<-db.lmod[full_ff$BPFI[k],]%>%apply(.,2,rms)%>%log%>%t()%>%cbind(.,j)%>%as.tibble()
  db_ftf[j,] <-db.lmod[full_ff$FTF[k],] %>%apply(.,2,rms)%>%log%>%t()%>%cbind(.,j)%>%as.tibble()
  db_bsf[j,] <-db.lmod[full_ff$BSF[k],] %>%apply(.,2,rms)%>%log%>%t()%>%cbind(.,j)%>%as.tibble()
  j=j+1
}
db_bpfo%<>%cbind(.,mc_state)
db_bpfi%<>%cbind(.,mc_state)
db_ftf%<>%cbind(.,mc_state)
db_bsf%<>%cbind(.,mc_state)
#
#
# Grid plots
# install.packages("ggpubr")
library(ggpubr)
p1=filter(db_bpfo,mc_state==2)%>%select("X-axis","idx")%>%melt(.,("idx"))%>%
  ggplot(aes(idx,value,colour=variable))+geom_point()+
  theme_bw()+
  ylab("Log scale")+
  xlab("Measurement sample ")
p2=filter(db_bpfi,mc_state==2)%>%select("X-axis","idx")%>%melt(.,("idx"))%>%
  ggplot(aes(idx,value,colour=variable))+geom_point()+
  theme_bw()+
   ylab("Log scale")+
  xlab("Measurement sample ")
 p3=filter(db_ftf,mc_state==2)%>%select("X-axis","idx")%>%melt(.,("idx"))%>%
  ggplot(aes(idx,value,colour=variable))+geom_point()+
  theme_bw()+
  ylab("Log scale")+
  xlab("Measurement sample ")
p4=filter(db_bsf,mc_state==2)%>%select("X-axis","idx")%>%melt(.,("idx"))%>%
  ggplot(aes(idx,value,colour=variable))+geom_point()+
  theme_bw()+  
   ylab("Log scale")+
  xlab("Measurement sample ")
#
plt<-ggarrange(p1,p2,p3,p4,ncol=2,nrow=2, hjust = -2.2,labels = c("BPFI","BSFO","FTF","BSF"),legend="top",common.legend = TRUE)
annotate_figure(plt,
                top=text_grob(paste("Comparison of RMS=",rpm,"Trends in bearing Fault Freq energy for X-axis")),
                bottom=text_grob(paste("Dataset",sensor_id, "May-June 2018"),
                                 hjust=1,x=1,face="italic",size=10))
ggsave(paste("Harmonic Energy ALL faults",sensor_id,".jpg",sep=""),plot=last_plot(),dpi=900,width=25,height=22,units=c("cm"),path=root_dir)
h=h+1
}
#
#
### look for ON/ OFF state using std_dev ###############
#
#
library(jsonlite)
library(tibble)
use_dir_summ<-tot_dir[grep(tot_dir, pattern = "summary")]
#### RAW plots
# create a dummy version of infiles for summary
infiles=list.files(use_dir_summ,pattern="*.json",full.names = T,recursive = T)
#
no_col <- length(infiles)
db_full  <- fromJSON(infiles[1])%>%as_data_frame()%>%select(-(c("build","channel","lsb_per_g","sample_size")))
for(l in 2:no_col){
  db_full[l,]  <- fromJSON(infiles[l])%>%as_data_frame()%>%select(-(c("build","channel","lsb_per_g","sample_size")))
}
db_full$timestamp%<>%as.POSIXct(.,tz='UTC',origin="1970-01-01")
db_full$machine_state%<>%as.factor(.)
write_rds(db_full,"db_summary.rds")
#######  Mahalano Plots , Box plots, and Density plots
library(viridis)
library(ggthemes)
library(dplyr)
sensor_list <- db_full$sensorid%>%as.factor%>%levels
for (i in seq_along(sensor_list)) {
  y <-db_full%>%filter(sensorid==sensor_list[i])%>%filter(machine_state==1)
  y.mah<- y%>%select(c(kurtosis,rms))%>%mah
  y<-cbind(y,y.mah)
  y.long <- y%>%select(stdDev,skewness,shape,rms,sensorid,kurtosis,impulse,crest,clearance,y.mah)%>%melt(id.vars=c("sensorid"))

# BOx plots to show outliers  
    ggplot(y.long,aes(x=variable,y=value,col=variable)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 1)+
    facet_wrap(~variable,scale="free") +
    xlab("Condition Indicators  ")+
    ylab("Value")+
    ggtitle(paste(sensor_list[i],"Box Plots of Condition Indicators machine_state = ON ") )
    ggsave(file=paste(sensor_list[i],"_box.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
# Histogram or Density plots to show Distribution
  ggplot(y.long,aes(value,col=variable))+
    geom_histogram(aes(value),bins = 90) +
    facet_wrap(~variable,scale="free")+
    xlab("Selection of Condition Indicators ")+
    ylab("Density") +
    ggtitle((paste(sensor_list[i],"Density plots of select Condition Indicators machine_state= ON")))
    ggsave(file=paste(sensor_list[i],"_density.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))

# TimeSeries plots to show trends 
y%>%select(stdDev,timestamp,skewness,shape,rms,sensorid,kurtosis,impulse,crest,clearance,y.mah)%>%
    melt(id.vars=c("sensorid","timestamp"))%>%
    ggplot(aes(timestamp,value,col="Conditions"))+
      geom_line(aes(timestamp,value,color=variable))+
      facet_wrap(~variable,scales="free")+
      xlab("Trend of Condition Indicators ")+
      ggtitle((paste(sensor_list[i],"TimeSeries plots of select Condition Indicators machine_state= ON")))
    ggsave(file=paste(sensor_list[i],"_timeseries.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))

#Trend for 2 indicators Mahalanobsis and StdDev
  y%>%filter(sensorid==sensor_list[i] )%>%filter(timestamp>"2018-04-21")%>%
   select(timestamp,rms,kurtosis,crest,rms,stdDev,y.mah,skewness)%>%
  melt(id.vars=c("timestamp"))  %>%
  ggplot(aes(timestamp,value,col=variable))+  geom_point() +
  geom_smooth(method = 'gam')+#, span=0.5, level=0.99)+
  facet_wrap(~variable,scale="free") +
  xlab("Trend of Mahalanobsis Indicators ")  +
  ggtitle((paste(sensor_list[i],"TimeSeries plots of select Condition Indicators")))
  ggsave(file=paste(sensor_list[i],"_timeseries.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
# individual STD-DEV models
  db_full%>%filter(machine_state==1)%>%filter(sensorid==sensor_list[i])%>%filter(timestamp>"2018-05-01")%>%
    select(timestamp,stdDev,sensorid)%>%
    ggplot(aes(timestamp,stdDev,col=sensorid))+  geom_point() +
    geom_smooth(method = 'auto')+#, span=0.5, level=0.99)+
    xlab("StdDev Indicator ")  +
    ggtitle(paste(sensor_list[i],"plots of StdDev Indicator"))
  ggsave(file=paste(sensor_list[i],"_StdDev_timeseries.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
  
}
# All sensors together
db_full%>%filter(machine_state==1)%>%filter(timestamp>"2018-05-01")%>%
  select(timestamp,rms,kurtosis,stdDev,sensorid)%>%
  melt(id.vars=c("timestamp","sensorid"))  %>%
  ggplot(aes(timestamp,value,col=variable))+  geom_point() +
  geom_smooth(method = 'gam')+#, span=0.5, level=0.99)+
  facet_wrap(~sensorid,scale="free") +
  xlab("Trend Indicators ")  +
  ggtitle("TimeSeries plots of select Condition Indicators")
ggsave(file=paste(sensor_list[i],"_timeseries.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))




db_full%>%filter(machine_state==1)%>%filter(timestamp>"2018-05-01")%>%
  select(timestamp,stdDev,sensorid)%>%
  melt(id.vars=c("timestamp","stdDev"))  %>%
  ggplot(aes(timestamp,stdDev,col=value))+  geom_point() +
  geom_smooth(method = 'auto')+#, span=0.5, level=0.99)+
  facet_wrap(~value,scale="free") +
  xlab("StdDev Indicator ")  +
  ggtitle("TimeSeries plots of StdDev Indicators for all sensor_id")






