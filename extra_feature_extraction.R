#Purpose of this script is
## Libraries and helper functions
#.libPaths(Sys.getenv('RLIBPATHS'))
#install.packages(c("moments","ggplot2","robustbase","reshape2","tibble","dplyr","readr","magrittr","spectral","seewave","viridis","lubridate","stringr","e1071","readr","ggthemes","jsonlite","BayesLCA"))

#library(moments)
library(ggplot2)
library(robustbase)
library(reshape2)
library(tibble)
library(dplyr)
library(readr)
library(magrittr)
library(spectral)
#library(seewave)
library(viridis)
library(lubridate)
library(stringr)
library(e1071)
library(readr)
#library(ggthemes)
#library(jsonlite)
#library(BayesLCA)
#
options(mc.cores = parallel::detectCores())
#  One liner helper Functions
abs_max <- function(xx){ max(abs(xx))}   ## F1  ( need to auto-scale it)
ave_abs <- function(xx){ sum(abs(xx))/length(xx)}  ##F2 (need to auto-scale)
pk_pk <- function(xx){ max((xx)-min(xx)) }  ## F3 (auto-scaled)
rms_x <- function(xx){( sqrt(sum(xx^2)/length(xx) ) )}
crest <- function(xx){ ( max(abs(xx)/(sqrt(sum(xx^2))/length(xx)))) }
clearence <- function(xx) { max(abs(xx))/rms_x(xx) }
impulse_f <- function(xx) { max(abs(xx))/ave_abs(xx) }
shape_f <- function(xx)  { rms_x(xx)/ave_abs(xx) }
ent_f <-function(xx){entropy(discretize(xx,numBins=length(xx))) }
mah   <-function(xx,mah_offset=0) {mahalanobis(xx,colMeans(xx),cov(xx))-mah_offset }
lmod  <-function(x,bw=24000){
  y <- scale(x,center=T,scale=FALSE)%>%fft%>%Mod
  return (y[1:bw])}
getmode <- function(v) {uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
acc2vel <-function(x,bw=1000) {
  # fill in the code to convert accel to Velocity .
  # Pick the mid point of the BW and set a 6db LP filter
  # then calcul the RMS_vel for that sample
  n<-length(x)
}
# helper function for FFT profiling or to use as features
# Returns the Amplitude values
# n is the number of FFT features of interest
fft.profile <- function (dataset, n=5,bw=1000){
  amplitude <- lmod(dataset,bw)%>%log%>%multiply_by(10)
  # Ignore the 2nd half, which are complex conjugates of the 1st half,
  # and calculate the Mod (magnitude of each complex number)
  # amplitude <- Mod(fft.data[1:bw])%>%log
  # Calculate the frequencies
  frequencies <- seq(1, bw) # because 1 sec sampling
  # sorted  <- sort.int(amplitude[2:100], decreasing=TRUE, index.return=FALSE)
  # check for each element in top_peak is it within +/-2 Hz of any other in the array
  sorted  <- sort.int(amplitude[1:500], decreasing=TRUE, index.return=TRUE)
  top_idx <- sorted$ix[1:n] # indexes of the largest n components
  parp <- amplitude[top_idx]%>%sum%>%divide_by(length(amplitude[top_idx]))%>%
    divide_by(amplitude[-top_idx]%>%sum%>%divide_by(length(amplitude[-top_idx])))
  rpm_e   <- amplitude[1:100]%>%sort.int(.,decreasing = TRUE,index.return = TRUE)%$%ix[1:5]%>%getmode%>%divide_by(2)
  y       <- c(frequencies[top_idx[1:5]],parp,rpm_e)
  return (y) # convert indexes to frequencies
}
# Function to return the harmonic energy in each ofthe 5 freq bands from Brendan
# inputs needed are the RPM estimate, and the data sample
# outputs will be a  element array representing
# bal_en = energy +/- 30% centered on the 1st RPM harmonic
# aln_en = energy +/- 50% centered on the 2nd harmonic
# loo_en = energy from -50% 3rd to 5th + 50%
# brn_en = energy from -50% 5th to 25th
# gbx_en = energy above 25th
csi.profile <- function(dataset, bw=12000, parp_th=4){
  amplitude <- lmod(dataset,bw)%>%log%>%multiply_by(10)
  sorted  <- sort.int(amplitude[1:190], decreasing=TRUE, index.return=TRUE)
  top_idx <- sorted$ix[1:n] # indexes of the largest n components
  parp <- amplitude[top_idx]%>%sum%>%divide_by(length(amplitude[top_idx]))%>%
    divide_by(amplitude[-top_idx]%>%sum%>%divide_by(length(amplitude[-top_idx])))

  #plot.ts(amplitude[2:120]);abline(v=(top_idx-1),col="red")


  guess   <- amplitude[1:90]%>%sort.int(.,decreasing = TRUE,index.return = TRUE)%$%ix[1:5]
  rpm_e <- if ((parp >parp_th)&(guess[1] > 5)) guess[1] else 0

  bal_bw <- c(trunc(0.3*rpm_e):trunc(1.5*rpm_e))
  aln_bw <- c(trunc(1.5*rpm_e):trunc(2.5*rpm_e))
  loo_bw <- c(trunc(2.5*rpm_e):trunc(5.5*rpm_e))
  brn_bw <- c(trunc(5.5*rpm_e):trunc(25*rpm_e))
  gbx_bw <- c(trunc(25*rpm_e):bw)

  bal_en <- sum(amplitude[bal_bw])/(length(bal_bw))
  aln_en <- sum(amplitude[aln_bw])/(length(aln_bw))
  loo_en <- sum(amplitude[loo_bw])/(length(loo_bw))
  brn_en <- sum(amplitude[brn_bw])/(length(brn_bw))
  gbx_en <- sum(amplitude[gbx_bw])/(length(gbx_bw))
  y     <- c( bal_en,aln_en,loo_en,brn_en,gbx_en)
  return(y)
}


library(signal)

# Interquartile
#This quantifies the range between the 25th and 75th percentile and
# can be explored as a measure of variability of the underlying data.
# Percentiles (Yang et al., 2011): These quantify various positions of the underlying
# distribution can be employed to assist in the characterisation of extreme or anomalous points.
# Simple display of the Frequency content for a specific sample
#
# Peak-to-peak interval: This quantifies the time between peaks within a
# signal and can be used to assist in quantifying deviations
# from periodicity of an asset or component.
# Prior Knowledge is the RPM speed
##### Select which ACCOUNT
account= "boortmalt"
account= "analog"
account= "boortmalt"
account= "merck"
#account= "dairygold"
# Summary files first to generate TRENDS
root_dir  <- paste("C:/Users/dredmond/Desktop/DataSet/cache/data/",account, sep="")
tot_dir =list.dirs(root_dir,full.names=TRUE,  recursive=T)
use_dir_summ<-tot_dir[grep(tot_dir, pattern = "summary")]
#### RAW plots
# create a dummy version of infiles for summary
#on_off_State_files=list.files(use_dir,pattern="on_off_model.json",full.names = T,recursive = T)
infiles=list.files(use_dir_summ,pattern="*15........",full.names = T,recursive = T)
#
no_col  <- length(infiles)
print(no_col)
db_sum  <-  fromJSON(infiles[100])%>%as.data.frame()
db_sum %<>% slice(rep(1:n(),each=no_col))
#
sensor_list <- NULL
#readingin in JSON files
# each file si a single row entry of the summary statistics
for(j in 1:no_col)
{
  db_sum[j,] <-     fromJSON(infiles[j])%>%as.data.frame()
  sensor_list[j] <-fromJSON(infiles[j])$sensorid%>%as.character
}
# tidy up the DD-MM-YY formats before saving only the significant columns
db_sum$timestamp <- as.POSIXct(db_sum$timestamp,tz='UTC',origin="1970-01-01")
db_sum<-db_sum%>%select(-c(build,channel,lsb_per_g,sample_size,sample_rate,sensorid,vbat))%>%cbind(sensor_list)#%>%rename("sensorid",value)
write_rds(db_sum,paste0(root_dir,"db_summary.rds"),compress=c("none"))
sensor_id <- db_sum$sensor_list%>%as.factor%>%levels
#
for (i in seq_along(sensor_id))
  {
  y <-db_sum%>%filter(sensor_list==sensor_id[i])%>%filter(timestamp> "2018-07-02")%>%
    filter(stdDev>0.01)
  y.mah<- y%>%select(c(kurtosis,rms,stdDev,skewness))%>%covMcd%$%mah
  y<-cbind(y,y.mah)
  y%>%select(stdDev,rms,timestamp,sensor_list)%>%
    melt(id.vars=c("sensor_list","timestamp"))%>%
    ggplot(aes(timestamp,value,col="Conditions"))+
    geom_point(aes(timestamp,value,color=variable))+
    # stat_smooth(method="glm",se=FALSE)+
    facet_wrap(~variable,scales="free")+
    #theme_stata()+
    xlab("Trend of Condition Indicators ")+
    ggtitle((paste(sensor_id[i],"StdDev and RMS plots of select Condition Indicators")))
  ggsave(file=paste(sensor_id[i],"_timeseries.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
  #Trend for 2 indicators Mahalanobsis and StdDev
  y%>%filter(sensor_list==sensor_id[i] )%>%filter(timestamp>"2018-04-21")%>%
    select(timestamp,rms,kurtosis,crest,rms,stdDev,y.mah,skewness)%>%
    melt(id.vars=c("timestamp"))  %>%
    ggplot(aes(timestamp,value,col=variable))+
    theme_stata()+
    geom_point() +
    geom_smooth(method = 'loess', se=FALSE)+#, span=0.5, level=0.99)+
    #facet_wrap(~variable,scale="free") +
    ggtitle((paste(sensor_id[i],"TimeSeries plots of select Condition Indicators")))
  ggsave(file=paste(sensor_id[i],"_timeseries.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
  #Trend for Mahalanobsis and StdDev
  y%>%filter(sensor_list==sensor_id[i] )%>%filter(timestamp>"2018-04-21")%>%
    select(timestamp,y.mah)%>%
    melt(id.vars=c("timestamp"))  %>%
    ggplot(aes(timestamp,value,col=variable))+
    theme_stata()+
    geom_point() +
    geom_smooth(method = 'loess', se=FALSE)+#, span=0.5, level=0.99)+
    #facet_wrap(~variable,scale="free") +
    xlab("Trend of Mahalanobsis Indicators ")  +
    ggtitle((paste(sensor_id[i],"Mahalanobsis plots as a CI")))
  ggsave(file=paste(sensor_id[i],"_Mahalanobsis.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
  # individual STD-DEV models
  y%>%filter(sensor_list==sensor_id[i])%>%filter(timestamp>"2018-05-01")%>%
    select(timestamp,stdDev,sensor_list)%>%
    ggplot(aes(timestamp,stdDev,col=sensor_list))+
    theme_stata()+
    geom_point() +
    geom_smooth(method = 'loess',se=FALSE)+#, span=0.5, level=0.99)+
    xlab("StdDev Indicator ")  +
    ggtitle(paste(sensor_id[i],"plots of StdDev Indicator"))
  ggsave(file=paste(sensor_id[i],"_StdDev_timeseries.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
}
#

#
###deeper dive
# sensor_list == MCSVSZ0PV , MCSVSS02J, MCSVE6B8S, MCSVX5074
###
if (account== "merck") {
sub_set = c("MCSV10DS5","MCSVDOVBG") # Boiler
#for (i in seq_along(sub_set)) {

  db_sum%>%filter(sensorid==sub_set)%>%filter(timestamp>"2018-05-02")%>%filter(machine_on == 1) %>%
  ggplot(aes(timestamp,stdDev,colour=2))+
  geom_point()+
  geom_point(aes(timestamp,kurtosis,colour=4))+
  geom_smooth(method = "glm") +
  theme_stata()+
  theme(legend.position="none") +
  facet_wrap(~sensorid)+
  xlab("Trend of StdDev and kurtosis  ")+
  ylab("max value of sample measurement")+
  ggtitle("Steam Generator")
ggsave(file=paste(sensor_id[i],"_steam_generator.jpg"),dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
}
#
######################### DEEPER INVESTIGATION
#
#
root_dir  <- paste("C:/Users/dredmond/Desktop/DataSet/cache/data/",account, sep="")
tot_dir =list.dirs(root_dir,full.names=TRUE,  recursive=T)
use_dir_raw <-tot_dir[grep(tot_dir, pattern = "./csv_g")]
infiles     <-list.files(file.path(use_dir_raw),pattern="15.........*", full.names=TRUE)
########
freqData <- data.frame(1) ## Will store peak frequency data for each file
bw<- 500
j=1
i=1
h=1
n=5
bw= 1200
sensor_id <- NULL
#
#### RAW plots for each sensor
nncol= 9 #
#
raw_header <- c("sensor","timestamp","fft.1","fft.2","fft.3","fft.4","fft.5","parp","rpm_e")
db_fft <- matrix(NA,nrow=length(infiles),ncol= length(raw_header))
colnames(db_fft)    <-  raw_header%>%as.character()
#Set up the loop for each file within that sensor category
for (i in 1:length(infiles))
{
  db          <-  read_csv(infiles[i], col_name="X")
  timestamp   <-  infiles[i]%>%str_extract("15........")%>%as.numeric  #%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  ff_pro      <-   fft.profile(db,n,bw=800)
  db_fft[i,1]  <- str_extract(infiles[i], "MCS......")
  db_fft[i,2]   <- timestamp
  db_fft[i,3:9] <- ff_pro
}
#simple check to make sure everything is lined-up correctly
#
db_fft%<>%as.tibble
db_fft$timestamp    <- db_fft$timestamp%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
db_fft$sensor       <- db_fft$sensor%>%as.factor
db_fft[,3:9]        <- lapply(db_fft[,3:9],as.numeric,2)%>%lapply(.,round,1)
write_csv(path=paste0(root_dir,"_fft_raw.csv"),db_fft )

#db_fft<- read_csv(paste0(root_dir,"_fft_raw.csv") )
#
parp_th   <- 4 # Minimim Peak to Average ratio
sensor_id <- db_fft$sensor%>%as.factor%>%levels
for (i in seq_along(sensor_id))
{
  y <- db_fft%>%filter(parp>parp_th)%>%filter(sensor==sensor_id[i])
  y%>%melt(id.vars=c("timestamp","sensor"))%>%
    ggplot(aes(x=timestamp,y=value,fill=variable,col=variable))+
    geom_point()+
      #facet_wrap(~variable,scale="free_y")
      ylab("Frequency Hz")+
      xlab("Date ")+
      geom_smooth( method="glm",se=FALSE)+
      theme_stata()+
      #scale_y_continuous("Frequency Hz",sec.axis = sec_axis(~. *.1, name = "Amplitude dB"))+
      #theme_fivethirtyeight()+
      ggtitle(paste0(sensor_id[i]," Peak harmonics"))
      ggsave(file=paste0(account,sensor_id[i],"_Peak_harmonics.png"),
             dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
  }
#
#
#########################################################################
#########################################################################
#### Next Go To the extended_features section #######
#############Specific for MERCK
psd_sum <- function(amplitude,bw_array){  amplitude[bw_array]%>%sum%>%divide_by(length(bw_array))
}

extra.feature <- function(db,bw=1200,parp_th=5) {
  # Statistical features
  # Key frequencies
  # split the spectrum into 4 bins then calc. total energy / bin size
  # FFT calculation
  amplitude <- lmod(db,bw)%>%log%>%multiply_by(10)
  #
  sorted <- sort.int(amplitude, decreasing=TRUE, index.return=TRUE)
  top_idx <- sorted$ix[1:n] # indexes of the largest n components
  parp <- amplitude[top_idx]%>%sum%>%divide_by(length(amplitude[top_idx]))%>%
    divide_by(amplitude[-top_idx]%>%sum%>%divide_by(length(amplitude[-top_idx])))
  guess <- amplitude[1:100]%>%sort.int(.,decreasing = TRUE,index.return = TRUE)%$%ix[1:5]%>%divide_by(2)
  # need the rpm_e estimate to calculate the Fault freq with overlap
  rpm_e <- if ((parp >parp_th)&(guess[1] > 5)) guess[1] else 50000
  # dc to Sub-hamonics : Train fault, or imbalance
  # 1st harmonic to 2nd  Mis-align
  # 2nd to 4th  Bearing Spin
  # 3-6 X harmonic ~ Bearing Fault BPFO/i
  band_f <- if(is.na(rpm_e)) NA else c(2,1.05*rpm_e,2.05*rpm_e,6.05*rpm_e)%>%trunc()
  # Estimates the energy DENSITY with overlap
  vhf   <- psd_sum(amplitude=amplitude, bw_array= (band_f[4]-1):bw)
  hf    <- psd_sum(amplitude=amplitude, bw_array= (band_f[3]-1):(band_f[4]+1))
  mf    <- psd_sum(amplitude=amplitude, bw_array= (band_f[2]-1):(band_f[3]+1))
  lf    <- psd_sum(amplitude=amplitude, bw_array= (band_f[1])  :(band_f[2]+1))
  #
  bpfo  <-  psd_sum(amplitude=amplitude, bw_array=(floor(3.62*rpm_e-1):floor( 3.62 * rpm_e+1)))
  bpfi  <-  psd_sum(amplitude=amplitude, bw_array=(floor(5.38*rpm_e-1):floor( 5.38 * rpm_e+1)))
  ftf   <-  psd_sum(amplitude=amplitude, bw_array=(floor(0.40*rpm_e-1):floor( 0.40 * rpm_e+1)))
  bsf   <-  psd_sum(amplitude=amplitude, bw_array=(floor(2.46*rpm_e-1):floor( 2.46 * rpm_e+1)))
  #
  # Statistical features
  stat_feat       <- c(quantile(db$X, names=FALSE), mean(db$X), sd(db$X), skewness(db$X), kurtosis(db$X), rms_x(db$X))
  fault_freq_amp  <- c(ftf,bpfi,bpfo,bsf)
  psd             <- c(vhf,hf,mf,lf)  #psd<- c(sum(vhf)%>%divide_by(length(vhf)), sum(hf)%>%divide_by(length(hf)), sum(mf)%>%divide_by(length(mf)), sum(lf)%>%divide_by(length(lf)))
  features        <- c(stat_feat, fault_freq_amp, psd)
  #  Returns the following ....
  #  ff_pro   <- "fft[1:5]","parp","rpm_e" ;;7
  # stat_feat <- "Min", "Qu.1", "Median", "Qu.3", "Max", "Mean", "stdDev", "Skew", "Kurt", "RMS", ;; 10
  # fault_freq_amp <-   "FTF", "BPFI", "BPFO.x", "BSF.x",  ;; 4
  # Power Spec. Density <- "psd[1:4]" ;; 4
  return(features)
}
#
extra_db  <- matrix(nrow=length(infiles), ncol=32)
#Set up the loop for each file within that sensor category
for (i in 1:length(infiles) ) {
  db            <- read_csv(infiles[i], col_name="X")
  ff_pro        <- fft.profile(db,n,bw=bw)
  ex_pro      <- extra.feature(db,bw=1000,parp_th = parp_th)
  csi_pro     <- csi.profile(db,bw=10000)
  extra_db[i,2]     <-  infiles[i]%>%str_extract("15........")%>%as.numeric
  extra_db[i,3:9]   <- ff_pro%>%as.vector()
  extra_db[i,10:27] <- ex_pro%>%as.vector()
  extra_db[i,28:32] <- csi_pro%>%as.vector()
  extra_db[i,1]     <- str_extract(infiles[i], "MCS......")
   }
ex.cnames <- c("sensor", "timestamp","fft.1","fft.2","fft.3","fft.4","fft.5","parp","rpm_e",
               "Min", "Qu.1", "Median", "Qu.3", "Max",
               "Mean", "stdDev", "Skew", "Kurt", "RMS",
               "FTF", "BPFI", "BPFO", "BSF",
               "VHF.psd","HF.psd","MF.psd","LF.psd",
               "bal_en","aln_en","loo_en","brn_en","gbx_en")

#
extra_db%<>%as.tibble
colnames(extra_db)    <-  ex.cnames%>%as.character()
extra_db$timestamp    <- extra_db$timestamp%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
extra_db$sensor       <- extra_db$sensor%>%as.factor
extra_db[,3:32]       <- lapply(extra_db[,3:32],as.numeric,2)%>%lapply(.,round,2)
write_csv(path=paste0(root_dir,"_extra"),extra_db )

#extra_db$timestamp   <- seq(1,length(infiles)*4800,4800)%>%as.POSIXct(tz="UTC",origin="2018-05-01")


#
sensor_id <- extra_db$sensor%>%as.factor%>%levels
  for (i in seq_along(sensor_id))
  {
    # Harmonics
    extra_db%>%filter(parp>parp_th)%>%filter(sensor==sensor_id[i])%>%
      mutate(rpm_x=getmode(rpm_e))%>%
      select(c(fft.1,fft.2,fft.3,fft.4,fft.5,rpm_x,timestamp,sensor))%>%
      melt(id.vars=c("timestamp","sensor"))%>%
      ggplot(aes(x=timestamp,y=value,fill=variable,col=variable))+
      geom_point()+
      ylab("Frequency Hz")+
      xlab("Date ")+
      #geom_smooth( method="glm",se=FALSE)+
      theme_stata()+
      ggtitle(paste0(sensor_id[i]," Harmonics_freatures"))
      ggsave(file=paste0(account,sensor_id[i],"_Harmonics_features.png"),
           dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
    #
    ## Statistical features
    extra_db%>%filter(parp>parp_th)%>%filter(sensor==sensor_id[i])%>%
      select(c(Max,Mean,stdDev,timestamp,sensor))%>%
      melt(id.vars=c("timestamp","sensor"))%>%
      ggplot(aes(x=timestamp,y=value,fill=variable,col=variable))+
      geom_point()+
      ylab(" Summary values")+
      xlab("Date ")+
     # geom_smooth( method="glm",se=FALSE)+
      theme_stata()+
      ggtitle(paste0(sensor_id[i]," Statistical_freatures"))
      ggsave(file=paste0(account,sensor_id[i],"_Statistical_features.png"),
           dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
    #
    ## PSD features
    extra_db%>%filter(parp>parp_th)%>%filter(sensor==sensor_id[i])%>%
      select(c(VHF.psd,HF.psd,MF.psd,LF.psd, timestamp,sensor))%>%
      melt(id.vars=c("timestamp","sensor"))%>%
      ggplot(aes(x=timestamp,y=value,fill=variable,col=variable))+
      geom_point()+
      ylab("Amplitude dB")+
      xlab("Date ")+
      geom_smooth( method="glm",se=FALSE)+
      theme_stata()+
      ggtitle(paste0(sensor_id[i]," PSD_freatures with Bandwidth set to",bw,"Hz"))
    ggsave(file=paste0(account,sensor_id[i],"_PSD_features.png"),
           dpi=900,path=root_dir,width=25,height=22,units=c("cm"))

    ## Fault freq features
    extra_db%>%filter(parp>parp_th)%>%filter(sensor==sensor_id[i])%>%
      select(c(FTF, BPFI, BPFO, BSF,  timestamp,sensor))%>%
      melt(id.vars=c("timestamp","sensor"))%>%
      ggplot(aes(x=timestamp,y=value,fill=variable,col=variable))+
      geom_point()+
      ylab("Amplitude dB")+
      xlab("Date ")+
      geom_smooth( method="glm",se=FALSE)+
      theme_stata()+
      ggtitle(paste0(sensor_id[i]," Fault_freq with Bandwidth set to",bw,"Hz"))
    ggsave(file=paste0(account,sensor_id[i],"_Fault_freq.png"),
           dpi=900,path=root_dir,width=25,height=22,units=c("cm"))

    extra_db%>%select(c(parp,gbx_en,brn_en,loo_en,aln_en,bal_en,rpm_e,timestamp))%>%
      filter(parp>parp_th)%>% melt(id.vars=c("timestamp"))%>%
      ggplot(aes(x=timestamp,y=value,fill=variable,col=variable))+
      geom_point() +
      ylab("Amplitude dB /root(Hz)")+
      xlab("Date ")+
      geom_smooth( method="gam",se=TRUE)+
      facet_wrap(~variable,scales="free") +
      theme_linedraw()+
       ggtitle(paste0(sensor_id[i]," Coarse fault indicators with 1st harmonic"))
    ggsave(file=paste0(account,sensor_id[i],"_Coarse_Fault_freq.png"),
           dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
    #

  }
    ## Guideline features

#
# lets begin to look for a regular regression on ALL the features for each sensor
#

# Then look for clusters - initial for each sensor, but then we could try across sensor_id
#
# Check the strength of the feature correlation
#
# Isolate one sensor with a fault condition, and model this condition using all the indicators
# see what does ti take to overfit
#
# explore if we could use TIMESERIES forecasting on the fault freq. to identify trends

#install.packages("corrplot")
#
#

library(corrplot)

pc <- princomp(~., data=db_sum[,-c(13:15)] , cor = TRUE, score = TRUE)
summary(pc)
plot(pc)
plot(pc, type = "l")
pc$loadings   ##influence or eigen values
#
db_cor<- cor(db_sum[,-c(13:15)], use="pairwise.complete.obs")
corrplot(db_cor,method = "color",order = "FPC", tl.srt = 45, mar = c(0.01, 0.01, 0.01, 0.01))
ggsave(file=paste0(account,sensor_id[i],"_corrplot.jpg"),
       dpi=900,path=root_dir,width=25,height=22,units=c("cm"))
# lets do a little PCA test
y.clust<-db_sum[,-c(13:15)]%>%scale(.,scale = FALSE)
#select(-c(sensor_list,timestamp,clearance,impulse,min,max,mean,variance,shape))%>%scale(.,scale=FALSE)
library(BayesLCA)
library(mclust)
y.mod <- Mclust(y.clust)
summary(y.mod)
plot.Mclust(y.mod)







db_sum%>%select(c(timestamp,stdDev,sensor_list))%>%melt(id.vars=c("sensor_list","timestamp"))%>%
  ggplot(aes(x=timestamp,value,color=sensor_list))+geom_point()+
  facet_wrap(~variable)


