#Purpose of this script is
## Libraries and helper functions
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
library(e1071)
library(readr)
library(ggthemes)
#
options(mc.cores = parallel::detectCores())
# Functions
#
rms_x     <- function(xx){( sqrt(sum(xx^2)/length(xx) ) )}
rms_diff  <- function(xx){( sqrt(sum( )))}
crest_x   <- function(xx){( max(abs(xx)/(sqrt(sum(xx^2))/length(xx)))) }

mah   <-function(xx,mah_offset=0) {mahalanobis(xx,colMeans(xx),cov(xx))-mah_offset }
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#####----
# helper function for FFT profiling or to use as feature
##### Select which ACCOUNT
account= "analog"
account= "boortmalt"
account= "merck"
#account= "dairygold"
#account= "analog"

# Summary files first to generate TRENDS
db_raw  <-read_csv(paste0("C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/dataframes_csv/",account,"_raw.csv"))
# pick a sensor
sensor_list <- db_raw%>%colnames%>%str_extract("MCSV.....")%>%as.factor
timestamp   <- db_raw%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
# split the sensor_list into 4


db_sum  <-read_csv(paste0("C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/dataframes_csv/",account,"_summary.csv"))
db_sum$timestamp%<>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
db_sum$sensorid%<>%as.factor()
db_sum$rmsvel  %<>%as.numeric()
db_sum$rpm  %<>%as.numeric()

sensor_list <- db_sum$sensorid%>%as.factor%>%levels
list_a <- sensor_list


peak_freq = 100
scale_factor  = 1/ (2 * pi * peak_freq * sqrt(2))*9.8*1000
velocity_cal <- (1561.866*Acc)/freq.


list_c <-c("MCSV8LIAW",  "MCSVE76FY" , "MCSVEJ945", "MCSVEYY7VG", "MCSVPNH5Z","MCSVW3B28","MCSVZFDRK","MCSVXLOCE" )
list_a <-c("MCSV82BP4","MCSVDOVBG")

list_a <-c("MCSV8LIAW",  "MCSVE76FY" , "MCSVEJ945", "MCSVEYY7VG", "MCSVPNH5Z","MCSVW3B28","MCSVZFDRK","MCSVXLOCE" )
#merck list
#list_a <- c( "MCSV10DS5", "MCSV82BP4", "MCSVD0VBG", "MCSVDOVBG", "MCSVE6B8S" ,"MCSVFS05V", "MCSVLXBBD", "MCSVPCDEW", "MCSVRQ00M", "MCSVSS02J", "MCSVSZ0PV", "MCSVX5074")
list_a <- c(  "MCSVD0VBG", "MCSVDOVBG" ,"MCSVLXBBD" )

vel_only<- db_sum%>%filter(machine_on ==1)%>%
  filter(rmsvel < 9.50)%>%filter(rmsvel >1)%>%
  filter(timestamp >"2018-10-14 00:00:00 UTC")%>%
  filter(sensorid == list_a)%>%
  select(rmsvel,rpm,sensorid,timestamp) %>% na.omit


vel_only%>%select(-timestamp)%>%aggregate(.~sensorid,data=.,FUN="sd")

res<- lm(timestamp~ rmsvel,data=vel_only )
rstudent(res)
lev <- hat(model.matrix(res))
plot(lev)

rpm_est<-vel_only%>%select(rpm,sensorid)%>%aggregate(.~sensorid,data=.,FUN="getmode")
print(rpm_est)
ggplot(vel_only,aes(x=timestamp,y=rmsvel,col=sensorid))+geom_point() +
  facet_wrap(~sensorid,scales="free")+
  #stat_smooth(method="glm",se=F)+
  ggtitle("Trending RMS velocity for Suspicious assets")



# Boortmalt
fan_FO_8011 <-c("MCSV1EJJU ","MCSV7MVF1","MCSVDB6UQ")

vel_only<- db_sum%>%
  filter(machine_on ==1)%>%
  filter(timestamp >"2018-09-13 00:00:00 UTC")%>%
  filter(sensorid == fan_FO_8011 )%>%
  select(rpm,sensorid,timestamp) %>% na.omit
vel_only%>%aggregate(.~sensorid,data=.,FUN="getmode")


### RAW aADC files
db_raw  <-read_csv(paste0("C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/dataframes_csv/",account,"_raw.csv"))
# pick a sensor
sensor_list <- db_raw%>%colnames%>%str_extract("MCSV.....")%>%as.factor
timestamp   <- db_raw%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
# split the sensor_list into 4


# This data looks useful
# Next read in the RAW_ADC datafiles
use_dir_raw  <- paste0("C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/dataframes_csv")
db_adc_files <-list.files(file.path(use_dir_raw),pattern="*_48k.csv", full.names=TRUE)

account = "boortmalt"
db_brt_files <- db_adc_files[which(db_adc_files%>%str_extract(account)==account)]
upwork_dir  <- paste0("C:/Users/dredmond/Documents/CbM/Deployments/UpWork")

i=1
for ( i in i : length(db_brt_files))
  {
  db_brt_raw  <- read_csv(db_brt_files[i])
  timestamp   <- db_brt_raw%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  sensor_id <- db_brt_raw%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  idx <- which(timestamp > "2018-10-30")
  write_csv(db_brt_raw[,idx],file.path(upwork_dir, paste0(sensor_id,"_2018-10-30.csv")) )
}


#----------------------------------------------------
lmod<-function(x,bw=5000){scale(x,center=T,scale=FALSE)%>%fft%>%Mod}
q <-nrow(db_brt_raw)
q_4 <- q/4
db.fft <- matrix(nrow=,ncol=ncol(db_brt_raw))
 #
for( j in 1:ncol(db_brt_raw))     {
  temp    <- cbind(db_brt_raw[1:q_4,j], db_brt_raw[(q_4+1):(q_4*2),j], db_brt_raw[(q_4*2+1):(q_4*3),j], db_brt_raw[(q_4*3+1):q,j])%>%apply(.,1,mean)%>%as.tibble
  db.fft  <-  temp%>%lmod%>%divide_by(1561)%>%as_tibble()
  all_fft <-  temp%>%lmod
  peak_fft <- whichall_fft[1:99]
  peak_acc <- max(db.fft$value[1:99])

}

band_wd = 4000
frequency<- seq(1 : band_wd)
harmonics <- peak_fft*(seq(from=1,to=20, by=2))
fft_plt <- cbind(frequency,db.fft[1:band_wd,1])
j
ggplot(fft_plt,aes(frequency,value))+geom_line()+
  geom_vline(xintercept=harmonics,col="green")+
  geom_hline(yintercept=peak_acc,col="blue")+
  ylab("mm/sec")+
  ggtitle(paste0("Spectrum plot for ",sensor_id," showing Odd harmonics ",peak_fft*60," RPM"))





temp[1:500,]%>%aggregate(.,data=.,FUN="getmode")






psd_features <- function(db)
{
  amplitude        <- lmod(db)## Log of the Mod with filtering
  y_psd <- seq(1,9)
  for( i in 1 :9) {
    f           <- amplitude[freq_interval[i]:freq_interval[(i+1)]]
    y_psd[i]    <- amplitude[f]%>%psd_fn
  }
  return(y_psd)
}
# --
lmod    <- function(x,bw=20000){return(scale(x,center=T,scale=FALSE)%>%fft%>%Mod)}
psd_fn  <- function(amplitude){return(amplitude%>%na.omit%>%rms)}
stat_features <-function(db) {stat_feat<-c(  mean(db), sd(db), e1071::skewness(db), e1071::kurtosis(db,na.rm=TRUE,type=2), rms_x(db), crest_x(db) )
return(stat_feat)
}

# --
#
freq_interval<-c(0,25,50,100,200,350,500,800,1200,3000)

db_freq_name<-c("f1","f2","f3","f4","f5","f6","f7","f8","f9")
db_stat_name<-c("mean","stdDev","skew","kurt","rms","crest")

db_psd  <- matrix(nrow=length(db_raw),ncol=9)
db_stat <- matrix(nrow=length(db_raw),ncol=6)

for (j in 1: length(db_raw)) {
     db          <- db_raw[,j]
     db_psd[j,]  <- apply(db,2,psd_features) %>%t   #%>%as.tibble()
     db_stat[j,] <- apply(db,2,stat_features)%>%t   #%>%as.tibble()
}
#
colnames(db_psd)  <- db_freq_name
colnames(db_stat) <- db_stat_name
db_all  <- cbind(db_psd,db_stat)
#
db_all<-cbind(db_psd,db_stat,sensor_list,timestamp)%>%as.tibble
db_all$sensor_list%<>% as.factor
db_all$timestamp%<>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
#
db_all%>% dplyr::select(c("kurt","rms","f1","f2","f3","f4","f5","timestamp","sensor_list")) %>%
  melt(id.vars=c("timestamp","sensor_list")) %>%
  ggplot(aes(x=timestamp,y=value,colour=variable))+geom_point() +
  facet_wrap(~sensor_list,scales = "free")  +
  geom_smooth(method="loess",se=F)
#







################  Visualisation for DESKTOP
glimpse(db_sum)

db_sum%>%select(machine_on,timestamp,stdDev,kurtosis,sensorid)%>%
  filter(sensorid=="MCSV10DS5")%>%
  filter(machine_on==1)%>%filter(timestamp>"2018-10-08")%>%
  filter(machine_on==1)%>%filter(timestamp<"2018-10-18")%>%
  ggplot(aes(x=timestamp,y=stdDev,col=sensorid,fill=sensorid))+
  geom_line()+ geom_point(colour="purple",size=0.75)+
  xlab("Measurement Date stamp")+
  ylab("Standard Deviation")+
  theme_stata()+
  theme(legend.position="none") +
  #annotate("rect", xmin = "2018-10-14 09:00:00 UTC", xmax = "2018-10-15 09:00:00 UTC", ymin = 0.1, ymax = 0.8, alpha = .2, col="pink")+
  #geom_rect(aes(xmin=159683.438, xmax=159684.186, ymin=0, ymax=Inf))+
  ggtitle("Trending Summary Statistic with Anomoly Detection")


db_fft<-read_csv

