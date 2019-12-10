source("./CbM/r-code/basics.R")
account="merck"
#pick out the data files
source("./CbM/r-code/data_str.R")
# Load data


psd_energy <- function(xx){
  psd_a <- sum(xx[2:50])    / length(xx[0:50])
  psd_b <- sum(xx[45:150])  / length(xx[45:150])
  psd_c <- sum(xx[120:400])  / length(xx[120:400])
  psd_d <- sum(xx[320:800]) / length(xx[320:800])
  psd_e <- sum(xx[720:6000]) / length(xx[720:6000])
  y <- cbind(psd_a,psd_b,psd_c,psd_d,psd_e)
  return(y)
}
on_off     <- function(xx){
  thres    <- (max(xx)-min(xx))/ 3
  yy      <-  which(xx > thres, arr.ind=TRUE)
  return(yy)
}


db_ffts     <- hdf_db(hdf_ffts[1])
db_samps    <- hdf_db(hdf_samps[1])
timestamp   <- db_ffts%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")

psd_db    <- apply(db_ffts[,-1],2,psd_energy)%>%t
colnames(psd_db) <- c("v_low_2to50","low45to150","mid120to400","high320to800","v_high720to6000")

psd_db<-cbind(timestamp[-1],psd_db)


psd_db%>%as.tibble%>%na.omit()%>%
#  select(timestamp,v_low,low,mid,high,v_high)%>%
  select(timestamp,v_low,low,mid)%>%
  filter(timestamp < "2018-12-01")%>%
  filter(timestamp > "2018-09-10")%>%
  melt(id.vars=("timestamp")) %>%
  ggplot(aes(timestamp,value,colour=variable))+geom_line()+
  facet_wrap(~variable)


             ?






for (i in 1:length(hdf_fft)){
  db_fft   <-  0
  db_std   <-  0
  db_psd   <-  0
  psd_plt  <-  0
  db_fft   <- hdf_db(hdf_fft[i])
  db_samps <- hdf_db(hdf_samp[i])
  db_std   <- apply(db_samps, 2, sd)%>%as.vector()%>%as.numeric()
  db_thre  <- db_std > quantile(db_std)[3]%>%as.logical()

  db_psd   <- apply(db_fft, 2, psd_energy)%>%t
  timestamp   <- db_samps%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  #timestamp   <- timestamp[-1]
  #
  sensor_list  <- db_fft%>%colnames%>%str_extract("MCS......")%>%as.factor%>%levels()
  #
  psd_plt    <- cbind(db_psd[-1,],db_thre)
  colnames(psd_plt) <- c( "vlow","low","mid","high","vhigh", "stdDev")
  psd_plt$stdDev%<>%as.factor()
  #
  psd_plt%>%mutate(timestamp)%>%
    filter(timestamp>"2018-09-12")%>%filter(stdDev == 1)%>%
    melt(id.vars=c("timestamp", "stdDev"))  %>% na.omit()%>%
    ggplot(aes(timestamp,value,col=variable))  + geom_point() +
    stat_smooth(method = "loess",se=F)  +
    ggtitle(paste(sensor_id,i,"PSD trends"))

  ggsave(paste(sensor_id,i,"_psd.pdf"),path=hdf_dir)

}

## build a GBM model on the std Dev threshold




index <- createDataPartition(psd_plt$stdDev, p = 0.7, list = FALSE)
train_data <- psd_plt[index, ]
test_data  <- psd_plt[-index, ]

model_gbm <-  train(stdDev ~ ., data = train_data, method = "gbm", preProcess = c("scale", "center"),  trControl = trainControl(method = "repeatedcv",  number = 5,repeats = 3,verboseIter = FALSE), verbose = 0)
model_gbm

confusionMatrix(data = predict(model_gbm,test_data),reference = test_data[,6]  )






  db_adc<- hdf_db(hdf_samp[i])
  db_psd <- apply(db_fft, 2,psd_energy)%>%t
  colnames(db_psd) <- c("vlow","low","mid","high","vhigh")
  timestamp   <- db_fft%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")

  sensor_id  <- db_fft%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  timestamp   <- db_fft%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  n_col <- ncol(db_fft)
  p_col <- n_col-200
# rpm_id      <- rpm_est[which(rpm_est$sensorid==sensor_id),]%$%rpm
# harmonics    <- seq(1:20)*rpm_id/60
  db_fft_plt<- cbind(db_fft[,1],db_fft[,c(p_col:n_col)])
  colnames(db_fft_plt) <- c("freq", timestamp[p_col:n_col])
  band_width  <- 550
#db_fft_plt%>%as.tibble()%>%melt(id.vars=c("freq"))%>%filter(freq<band_width)%>%
# ggplot(aes(freq,value,col=variable))+geom_line()+
#geom_vline(xintercept = rpm_s,col="red",linetype="dashed")+
#    facet_wrap(~variable)+
#  ggtitle(paste(sensor_id,"spectrum  plot of last 6 measurements"))
#  ggsave(file=paste(sensor_id,"spectra_plot_last_6.pdf"),dpi=900,path=hdf_dir,width=1080,height=182,units=c("mm"))
# now plot the residuals
  band_width  <- 550
  in_band_psd <- apply(db_fft_plt[0:band_width,-1],2,rms)%>%as.vector
  rpm_s  <- apply(db_fft_plt,2,which.max)-2
  db_fft_res <- db_fft_plt
  for (j in 2:length(rpm_s)){
    db_fft_res[,j]<-replace(db_fft_plt[,j],rpm_s[j]:(rpm_s[j]+5),0)
  }
  resid_psd <- apply(db_fft_res[0:band_width,-1],2,rms)%>%as.vector
  time_x    <- timestamp[p_col:n_col]%>%as.vector
  res_df    <- rbind(resid_psd,in_band_psd,time_x)%>%t%>%as.tibble()
  res_df$time_x%<>%as.POSIXct(tz='UTC',origin="1970-01-01")
  res_df%>%melt(id.vars=c("time_x"))%>%
    ggplot(aes(time_x,value,col=variable)) +geom_point()+
    scale_x_datetime()+
    ggtitle("Tracing PSD and Residual PSD for",sensor_id)
  #
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
db_adc<-hdf_db(hdf_account)
db_sum$timestamp%<>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
db_sum$sensorid%<>%as.factor()
db_sum$rmsvel  %<>%as.numeric()
db_sum$filt_rms  %<>%as.numeric()
db_sum$rpm  %<>%as.numeric()
sensor_list <- db_sum$sensorid%>%as.factor%>%levels
list_a <- sensor_list
rpm_est<-db_sum%>%filter(machine_on==1)%>%select(rpm,sensorid)%>%aggregate(.~sensorid,data=.,FUN="getmode"



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
###############
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
  adc_samps <- re
  adc_account[i])
  #glimpse(db_samps)
  timestamp   <- db_samps%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  sensor_id   <- db_samps%>%colnames%>%str_extract("MCSV.....")%>%as.factor%>%levels()
  time_axis   <- seq(1,nrow(db_samps))/48000
  db_ssamps   <-  cbind(time_axis[0:12000],db_samps[0:12000,])
  #db_1<- db_samps[,2707]
  s_ncol      <- ncol(db_samps)
  p_col       <- s_ncol-300
  db_sub      <- db_samps[,p_col:s_ncol]
  timestamp   <- db_sub%>%colnames%>%str_extract("15........")%>%as.numeric%>%as.POSIXct(tz='UTC',origin="1970-01-01")
  rms_dif     <-apply(db_sub,2,rms_d)%>%as.vector()
  std_vec     <-apply(db_sub,2,sd)%>%as.vector()
  krt_vec     <-apply(db_sub,2,kurtosis)%>%as.vector()
  df          <- cbind(timestamp,rms_dif,std_vec,krt_vec)%>%as.tibble()
  df$timestamp%<>%as.POSIXct(tz='UTC',origin="1970-01-01")
  df%>% melt(id.vars=c("timestamp"))%>%
    ggplot(aes(timestamp,value,col=variable))+geom_line()+
    scale_x_datetime()+
    ggtitle(paste(sensor_id,"trending features"))
  ggsave(file=paste(sensor_id,"tending_features.pdf"),plot=last_plot(),dpi=900,path=raw_dir,width=1080,height=182,units=c("mm"))
}
# facet_wrap(~variable)
# When did it happen?
timestamp[which.max(std_vec) ]
timestamp[which.max(std_vec[0:481])]
# Prior to that when was the last peak
#
#install.packages(c("rgenoud","parallel","furrr","tsibble","forecast"))

library(tidyverse)
library(forecast)
library(rgenoud)
library(parallel)
library(lubridate)
library(furrr)
library(tsibble)
#library(brotools)
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

avia_clean_monthly <- read_csv("https://raw.githubusercontent.com/b-rodrigues/avia_par_lu/master/avia_clean_monthy.csv")


order_list <- list("p" = seq(0, 3), "d" = seq(0, 2),"q" = seq(0, 3)) %>%
  cross() %>% map(lift(c))
  list("p" = seq(0, 3), "d" = seq(0, 2), "q" = seq(0, 3)) %>%
  cross() %>%
  head(3)


avia_clean_train <- avia_clean_monthly %>%
  select(date, passengers) %>%
  filter(year(date) < 2013) %>%
  group_by(date) %>%
  summarise(total_passengers = sum(passengers)) %>%
  pull(total_passengers) %>%
  ts(., frequency = 12, start = c(2005, 1))

avia_clean_validation <- avia_clean_monthly %>%
  select(date, passengers) %>%
  filter(between(year(date), 2013, 2016)) %>%
  group_by(date) %>%
  summarise(total_passengers = sum(passengers)) %>%
  pull(total_passengers) %>%
  ts(., frequency = 12, start = c(2013, 1))

avia_clean_test <- avia_clean_monthly %>%
  select(date, passengers) %>%
  filter(year(date) >= 2016) %>%
  group_by(date) %>%
  summarise(total_passengers = sum(passengers)) %>%
  pull(total_passengers) %>%
  ts(., frequency = 12, start = c(2016, 1))



logged_test_data <- ihs(avia_clean_test)
logged_validation_data <- ihs(avia_clean_validation)
logged_train_data <- ihs(avia_clean_train)

cost_function_rmse <- function(param, train_data, validation_data, forecast_periods){
  order <- param[1:3]
  season <- c(param[4:6], 12)
  model <- purrr::possibly(arima, otherwise = NULL)(x = train_data, order = order,seasonal = season,method = "ML")
  if(is.null(model)){
    return(9999999)
  } else {
    forecast_model <- forecast::forecast(model, h = forecast_periods)
    point_forecast <- forecast_model$mean
    sqrt(mean(point_forecast - validation_data) ** 2)
  }
}

library(forecast)
starting_model <- auto.arima(logged_train_data)
summary(starting_model)
cost_function_rmse(c(1, 0, 2, 2, 1, 0),
                   train_data = logged_train_data,
                   validation_data = logged_validation_data,
                   forecast_periods = 65)



domains <- matrix(c(0, 3, 0, 2, 0, 3, 0, 3, 0, 2, 0, 3), byrow = TRUE, ncol = 2)


library(parallel)


auto_arima_rmse <- genoud(cost_function_rmse,
                          nvars = 6,
                          data.type.int = TRUE,
                          starting.values = c(1, 0, 2, 2, 1, 0), # <- from auto.arima
                          Domains = domains,
                          #cluster = cl,
                          train_data = logged_train_data,
                          validation_data = logged_validation_data,
                          forecast_periods = length(logged_validation_data),
                          hard.generation.limit = TRUE)





best_model_rmse <- arima(logged_train_data, order = auto_arima_rmse$par[1:3],
                         season = list(order = auto_arima_rmse$par[4:6], period = 12),   method = "ML")
summary(best_model_rmse)


best_model_rmse_forecast <- forecast::forecast(best_model_rmse, h = 65)

best_model_rmse_forecast <- to_tibble(best_model_rmse_forecast)
## Joining, by = "date"
## Joining, by = "date"
starting_model_forecast <- forecast(starting_model, h = 65)

starting_model_forecast <- to_tibble(starting_model_forecast)
## Joining, by = "date"
## Joining, by = "date"
plo

avia_clean_monthly %>%
  group_by(date) %>%
  summarise(total = sum(passengers)) %>%
  mutate(total_ihs = ihs(total)) %>%
  ggplot() +
  ggtitle("Minimization of RMSE") +
  geom_line(aes(y = total_ihs, x = date), colour = "#82518c") +
  scale_x_date(date_breaks = "1 year", date_labels = "%m-%Y") +
  geom_ribbon(data = best_model_rmse_forecast, aes(x = date, ymin = lower95, ymax = upper95),
              fill = "#666018", alpha = 0.2) +
  geom_line(data = best_model_rmse_forecast, aes(x = date, y = point_estimate),
            linetype = 2, colour = "#8e9d98") +
  geom_ribbon(data = starting_model_forecast, aes(x = date, ymin = lower95, ymax = upper95),
              fill = "#98431e", alpha = 0.2) +
  geom_line(data = starting_model_forecast, aes(x = date, y = point_estimate),
            linetype = 2, colour = "#a53031") +
  theme_blog()




cost_function_bic <- function(param, train_data, validation_data, forecast_periods){
  order <- param[1:3]
  season <- c(param[4:6], 12)
  model <- purrr::possibly(arima, otherwise = NULL)(x = train_data, order = order,
                                                    seasonal = season,
                                                    method = "ML")
  if(is.null(model)){
    return(9999999)
  } else {
    BIC(model)
  }
}




cost_function_bic(c(1, 0, 2, 2, 1, 0),
                  train_data = logged_train_data,
                  validation_data = logged_validation_data,
                  forecast_periods = 65)
## [1] -184.6397
Let the genetic algorithm run again:

  cl <- makePSOCKcluster(8)
clusterExport(cl, c('logged_train_data', 'logged_validation_data'))

tic <- Sys.time()

auto_arima_bic <- genoud(cost_function_bic,
                         nvars = 6,
                         data.type.int = TRUE,
                         starting.values = c(1, 0, 2, 2, 1, 0), # <- from auto.arima
                         Domains = domains,
                         cluster = cl,
                         train_data = logged_train_data,
                         validation_data = logged_validation_data,
                         forecast_periods = length(logged_validation_data),
                         hard.generation.limit = TRUE)
toc_bic <- Sys.time() - tic
This time, it took 6 minutes, a bit slower than before. Let's take a look at the solution:

#  auto_arima_bic




########################

library(dlm)
library(tidyverse) # for plotting
# set parameters
mu_z <- 0.1
mu_y <- 0.5
rho1 <- 0.9
sig_z <- 1
sig_u1 <- 10
sig_u2 <- 0.1

GG <- matrix(c(1+rho1,-rho1,mu_z,1,0,0,0,0,1),nrow=3,byrow=TRUE)
FF <- matrix(c(1,0,mu_y,0,1,0,0,0,1),nrow=3,byrow=TRUE)
m0 <- matrix(c(0,0,0), nrow=3)
C0 <- matrix(rep(0,9),nrow=3)
W  <- matrix(c(sig_z,0,0,
              0,0,0,
              0,0,0),nrow=3,byrow=TRUE)

V <- matrix(c(sig_u1,0,0,
              0,sig_u2,0,
              0,0,0),nrow=3,byrow=TRUE)

my_dlm <- dlm(FF=FF,V=V,GG=GG,W=W,V=V,m0=c(0,0,1),C0=C0)

set.seed(20180610)  #set seed for reproducible results
y1 <- dlmForecast(my_dlm, nAhead=120,sampleNew=1)
df <- data.frame(y=y1$newObs, z=y1$newStates)
df$id <- seq.int(nrow(df))
ggplot(data=df, aes(x=id,y=y.1))+
  geom_line()+
  geom_line(linetype=2,aes(y=z.1),color="red")+
  labs(x="time",y="Z", title="Simulated dynamic linear model", subtitle="black solid line observation y1\nred dotted line fundamental z")


my_dlmfc <- function(par=c(rho1,sig_z,sig_u1,sig_u2,mu_z,mu_y)){
  rho1 = par[1]
  sig_z = par[2]
  sig_u1 = par[3]
  sig_u2 = par[4]
  mu_z = par[5]
  mu_y = par[6]
  GG <- matrix(c(1+rho1,-rho1,mu_z,1,0,0,0,0,1),nrow=3,byrow=TRUE)
  FF <- matrix(c(1,0,mu_y,0,1,0,0,0,1),nrow=3,byrow=TRUE)
  m0 <- matrix(c(0,0,0), nrow=3)
  C0 <- matrix(rep(0,9),nrow=3)
  W <- matrix(c(sig_z,0,0,
                0,0,0,
                0,0,0),nrow=3,byrow=TRUE)
  V <- matrix(c(sig_u1,0,0,
                0,sig_u2,0,
                0,0,0),nrow=3,byrow=TRUE)
  my_dlm3 <- dlm(FF=FF,V=V,GG=GG,W=W,V=V,m0=c(0,0,1),C0=C0)
  return(my_dlm3)
}



my_mle <- dlmMLE(matrix(unlist(y1$newObs),ncol=2),
                 rep(0.5,6),  #set initial values
                 my_dlmfc, lower=c(-1,0,0,0,-Inf,-Inf) # as dlmMLE uses optim L-BFGS-B method by default we can set lower/upper bounds to parameters to keep us out of trouble
)
my_comp <- data.frame(row.names=c("rho1","sig_z","sig_u1","sig_u2","mu_z","mu_y"),real=c(rho1,sig_z,sig_u1,sig_u2,mu_z,mu_y), mle=my_mle$par)
knitr::kable(my_comp,digits=3, caption = "Checking MLE against simulated data")
