library(moments)
library(tidyverse)
#library(ggplot2)
#library(robustbase)
library(reshape2)
#library(tibble)
library(entropy)
#library(readr)
#library(magrittr)
library(spectral)
library(seewave)
library(viridis)
library(lubridate)
#library(stringr)
library(e1071)
library(ggthemes)
# Source all the libraries and directory path names
# Summary files first to generate TRENDS
library(rhdf5)

# one-liner functions
crest_diff <-function(xx){seewave::crest(diff(xx),f=12000,plot=FALSE)$C[1]}

rms_d<-function(xx){
  yy<-diff(xx)
  sqrt(sum(yy^2)/length(yy))}
kurt_d<-function(xx){
  yy<-diff(xx)
  kurtosis(yy)}
rms_x = function(xx){( sqrt(sum(xx^2)/length(xx) ) )}
ent_f <-function(xx){ entropy(discretize(xx,numBins=length(xx))) }
hdf_fft_fn<-function(hdf_infile){
  h5ls(hdf_infile)
  #axis1   = h5read(file=hdf_account, name="main/axis1")
  #blocki  = h5read(file=hdf_account, name="main/block0_items") ## 21 labels  
  axis0   = h5read(file=hdf_infile, name="main/_i_table" ,compoundAsDataFrame = FALSE)  #27 labels
  block0  = h5read(file=hdf_infile, name="main/table", compoundAsDataFrame = FALSE) ## 21 columns 
  dataset<- t(block0$values_block_0)
  colnames(dataset)<-axis0
  return(dataset)
}

hdf_db<-function(hdf_infile){
  axis0 <- h5read(file=hdf_infile, name = "main/axis0")
  block0 <-h5read(file=hdf_infile,name = "main/block0_values")
  dataset<-t(block0)
  colnames(dataset)<-axis0
  return(dataset)
}