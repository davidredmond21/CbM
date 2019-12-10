#load neccessary libraries
account= "merck"
source("./CbM/r-code/MRA_libs.R")
source("./CbM/r-code/MRA_dirs.R")
source("./CbM/r-code/MRA_ETL.R")
####Fetch files from Zybo Spectraquest dataset

#dataset directory
# fetch a list of files

root_dir="C:/Users/dredmond/Desktop/DataSet/SpectraQuestData"
tot_dir =list.dirs(root_dir,full.names=TRUE,  recursive=T)
## pick out Zybo  Only
use_zybo_dir<-tot_dir[grep(tot_dir, pattern = "*_Zybo")]
#
infiles =list.files(file.path(use_zybo_dir),pattern=".csv", full.names=TRUE)
# there are 3 fault conditions
# and 5 RPM speeds
#
cond = c("OuterRaceFault", "NoFault", "innerRaceFault")
for (i in seq_along(infiles)) {
  fault_rpm    <- infiles[i]%>%str_extract("RPM.")%>%as.factor%>%levels()
  fault_cond   <- infiles[i]%>%str_extract(cond)%>%as.factor%>%levels()
  db_full      <- read_csv(infiles[i])%>%add_column(fault_rpm,fault_cond)
}


# compute wavelet coefficients using Discrete Wavelet Transform
# level of decomposition : log2(length of timeseries/length(filter)+1)
sample_size = nrow(db_full) # 500000
num_lvl <- floor(log2(((sample_size-1)/19)+1))   # 500000 samples/sec , length of filter is 20
# wont get any meaningful data if we reduce to length below 20. This should capture sub-sychronous energy
#matrix of size 570x14 to hold wavelet variances calculated from N_lvl=14 levels of signal decomposition; 570 - total no of records in the database

column_names <- c("Bin1","Bin2","Bin3","Bin4","Bin5","Bin6","Bin7","Bin8","Bin9","Bin10","Bin11","Bin12","Bin13","Bin14")#,"Condition", "RPM")
egy_dcoeff <- matrix(nrow = length(infiles)*3, ncol = (num_lvl))%>%as.data.frame%>%setNames(.,column_names)

for(i in seq_along(infiles))
  {
#  for(i in 1:4){
  #get 3 axis vibration data from the CSV file
  idx        <- c(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3)
  fault_cond[idx] <- infiles[i]%>%str_extract(cond)%>%as.factor%>%levels()%>%c(.,.,.)
  fault_rpm[idx]  <- infiles[i]%>%str_extract("RPM.")%>%as.factor%>%levels()%>%c(.,.,.)
  dt         <- read_csv(infiles[i],col_names=FALSE)
  df_wave    <- apply(dt[,1:3],2,dwt_full)%>%t%>%as.tibble
  egy_dcoeff[idx,] <- df_wave
}
# little tidy up needed
df_egy <- cbind(egy_dcoeff,fault_cond,fault_rpm)
#
glimpse(df_egy)

write_csv(df,"./df_wavelet.csv")
########## Start the ML process

# helper function useful for ML
one_hot <- function(df, key) {
  key_col <- dplyr::select_var(names(df), !! rlang::enquo(key))
  df <- df %>% mutate(.value = 1, .id = seq(n()))
  df <- df %>% tidyr::spread_(key_col, ".value", fill = 0, sep = "_") %>% select(-.id)
}


df_egy <- read_csv("./df_egy.csv")
df_egy$fault_cond%<>%as.factor()
df_egy$fault_rpm%<>%as.factor()

#
#
library(gbm)

#GBM: Generalized Boosted Models
#an ensemble of classification or regression trees
#with the gbm package, can do both a few different kinds including
#AdaBoost and Gradient Boosting

#advantages of GBM over logistic regression (parametric method):
#robust to outliers
#can still make predictions when an observation has missing data!
#handles unequal class sizes and unbalanced predictor variables well
#(logistic regression isn't as good at this)
#you don't need to specify interaction terms with tree models!
#usually have greater predictive ability

#potential drawbacks:
#trees can overfit, especially if the number of
#ending nodes is too small or the number of trees is too large
#definitely want to use CV, can't use in-sample
#prediction rate as a measure of goodness of fit!





## STEP 1 drop the RPM classifier column - so we are ignoring the SPEED
df_egy %<>%select(-c(fault_rpm))
# One-hot encode fault
df <- one_hot(df_egy,fault_cond)
glimpse(df)
df <- as.data.frame(df)
# prep a TRAINING and VALIDATION set
set.seed(121)


idx <- seq(1:nrow(df))
train<-sample(idx,324)  # use 80% training and 20% testing
df_train <-df[train,]
df_test  <-df[-train,] # use this as the test data set
ntrees =2000 # sufficient large
#
#### Lets set up the GBM model hyparameters
# Can run this on the FULL data - it will do its own CV


full_case <- colnames(df[,-c(1:14)])
X_var     <- colnames(df[,1:14])
# in this trial we will have 3 conditions.
# Step #1 is to see how each modelling exercise performs
for (i in full_case) {}


  print (paste0("df$",i))
  response <- df_train[,i]
  y_test   <- df_test[,i]




###
library(ROCR)
library(mlr)
library(randomForest)
library(nnet)
#install.packages('mboost')
#install.packages('party')
#install.packages('tgp')
#install.packages('RRF')
#install.packages('xgboost')
library(mboost)
library(party)
library(RRF)
library(tgp)
library(xgboost)

df_rf   <-makeLearner('regr.randomForest')
df_task <- makeRegrTask(id="df",data=df,target="fault_cond_innerRaceFault")
summary(df_task)






####Classification
class_lrn = makeLearner("classif.randomForest",
                        predict.type = "prob", # Get predicted probabilities
                        fix.factors.prediction = TRUE)



# Set up 10-fold CV
desc = makeResampleDesc("CV", iters = 10)
# See also crossval which does the same thing (with less flexibility)
## MLR doesnt like tibbles
df_train$fault_cond_innerRaceFault%<>%as.factor()
df_task <- makeClassifTask(id = "fault", data = df_train[,1:15], target = "fault_cond_innerRaceFault")

class_rf = resample(class_lrn,
                    task = df_task,
                    resampling = desc,
                    measures = list(mmce))

df_gam = makeLearner('regr.gamboost')
df_tgp = makeLearner('regr.cforest')
df_RRF = makeLearner('regr.RRF')
df_xgb = makeLearner('regr.xgboost')
df_nn = makeLearner('regr.nnet')
df_svm = makeLearner('regr.svm')

