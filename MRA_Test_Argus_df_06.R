#load neccessary libraries
#
source("C:/Users/dredmond/OneDrive - Analog Devices, Inc/Documents/CbM/r-code/MRA_libs_06.R")

####Fetch files from Zybo Spectraquest dataset
#dataset directory
# YOU will need to EDIT the "root_dir" to match your location
#

# fetch a list of files

root_dir    <- "C:/Users/dredmond/OneDrive - Analog Devices, Inc/Data/SpectraQ"
tot_dir     <- list.dirs(root_dir, full.names=TRUE,  recursive=T)
use_Arg_dir <- tot_dir[grep(tot_dir, pattern = "ArgusM")]
infiles     <- list.files(file.path(use_Arg_dir),pattern="*.csv",recursive=T, full.names=TRUE)

# Fault conditions are in one of 30 dir
# there are many files within these dir

#fault_cond   <- infiles%>%str_extract("Condition..")%>%as.factor%>%levels()
fault_Arg    <- infiles%>%str_extract(".....Fault")%>%as.factor%>%levels()

## for each infile with fault_cond[i]
# each file is a single measurement lasting 1 sec with 500k sample points
# each file has  X-Y-Z axis columns and add the label of condition
# stack the data into a long skinny set

for (i in 1:3 ){ # seq_along(infiles)) {
  fault_Arg    <- infiles[i]%>%str_extract(".....Fault")%>%as.factor%>%levels()
  db_full      <- read_csv(infiles[i], col_names = FALSE)%>%add_column(fault_Arg)
}

#compute wavelet coefficients using Discrete Wavelet Transform
#level of decomposition : log2(length of timeseries/length(filter)+1)
sample_size = nrow(db_full) #
num_lvl <- floor(log2(((sample_size-1)/19)+1))   #  samples/sec , length of filter is 20
# wont get any meaningful data if we reduce to length below 20. This should capture sub-sychronous energy
#matrix of size 570x14 to hold wavelet variances calculated from N_lvl=14 levels of signal decomposition; 570 - total no of records in the database
egy_dcoeff <- matrix(data=0, nrow = length(infiles), ncol = (num_lvl*3))%>%as.tibble()
Cond = 0
# in the loop we perform the wavelet decomposition, storing both W and V coefficients
# next step is to INTEGRATE over the band, and reduce to a single number  representing the energy with that wavelet band (octate)
# Finally we calculate the LOG of the VARIANCE of the W coeffieients.
# Literature tells us this is a good indicator when coupled with the specific wavelet levels to identify bearing faults
#
for(i in seq_along(infiles)){
  #for(i in 1:20){
  #get 3 axis vibration data from the CSV file
  fault_cond <- infiles[i]%>%str_extract(".....Fault")%>%as.factor%>%levels()
  #fault_cond <- c(fault_cond,fault_cond,fault_cond)  #  repeated fault for 3 -axis
  Cond[i]    <- infiles[i]%>%str_extract(".....Fault")%>%as.factor%>%levels()
  dt         <- read_csv(infiles[i],col_names=FALSE)
  df_wave    <- data.frame(apply(dt,2,dwt_full)[[1]])
  #idx        <- c(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3)
  egy_dcoeff[i,] <- (df_wave)
}
#low frequency coefficients in the bandwidth (0-30Hz)

db_dwt            <- data.frame(egy_dcoeff,Cond)
colnames(db_dwt)  <- c ("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11"
                        ,"W1","W2","W3","W4","W5","W6","W7","W8","W9","W10","W11",
                        "V_W1","V_W2","V_W3","V_W4","V_W5","V_W6","V_W7"
                        ,"V_W8","V_W9","V_W10","V_W11","Condition")

write_csv(db_dwt, "./CbM/db_Arg.csv")
glimpse(db_dwt)




