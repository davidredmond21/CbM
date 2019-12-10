data_str = "C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data"
raw_hdf_dir="C:/Users/dredmond/Documents/CbM/Python/cbmsignal/data/raw_adc_data/"
#
csv_files <- list.files(file.path(paste0(data_str,"/",account)),pattern="*csv", full.names=TRUE)
#
plots_dir <- paste0(data_str,"/",account,"/plots")


if (!file.exists(plots_dir)){
  dir.create(file.path(plots_dir))}

#
hdf_raw_files <- list.files(file.path(raw_hdf_dir),pattern=paste0(account,"_........._adc_....."), full.names=TRUE)
#
#