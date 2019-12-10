library(EBImage)
library(keras)

#read images
setwd("C:/Users/dredmond/Downloads/data")
infiles.pic   <- list.files(pattern="*.jpg",recursive=T, full.names=TRUE)
#
my.pic <- list()
for (i in 1:length(infiles.pic)){
  my.pic[[i]] <- readImage(infiles.pic[i])
}
display(my.pic[[i]])
#
summary(my.pic[[1]])
str(my.pic)
# resize to same dimensions
for( j in 1:length(infiles.pic)){
  my.pic[[j]]<- resize(my.pic[[j]], 28,28)
}
# str(my.pic)
# prepare data for training
# re-shape to a 1-dim
# display(my.pic[[2]])
for ( i in 1:length(infiles.pic)){
  my.pic[[i]] <- keras::array_reshape(my.pic[[i]],c(28,28,3))
}
str(my.pic)
# everything should be uniform
trainx <- NULL
for ( i in 1:length(infiles.pic)) {
  trainx<- rbind(trainx,my.pic[[i]])
}


