# ***Script by Rory R. Koenen, January 2020, version V1.7
# ***Correction method and models by Jan Rosing, DOI: 10.1160/TH15-04-0354 and from Morrison JF, DOI: 10.1016/0968-0004(82)90157-8
# Perform substrate correction and kinetic fitting using slow, tight-binding on Xa inactivation by TFPI
library(here)
require(stats)

#*******************************************
#*** section 1, data import and cleaning ***
#*******************************************
rm(list = ls()) # clear workspace
setwd(here())
date.1<-format(Sys.time(), "%Y-%m-%d")
time.1<-format(Sys.time(), "%H-%M")
filepath<-here(date.1)
dir.create(filepath)
#Xa.df<- read.table(pipe("pbpaste"), sep="\t", header=T) #use this when running on mac
Xa.df<- read.table(file="clipboard", sep="\t", header=T) #use this when running on windows
dim.1<-dim(Xa.df) #dimensions of the data
dim.2<-dim.1[2]-1 # take dimension of data frame minus the time vector.
names(Xa.df)<-make.names(colnames(Xa.df), unique = TRUE) #make names syntactically correct 
names.old <- names(Xa.df) #store the old names
names.1<-paste("XaCurve", 1:dim.2, sep="_") #make a nice list of universal names
names(Xa.df)<- c("time", names.1) #set correct labels as time and "XaCurve_n"
Xa.df<-as.data.frame(Xa.df, row.names = NULL) # convert into a dataframe
#the following lines check for empty cells, rows, columns:
if (any(is.na(Xa.df))){
Xa.df<-Xa.df[, colSums(is.na(Xa.df)) != nrow(Xa.df)] # remove columns only containing NA
Xa.df <- Xa.df[-which(!complete.cases(Xa.df)),] # remove rows containing NA
#Xa.df<-apply(Xa.df, 2, function(x) {ifelse(is.na(x), 0, x)})# replace missing values with zeros
}
if (sum(Xa.df[1,])!=0){
  Xa.df<- rbind(seq(0,0, length.out = dim.2), Xa.df)} # add zeros to first row if row does not contain zeros
Xa.df<-apply(Xa.df, 2, function(x) {ifelse(x < 0, 0, x)})  # replace any negative values with zeros
dim.1<-dim(Xa.df)#update dimensions as columns might have been removed
dim.2<-dim.1[2]-1
names.old <- names.old[1:dim.1[2]] #shorten the names.old vector if columns were removed.
Xa.df<-as.data.frame(Xa.df)  # convert data into a dataframe ("apply" turned it into a matrix)
Xa_cor<-as.data.frame(Xa.df[,1:2])

#***************************************
#*** section 2, substrate correction ***
#***************************************

#fit initial slope
initial_slope <-seq(0,0, length.out = dim.2)
x<-Xa_cor[4:13,1]
y<-Xa_cor[4:13,2] 
initial_slope<-lm(y~x) # we only take the first part (10 points) of the curve and ignore the first 3 timepoints because they sometimes cause trouble.

blank<- initial_slope$coefficients[1] #blank
slope.1<-initial_slope$coefficients[2] #slope

Xa_cor$corr<-Xa_cor[,2]-blank #correct all values with the blank
Xa_cor[1,3]<-0 #make the "zero" zero again

x<-Xa_cor[2:11,1]
y<-Xa_cor[2:11,3] # now the blank-corrected values are y 
initial_slope<-lm(y~0+I(x)) # we force the trendline through zero now.
slope.2<-initial_slope$coefficients[1] #new slope = should not differ from old slope

Xa_cor$linfit<- (Xa_cor[,1]*slope.2)
x <- Xa_cor$time
y <- Xa_cor$corr
l <- Xa_cor$linfit
plot(x,y,main = "Xa_curve + slope_fit", pch="°", xlab = "time", ylab = "mOD", xlim = c(0,50), ylim = c(0,600))
lines(x,l, col="red")

n<- length(Xa_cor[,1]) # set length of vector

#fit using 4-degree polynomial 
l<-as.integer(n*1) #this value, which should be 1 or less is critical. It should neither be too high or too low and can be empirically determined.
x<-Xa_cor$time[1:l] #set x values
y<-Xa_cor$corr[1:l] #set y values
polyfit.1<-lm(y ~ poly(x,4, raw=TRUE))
sum.fit.1<-summary(polyfit.1)
params.1<-as.numeric(c(coefficients(polyfit.1),0,0)) #store model parameters in a vector to use for the calculation below
coef.diff<-abs(slope.2-params.1[2])
coef.diff
#plot(polyfit.1) unhash to plot model details

#predict using the parameters from the 4-degree polynomial fit (polyfit.1)
#below is the function "fourth_cal" of the 4th degree polynomial to be fitted:
fourth_cal<-function(x){
  a<-params.1[1]
  b<-params.1[2]
  c<-params.1[3]
  d<-params.1[4]
  e<-params.1[5]
  y=a+(b*x)+(c*x^2)+(d*x^3)+(e*x^4)
  return(y)
}
Xa_fit.df<-as.data.frame(cbind(Xa_cor$time[1:l],Xa_cor$corr[1:l]))
times<- Xa.df$time[1:l]
datafit<-sapply(times, fourth_cal) #the sapply command runs the polynomial function over all values
Xa_fit.df$datafit<-datafit
names(Xa_fit.df)<- c("t", "data", "fit")

#make a plot of the fitted uncorrected Xa-curve
Xa_curve<-as.data.frame(Xa_fit.df)
names(Xa_curve)<- c("t", "data", "fit")
x<-Xa_curve$t
y<-Xa_curve$data
f<-Xa_curve$fit
plot(x,y,main = "Xa_curve + 4th degree fit", pch="°", xlab = "time", ylab = "mOD", xlim = c(0,50), ylim = c(0,600))
lines(x,f, col="red")

#now calculate the residuals from the linear fit of the initial slope and plot the result
Xa_fit.df$linfit<- (Xa_cor$linfit[1:l])
Xa_fit.df$res<-Xa_fit.df$linfit-Xa_fit.df$fit

x<-Xa_fit.df$data
y<-Xa_fit.df$res
ylim.1<-1.1*max(Xa_fit.df$res)
ylim.1b<- 1.1*min(Xa_fit.df$res)
xlim.1<-1.2*max(Xa_fit.df$data)
plot(x,y,main = "residuals vs corrected mOD", pch="°", xlab = "mOD", ylab = "residuals", xlim = c(0,xlim.1), ylim = c(ylim.1b,ylim.1))
#lines(x,f, col="blue")
#as one can see, the residuals increase as mOD increases, due to substrate consumption.

#fit the curve residuals versus mOD using a 6th degree polynomial
x<-Xa_fit.df$data[1:l]
y<-Xa_fit.df$res[1:l]
polyfit.2<-lm(y ~ poly(x,6, raw=TRUE))
sum.fit.2<-summary(polyfit.2)

#predict residuals using a defined 6-degree polynomial function
params.2<-as.numeric(coefficients(polyfit.2))
sixt_cal<-function(x){
  a<-0
  b<-params.2[2]
  c<-params.2[3]
  d<-params.2[4]
  e<-params.2[5]
  f<-params.2[6]
  g<-params.2[7]
  y=a+(b*x)+(c*x^2)+(d*x^3)+(e*x^4)+(f*x^5)+(g*x^6)
  return(y)
}

#predict data using 6-degree polynomial
corr_mOD<- Xa_fit.df$data
Xa_fit_val<-sapply(corr_mOD,sixt_cal) #the sapply command runs the polynomial function over all values. Note: "predict()" returns erratic results
Xa_fit.df$val<-Xa_fit_val
Xa_fit.df$cor<-Xa_fit.df$val+Xa_fit.df$data
res.cor<-Xa_fit.df$cor-Xa_fit.df$linfit
sumsq_rescor<-sum( (res.cor - mean(res.cor) )^2 )

x<-Xa_fit.df$data
y<-Xa_fit.df$res
f<-Xa_fit.df$val
ylim.1<-0.5*max(Xa_fit.df$res)
ylim.1b<- 0.1*min(Xa_fit.df$res)
xlim.1<-0.8*max(Xa_fit.df$data)
plot(x,y,main = "residuals  + 6th degree fit", pch="°", xlab = "mOD", ylab = "residuals", xlim = c(0,xlim.1), ylim = c(ylim.1b,ylim.1))
lines(x,f, col="blue")

#plot to see if the corrected curve now lies over the ideal fit

x<-Xa_fit.df$t
y<-Xa_fit.df$cor
f<-Xa_fit.df$linfit
ylim.1<-1.1*max(Xa_fit.df$linfit)
xlim.1<-1.1*max(Xa_fit.df$t)
plot(x,y,main = "cons. corrected mOD + initial fit", pch="°", xlab = "time", ylab = "mOD", xlim = c(0,xlim.1), ylim = c(0,ylim.1))
lines(x,f, col="blue")

params<-rbind(params.1,params.2)
colnames(params)<-c("intercept","X1","X2","X3","X4","X5","X6")
filename<-paste(filepath,"/","5_model_parameters",time.1,".csv", sep = "")
write.csv(params, file = filename) #saving a csv with the fitted parameters
rm(Xa_cor) #clean up obsolete dataframes
rm(Xa_curve) #clean up obsolete dataframes

#************************************************
#*** section 3, correcting experimental data  ***
#************************************************

#import the data to be corrected for substrate consumption from a new Excel file
TFPI_raw.df<-Xa.df

#fit all data in the imported dataset
TFPI_cor.df<-as.data.frame(TFPI_raw.df[,1])
for (i in 2:dim.1[2]){
corr_mOD<- TFPI_raw.df[,i]-blank
corr_mOD<-sapply(corr_mOD, function(x) {ifelse(x < 0, 0, x)}) #make any negatives zero
TFPI_cor.1<-ifelse(corr_mOD>blank,sapply(corr_mOD,sixt_cal),0) #only correct values that are over the blank
TFPI_cor.df[,i]<-corr_mOD+TFPI_cor.1
}

names.3<-paste("TFPI_cor", 1:dim.2, sep="_") #make a nice list of universal names
names(TFPI_cor.df)<- c("time", names.3) #set new labels as time and "TFPI_fit_n"

for(i in 2:dim.1[2]){
  x<-TFPI_raw.df[,1]
  y<-TFPI_raw.df[,i]
  xc<-TFPI_cor.df[,1]
  yc<-TFPI_cor.df[,i]
  ylim.1<-1.2*max(TFPI_cor.df[,i])
  plot(x,y,main = paste("corrected vs raw curves","\n",names.old[i]), xlab = "time", ylab = "mOD", ylim = c(0,ylim.1))
  par(new=TRUE) 
  lines(xc,yc,col="blue")
}

names<-gsub("[X]","",names.old)#remove the "X" before the column titles, may yield unwanted resuls
names(TFPI_cor.df)<-names
filename.1<-paste(filepath,"/","1_corrected TFPI curves",time.1,".csv", sep = "")
write.csv(TFPI_cor.df, file=filename.1, row.names = FALSE)
rm(TFPI_raw.df) #clean up obsolete dataframes

#*********************************************
#*** section 4, fitting experimental data  ***
#*********************************************

#fit the corrected data with the kinetic equation

#first define the function
#initial parameters you can adjust:
A0<-0
kobs<-0.1
V0<-10.3
Vs<-0

params.3<-list(A0=A0,kobs=kobs,V0=V0,Vs=Vs) #place the parameters in a list

TFPI_kinfit<-function(t,A0,kobs,V0,Vs){
  A=A0+(Vs*t)+(V0-Vs)*((1-exp(-kobs*t))/kobs)
  return(A)
} 
TFPI_kinfit(2.08,params.3[[1]],params.3[[2]],params.3[[3]],params.3[[4]]) #I found the function to be correct

#fit the corrected data

dim.1<-dim(TFPI_cor.df) #dimensions of the data
dim.2<-dim.1[2]-1 # take dimension of data frame minus the time vector.
names.2<-paste("TFPI_raw", 1:dim.2, sep="_") #make a nice list of universal names
names(TFPI_cor.df)<- c("time", names.2) #set correct labels as time and "TFPI_curve_n"
test_raw.df<-as.data.frame(TFPI_cor.df, row.names = NULL) # convert into a dataframe
test_raw.df[1,]<-0

#first fit initial slopes (V0) and calculate correlation coefficients
initial_slope.V0 <-seq(0,0, length.out = dim.2)
intercept.V0<-seq(0,0, length.out = dim.2)
corcoefs<-NULL
for(i in 2:dim.1[2]){
x<-test_raw.df[2:8,1]
y<-test_raw.df[2:8,i] 
linearfit<-lm(y~x)
coefs.1<-coefficients(linearfit)
initial_slope.V0[i]<-coefs.1[2]
intercept.V0[i]<-coefs.1[1]
corcoefs[i]<-(cor(x,y, method = "pearson"))^2
}

total_slope <-seq(0,0, length.out = dim.2)
corcoefs.total<-NULL
for(i in 2:dim.1[2]){
  x<-test_raw.df[,1]
  y<-test_raw.df[,i] 
  linearfit<-lm(y~x)
  coefs.1<-coefficients(linearfit)
  total_slope[i]<-coefs.1[2]
  corcoefs.total[i]<-(cor(x,y, method = "pearson"))^2
}

#now fit the curves using this TFPI kinfit function, with a logical test for linearity
#linear curves return infinity and are not fitted, but only return slope (calculated above)
#unfittable curves return an error and all parameters are 0
test_fit.df<-as.data.frame(test_raw.df[,1])
kinparams.df<-as.data.frame(params.3)
for (i in 2:dim.1[2]){
  y<-test_raw.df[,i]
  t<-test_fit.df[,1]
  if (corcoefs.total[i]<0.99){
  TFPI.fit.1<-tryCatch({nls(y ~ TFPI_kinfit(t,A0,kobs,V0,Vs), start = list(kobs=params.3[[2]],V0=params.3[[3]],Vs=params.3[[4]]))},error=function(e){cat("ERROR :","curve",i,conditionMessage(e), "\n")})
  kinparams_int<-c(0,coefficients(TFPI.fit.1))
  kinparams.df[i,]<-kinparams_int
  test_fit.df[,i]<-TFPI_kinfit(t,0,kinparams.df[i,2],kinparams.df[i,3],kinparams.df[i,4])
  } else {
    TFPI.fit.1<- t*initial_slope.V0[i]+intercept.V0[i]
    kinparams_int<-c(0,0, total_slope[i],0)
    kinparams.df[i,]<-kinparams_int
    test_fit.df[,i]<-TFPI.fit.1
  }
}

#plot the curve with kinetic fit
fit_res<-NULL
fit_res.df<-as.data.frame(test_raw.df[,1])
sumsq_fitres.df<-as.data.frame(max(test_raw.df[,1]))
for (i in 2:dim.1[2]){
t<-test_raw.df[,1]
A<-test_fit.df[,i]
ylim.2<-1.2*max(test_raw.df[,i])
plot(test_raw.df[,1],test_raw.df[,i], main = paste("corrected vs fitted curves","\n",names.old[i]), xlab = "time", ylab = "mOD", ylim = c(0,ylim.2)) 
par(new=TRUE) 
lines(t,A, col="red")
fit_res<-(test_raw.df[,i]-A)
fit_res.df[,i]<-fit_res
sumsq_fitres.df[,i]<-sum( (fit_res - mean(fit_res) )^2 )
}
colnames(fit_res.df)<-names.old
colnames(sumsq_fitres.df)<-names.old

#define a simple numeric 1st derivative function
first.deriv <- function(x,y){
  d<-(diff(y)/diff(x))
  return(d)
}

#take first derivative of fitted inactivation curves, normalize to 100% and plot 
#first test for "double zeros" in time, which makes first.deriv() return Inf. If "TRUE" rows are removed.
if(min(diff(test_fit.df[,1]))==0){
  test_fit.df <- test_fit.df[-which(diff(test_fit.df[,1])==0),]
}

TFPI_fitderiv.df<-as.data.frame(test_fit.df[-1,1])
for (i in 2:dim.1[2]){
  y<-test_fit.df[,i]
  x<-test_fit.df[,1]
TFPI.fitderiv.1<- first.deriv(x,y)
TFPI.fitderiv.1 <-100*(TFPI.fitderiv.1/TFPI.fitderiv.1[1])
TFPI_fitderiv.df[,i]<-TFPI.fitderiv.1
color.1<-as.hexmode(5*i)
plot(TFPI_fitderiv.df[,1],TFPI_fitderiv.df[,i], main = "Xa inactivation",type = "l", col=color.1 ,xlab = "time", ylab = "%Xa",ylim = c(0,100))
legend(40, (120-(13*(i-1))), legend= names.old[i], col=color.1, lty=1:2, cex=0.8, box.lty = 0)
par(new=TRUE)
}

#save all data in time-labelled files.
names<-gsub("[X]","",names.old)#remove the "X" before the column titles, may yield unwanted resuls
rownames(kinparams.df)<-names
kinparams.df<-kinparams.df[-1,]
names(test_fit.df)<-names
names(TFPI_fitderiv.df)<-names
filename.2<-paste(filepath,"/","4_kinetic_parameters_",time.1,".csv",sep = "")
filename.3<-paste(filepath,"/","2_fitted_curves_",time.1,".csv",sep = "")
filename.4<-paste(filepath,"/","3_Xa_inactivation_curves_",time.1,".csv",sep = "")
write.csv(kinparams.df, file=filename.2, row.names = TRUE)
write.csv(test_fit.df, file=filename.3, row.names = FALSE)
write.csv(TFPI_fitderiv.df, file=filename.4, row.names = FALSE)

#unhash to clean up workspace except dataframes after analysis
#rm(list = grep("df", ls(), value = TRUE, invert = TRUE))# clear all workspace except the dataframes
#rm(test_raw.df)





