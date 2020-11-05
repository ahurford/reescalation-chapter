# Title: Linear SIR model fitting to data
# Author: Amy Hurford
# Date: November 4, 2020
# ---------

rm(list = ls())
# load packages 
require(chron)
require(deSolve)
require(bbmle)

# set working directory
setwd('~/Desktop/HW-re-escalation/code/')

##-------- DATA CLEANING

# Pull in Canadian COVID data from Isha Berry repository
COVID.data <- read.csv('https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/timeseries_prov/active_timeseries_prov.csv', fill = TRUE)

# CLEANING THE NL DATA
NL.data <- COVID.data[COVID.data$province=="NL",]
# Begin when >= 10 cumumulative cases
NL1 <- min(which(NL.data$cumulative_cases>=10))
# Convert dates to days since a reference date
days.since <- julian(as.Date(NL.data$date_active, format = "%d-%m-%Y"),origin=as.Date(NL.data$date_active[NL1], format = "%d-%m-%Y"))
NL <- length(days.since)
start <- which(days.since==0)
# Day to end on
NL <-65
# Active cases for NL on the timeframe of interest
active.cases.NL <- NL.data$active_cases[start:(NL+start)]
# days since t=0 on timeframe of interest
days.since.NL <-days.since[start:(NL+start)]
# Final size from the data
final.size.NL <-NL.data$cumulative_cases[(NL+start)]
#plot(days.since.NL, active.cases.NL, ylab = "active cases", xlab = "days since March 23, 2020", main = "Newfoundland and Labrador", typ = "l")
# days until first less than 5 active cases from NL data
NL.end=days.since.NL[max(which(active.cases.NL>=5))]
# The maximum number of active cases from the NL data
max.NL = max(active.cases.NL)
# The days since t=0 when the peak occurs (from the NL data)
SDstart.NL <- days.since.NL[which(active.cases.NL==max(active.cases.NL))]
# Output a summary of the quantities of interst from the data
NL.data.summary = c("start date" = NL.data$date_active[NL1], "max" = max.NL, "final size" = final.size.NL, "duration" = NL.end - SDstart.NL)

## CLEANING THE NB DATA
NB.data <- COVID.data[COVID.data$province=="New Brunswick",]
# Begin the data when there are first >= 10 active cases
NB1 <- min(which(NB.data$cumulative_cases>=10))
# Days since a reference date
days.since <- julian(as.Date(NB.data$date_active, format = "%d-%m-%Y"),origin=as.Date(NB.data$date_active[NB1], format = "%d-%m-%Y"))
NB <- length(days.since)
start <- which(days.since==0)
# Day to end the data set
NB <-50
# active cases on the timeframe of interest
active.cases.NB <- NB.data$active_cases[start:(NB+start)]
# days since on the timeframe of interest
days.since.NB <-days.since[start:(NB+start)]
# final size from the NB data
final.size.NB <-NB.data$cumulative_cases[(NB+start)]
#plot(days.since.NB, active.cases.NB, ylab = "active cases", xlab = "days since March 12, 2020", main = "New Brunswick", typ="l")
# days until first less than 5 active cases from NB data
NB.end=days.since.NB[max(which(active.cases.NB>=5))]
# Maximum number of active cases in the data for NB
max.NB = max(active.cases.NB)
# The day corresponding to the maximum number of active cases
SDstart.NB <- min(days.since.NB[which(active.cases.NB==max(active.cases.NB))])
# Output a summary of the NB data
NB.data.summary = c("start date" = NB.data$date_active[NB1], "max" = max.NB, "final size" = final.size.NB, "duration" = NB.end - SDstart.NB)

#------------- FUNCTIONS
# negative log-likelihood function for the period
# prior to restrictions
nLL1 <- function(lambda,m){
  # preallocate model predictions at the times corresponding to the data
  Model.pred=NULL
  # I(t) - the model predicted number of active cases
  active.out <- (Istart+m/lambda)*exp(lambda*tvec) - m/lambda
  for(i in seq(1,length(time.data))){
    ii = min(which(tvec>=time.data[i]))
    # Need the model predictions for the times corresponding to the
    # data
    Model.pred[i] = active.out[ii]
  }
  # The negative loglikelihood. The RHS has dropped constants
  # that don't affect the location of the maximum
  # This assumes a Gaussian error distribution.
  res = sum((Model.pred - data)^2)
}

# The negative log-likelihood function, but for after when
# restrictions are enacted. Here it is assumed that there are no
# importations
nLL2 <- function(lambda){
  Model.pred=NULL
  # The equation for I(t) - the model predictions
  active.out <- Istart*exp(lambda*tvec)
  for(i in seq(1,length(time.data))){
    ii = min(which(tvec>=time.data[i]))
    # Extracting the model predictions for the time points
    # corresponding to those found in the data.
    Model.pred[i] = active.out[ii]
  }
  # The negative LL. The RHS has dropped constants
  # that don't affect the location of the maximum
  res = sum((Model.pred - data)^2)
}

# This function is I(t) - given globally defined lambda1, lambda2,
# and m1, this function plots the model predicted active number of
# infections
SIR <- function(){
  res1 <- (I0 + m1/lambda1)*exp(lambda1*tvec1) - m1/lambda1
  res2 <- (Ip+m2/lambda2)*exp(lambda2*(tvec2-SDstart)) - m2/lambda2
  res <- c(res1,res2)
}

# Function outputs all the model predicted quantities.
# lambda1: exponential growth rate prior to restrictions
# lambda2: exponential growth rate after restrictions
# m: importation rate prior to restrictions
# I0: number of active cases at t = 0
# tr: the time that restrictions are enacted
Model.pred <- function(lambda1, lambda2, m, I0,tr){
I.tr = (I0 + m/lambda1)*exp(lambda1*tr) - m/lambda1
# Assumed duration till reported recovered
gamma = 1/14
# Cumulative cases at restrictions
c.plus = gamma*((I0 + m/lambda1)*(exp(lambda1*tr)-1)-m*tr)/lambda1
# Cumulative cases after restrictions
c.minus = -(I.tr/lambda2)*gamma
# Model-predicted final size
final.size = c.plus + c.minus
# Model-predicted duration of restrictions
dur = log(5/I.tr)/lambda2
# Output all the model-predicted values
mod.summary = c("lambda1"=round(lambda1, digits =3), "lambda2" = round(lambda2, digits = 3), "m"=round(m,digits=3), "max" = round(I.tr, digits=3), "final size" = round(final.size,digits=3), "duration" = round(dur,digits=3))
return(mod.summary)
}

## ------ NL DATA FITTING
# Make time vector the the model predictions
incr <- 0.1
tvec1 = seq(0,SDstart.NL,incr)
tvec2 = seq(SDstart.NL,NL,incr)
# Date of restrictions
SDstart <- SDstart.NL
# prepare the data prior to restrictions (assumed to occur at the
# max number of active infections)
data<-active.cases.NL[which(days.since.NL <= SDstart.NL)]
time.data<-days.since.NL[which(days.since.NL <= SDstart.NL)]
# Global scope variables present in likelihood function
I0 <- active.cases.NL[1]
Istart <- I0
tvec <- tvec1
# Estimate lambda1 and m1 prior to restrictions
fit = mle2(nLL1, start=list(lambda=0.05, m=10))
lambda1 <- unname(coef(fit))[1]
m1 <- unname(coef(fit))[2]

# Prepare to use the nLL2 function - for fitting after restrictions
# i.e., after the peak in the number of active infections
time.data<-days.since.NL[which(days.since.NL >= SDstart.NL)]-SDstart.NL
data<-active.cases.NL[which(days.since.NL >= SDstart.NL)]
# Istart is cases at peak
Ip = (I0 + m1/lambda1)*exp(lambda1*SDstart.NL) - m1/lambda1
Istart <- Ip
tvec <- tvec2-SDstart.NL
# Estimate lambda2 - the exponential coefficient after restrictions
fit = mle2(nLL2, start=list(lambda=-0.1))
lambda2 <- unname(coef(fit))[1]
m2 <- 0
# Record the model predicted values
NL.pred = Model.pred(lambda1, lambda2, m1, I0,SDstart.NL)
# Print the model predictions and the analygous quantities in the data
print(NL.pred)
print(NL.data.summary)

# Make the NL part of the 2-panel figure
png(file = "Examples.png", width = 1000)
par(mfrow = c(1,2))
plot(c(tvec1,tvec2),SIR(), lwd = 3, col = "darkgrey",ylab = "active cases", xlab = "days since March 23, 2020", main = "Newfoundland and Labrador", bty="L", typ= "l", ylim = c(0,200), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
points(days.since.NL, active.cases.NL, pch=20)
lines(c(SDstart.NL,SDstart.NL+42.7), c(5,5), col = "darkgrey",lwd=3)
lines(c(SDstart.NL,NL.end), c(5,5), lwd=3)
text(25,12, labels="restrictions", col = "black", cex=1.2)

## ------ NB DATA FITTING
# Date when restrictions begin in NB - assumed to correspond to the
# maximum number of active cases
SDstart <- SDstart.NB
# Time vectors for the model predictions
tvec1 = seq(0,SDstart.NB,incr)
tvec2 = seq(SDstart.NB,NB,incr)
# Variables and data corresponding to prior to restrictions
I0 <- active.cases.NB[1]
data<-active.cases.NB[which(days.since.NB <= SDstart.NB)]
time.data<-days.since.NB[which(days.since.NB <= SDstart.NB)]
I0 <- active.cases.NB[1]
Istart <- I0
tvec <- tvec1
# Estimate the parameters prior to restrictions
fit = mle2(nLL1, start=list(lambda=0.1,m=1))
lambda1 <- unname(coef(fit))[1]
m1 <- unname(coef(fit))[2]

# Variables and data corresponding to after to restrictions
time.data<-days.since.NB[which(days.since.NB >= SDstart.NB)]-SDstart.NB
data<-active.cases.NB[which(days.since.NB >= SDstart.NB)]
Ip <- (I0+m1/lambda1)*exp(lambda1*SDstart.NB) - m1/lambda1
Istart <- Ip
tvec <- tvec2-SDstart.NB
# Estimate the parameter
fit = mle2(nLL2, start=list(lambda=-0.1))
lambda2 <- unname(coef(fit))[1]
NB.pred = Model.pred(lambda1, lambda2, m1, I0,SDstart.NB)
# Print out the model predictions and the corresponding quantities
# in the data
print(NB.pred)
print(NB.data.summary)

# Plot the NB part of the figure
plot(c(tvec1,tvec2),SIR(), lwd = 3, col = "darkgrey", xlab = "days since March 18, 2020", main = "New Brunswick", bty="L", typ = "l", ylim = c(0,85), ylab = "active cases", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
points(days.since.NB, active.cases.NB, pch = 20)
lines(c(SDstart.NB,SDstart.NB+30), c(5,5), col = "grey",lwd=3)
lines(c(SDstart.NB,NB.end), c(5,5),lwd=3)
text(25,8, labels="restrictions", col = "black", cex=1.2)
dev.off()

