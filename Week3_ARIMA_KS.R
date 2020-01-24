
load("~/Quantitative_Ecology/ARIMA_Workshop.RData")

library(zoo)
library(tseries)
library(forecast)
library(xts)

#Specify periodicity of the data (# observations = 30 per month)
nee <- ts(mangroves$nee, start=1, frequency=30)

#Plot the timeseries, examine for outliers/irregularity
par(mfrow=c(1,1), mai=c(1.25,0.8,0.1, 0.1))
plot(nee, typ="l", ylab= "NEE", xlab="")

#Remove outliers w/ tsclean()
plot(nee)
lines(tsclean(nee), col="red")
nee <- tsclean(nee)

#Decomposition of time series (observed, trend, seasonal, & random)
nee.d <- decompose(nee, 'multiplicative')
plot(nee.d)

#Test for stationarity, ARIMA model requires series to be stationary
    #Stationary = mean, variance, & autocovariance are time invariant
#Use Dickey-Fuller (ADF) to test for stationarity
    #p-value < 0.05 indicates the TS is stationary
adf.test(nee )

#Detecting autocorrelation:
    #Auto correlation function (ACF) plots display correlation between a series & its lags
    #Series correlated with lags indicates seasonal trends/components
acf(nee, lag.max=45)
pacf(nee, lag.max=45)

#Fitting ARIMA Model:
#auto.arima() helps generate optimal p, d, & q
arima.nee1 <-auto.arima(nee, trace=TRUE)

#Test out the auto-generated ARIMA model values:
    #Residuals should be normally distributed
tsdisplay(residuals(arima.nee1), lag.max=45)
#Resulting ACF & PACF plots show repeating pattern at lag 10
    #--> Try out p = 10 or q = 10
arima.nee2 <-arima(nee , order=c(10,1,3), seasonal= list(order=c(2,0,2)))
tsdisplay(residuals(arima.nee2), lag.max= 30)

#Compare observed & predicted values, try to minimize AIC
AIC(arima.nee1, arima.nee2)
par(mfrow=c(1,1))
plot(nee , typ="l"); lines(fitted(arima.nee2),col="red")

#Test for independence
    #Measuring for significant difference from white noise
    #Need p-value greater than 0.05
checkresiduals(arima.nee2, lag=36)

par(mfrow=c(1,1))
plot(nee , typ="l"); lines(fitted(arima.nee2),col="red")

#Forcast NEE for 30 days (prediction)
plot(forecast(arima.nee2, h=30))

############### Evaluating the potential drivers of NEE ################

#Create a time series object:
sal <- ts(mangroves$salinity.max, start= 1, frequency=30)

#Visualize the data:
par(mfrow=c(1,1), mai=c(0.25,0.8,0.1, 0.1))
plot(sal , typ="l", ylab= "Salinity", xlab="")

#Remove outliers:
plot(sal , typ="l", ylab= "Salinity", xlab="")
lines(tsclean(sal) , col="red")
sal <- tsclean(sal)

#Decompose the time series:
sal.d <- decompose(sal, 'multiplicative')
plot(sal.d)

#Test for stationarity:
    #P-value < 0.05 indicates the time series is stationary
adf.test(sal)
#If null hypothesis is not rejected, try differencing the time series w/ diff()
adf.test(diff(sal))

#Look for significant lags:
ccf( diff(sal),nee, na.action = na.pass, lag.max=40, plot=TRUE)

#Explore models of NEE:
#auto.arima() helps generate optimal p, d, & q
arima.nee3 <-auto.arima(nee, xreg=c(diff(sal),0), trace=TRUE)

#Compare to current model
#AIC score indicates which model is better, lower AIC = better model
AIC(arima.nee2, arima.nee3 )

#Adding salinity to NEE did not improve the model
    #--> Check if extreme salinity is more important
    #Create salinity index
sal.i <- sal
sal.i[sal.i < 25 ]<- 0
sal.i[sal.i >= 25 ]<- 1

plot(sal.i)

#Try adding extreme salinity indicator into the model & see if it improves
arima.nee4 <-auto.arima(nee, xreg=sal.i, trace=TRUE)

#Compare to current model
#AIC score indicates which model is better, higher AIC = better
AIC(arima.nee2,arima.nee4 )

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.nee4, lag=36)

par(mfrow=c(1,1))
plot(nee , typ="l"); lines(fitted(arima.nee4),col="red")

################# CAN YOU MAKE A BETTER MODEL? #################

#Try using PAR:


#Create timeseries object:
parlight<-ts(mangroves$par, start=1, frequency=30)

#Visualize data:
par(mfrow=c(1,1), mai=c(1.25,0.8,0.1, 0.1))
plot(parlight, typ="l", ylab= "PAR", xlab="")

#Remove outliers
plot(parlight , typ="l", ylab= "PAR", xlab="")
lines(tsclean(parlight) , col="red")

parlight <- tsclean(parlight)

#Decompose time series
parlight.d <- decompose(parlight, 'multiplicative')
plot(parlight.d)

#Test for stationarity:
#P-value < 0.05 indicates the time series is stationary
adf.test(parlight)
#If null hypothesis is not rejected, try differencing the time series w/ diff()
adf.test(diff(parlight))

#Look for significant lags:
ccf( diff(parlight),nee, na.action = na.pass, lag.max=40, plot=TRUE)

#Explore models of NEE:
#auto.arima() helps generate optimal p, d, & q
arima.nee5 <-auto.arima(nee, xreg=c(diff(parlight),0), trace=TRUE)

#Compare to current model
#AIC score indicates which model is better, lower AIC = better model
AIC(arima.nee4, arima.nee5 )


#________________________________________________________________
#Adding salinity to NEE did not improve the model
#--> Check if extreme salinity is more important
#Create salinity index
sal.i <- sal
sal.i[sal.i < 25 ]<- 0
sal.i[sal.i >= 25 ]<- 1

plot(sal.i)

#Try adding extreme salinity indicator into the model & see if it improves
arima.nee4 <-auto.arima(nee, xreg=sal.i, trace=TRUE)

#Compare to current model
#AIC score indicates which model is better, higher AIC = better
AIC(arima.nee2,arima.nee4 )

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.nee4, lag=36)

par(mfrow=c(1,1))
plot(nee , typ="l"); lines(fitted(arima.nee4),col="red")








