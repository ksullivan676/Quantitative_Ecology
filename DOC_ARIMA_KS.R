
Richness <- read.csv("~/Quantitative_Ecology/Richness.csv")
View(Richness)

library(zoo)
library(tseries)
library(forecast)
library(xts)

#Specify periodicity of the data (# observations = 12 per year)
doc <- ts(Richness$DOC, start=c(2006, 1), end=c(2018, 10), frequency=12)

#Plot the timeseries, examine for outliers/irregularity
par(mfrow=c(1,1), mai=c(1.25,0.8,0.1, 0.1))
plot(doc, typ="l", ylab= "[DOC]", xlab="")

#Remove outliers w/ tsclean()
plot(doc)
lines(tsclean(doc), col="blue")
doc <- tsclean(doc)

#Decomposition of time series (observed, trend, seasonal, & random)
doc.d <- decompose(doc, 'multiplicative')
plot(doc.d)

#Test for stationarity, ARIMA model requires series to be stationary
#Use Dickey-Fuller (ADF) to test for stationarity
#p-value < 0.05 indicates the TS is stationary
adf.test(doc)
#If null hypothesis is not rejected, try differencing the time series w/ diff()
adf.test(diff(doc))

#Detecting autocorrelation:
#Auto correlation function (ACF) plots display correlation between a series & its lags
#Series correlated with lags indicates seasonal trends/components
acf(diff(doc), lag.max=45)
pacf(diff(doc), lag.max=45)

#Fitting ARIMA Model:
#auto.arima() helps generate optimal p, d, & q
arima.doc1 <-auto.arima(diff(doc), trace=TRUE)

#Test out the auto-generated ARIMA model values:
#Residuals should be normally distributed
tsdisplay(residuals(arima.doc1), lag.max=45)

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.doc1, lag=36)

par(mfrow=c(1,1))
plot(diff(doc), typ="l"); lines(fitted(arima.doc1),col="red")

#Forcast NEE for 30 days (prediction)
plot(forecast(arima.doc1, h=6))

############### Evaluating the potential drivers of DOC ################

#Create a time series object:
wtable <- ts(Richness$wtr_table, start=c(2006, 1), end=c(2018, 10), frequency=12)

#Visualize the data:
par(mfrow=c(1,1), mai=c(0.5,0.8,0.1, 0.1))
plot(wtable , typ="l", ylab= "m AMSL", xlab="")

#Remove outliers:
plot(wtable , typ="l", ylab= "m AMSL", xlab="")
lines(tsclean(wtable) , col="blue")
wtable <- tsclean(wtable)

#Decompose the time series:
wtable.d <- decompose(wtable, 'multiplicative')
plot(wtable.d)

#Test for stationarity:
#P-value < 0.05 indicates the time series is stationary
adf.test(wtable)
#If null hypothesis is not rejected, try differencing the time series w/ diff()
adf.test(diff(wtable))

#Look for significant lags:
ccf(diff(wtable),doc, na.action = na.pass, lag.max=40, plot=TRUE)

#Explore models of NEE:
#auto.arima() helps generate optimal p, d, & q
arima.doc2 <-auto.arima(diff(doc), xreg=diff(wtable), trace=TRUE)

#Compare to current model
#AIC score indicates which model is better, lower AIC = better model
AIC(arima.doc1, arima.doc2)

#Test out the auto-generated ARIMA model values:
#Residuals should be normally distributed
tsdisplay(residuals(arima.doc2), lag.max=45)

#Try to make a better model
arima.doc3 <-arima(diff(doc), order=c(1,0,2), seasonal= list(order=c(0,0,2), period=12),  xreg=diff(wtable), method="ML")
tsdisplay(residuals(arima.rich3), lag.max= 20)

#Compare to current model
#AIC score indicates which model is better, lower AIC = better model
AIC(arima.doc2, arima.doc3)

#Test out the auto-generated ARIMA model values:
#Residuals should be normally distributed
tsdisplay(residuals(arima.doc3), lag.max=45)

#arima.doc2 had better AIC than arima.doc3, use this for model

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.doc2, lag=36)

par(mfrow=c(1,1))
plot(diff(doc), typ="l"); lines(fitted(arima.doc2),col="red")

#Forcast DOC for 6 months days (prediction)
plot(forecast(arima.doc2, xreg=wtable))


