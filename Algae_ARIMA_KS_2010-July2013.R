
Richness <- read.csv("~/Quantitative_Ecology/Richness.csv")
View(Richness)

library(zoo)
library(tseries)
library(forecast)
library(xts)

#Specify periodicity of the data (1 per month)
rich_data <- Richness[c(47:91), c(2)]
rich <- ts(rich_data, start=c(2009, 11), end=c(2013, 7), frequency =12)

#Plot the timeseries, examine for outliers/irregularity
par(mfrow=c(1,1), mai=c(1.25,0.8,0.1, 0.1))
plot(rich, typ="l", ylab= "Richness", xlab="")

#Remove outliers w/ tsclean()
plot(rich)
lines(tsclean(rich), col="blue")
rich <- tsclean(rich)

#Decomposition of time series (observed, trend, seasonal, & random)
rich.d <- decompose(rich, 'multiplicative')
plot(rich.d)

#Test for stationarity, ARIMA model requires series to be stationary
#Use Dickey-Fuller (ADF) to test for stationarity
#p-value < 0.05 indicates the TS is stationary
adf.test(rich)
#If null hypothesis is not rejected, try differencing the time series w/ diff()
adf.test(diff(rich))

#Detecting autocorrelation:
#Auto correlation function (ACF) plots display correlation between a series & its lags
#Series correlated with lags indicates seasonal trends/components
acf(diff(rich), lag.max=45)
pacf(diff(rich), lag.max=45)

#Fitting ARIMA Model:
#auto.arima() helps generate optimal p, d, & q
arima.rich1 <-auto.arima(diff(rich), trace=TRUE)

#Test out the auto-generated ARIMA model values:
#Residuals should be normally distributed
tsdisplay(residuals(arima.rich1), lag.max=45)

#Compare observed & predicted values, try to minimize AIC
#AIC(arima.rich1, arima.rich2)
par(mfrow=c(1,1))
plot(diff(rich), typ="l"); lines(fitted(arima.rich1),col="red")

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.rich1, lag=36)

#Forcast richness for 6 months (prediction)
plot(forecast(arima.rich1, h=6))

############### Evaluating the potential drivers of Richness ################

#Create a time series object:
doc_data <- Richness[c(47:91), c(4)]
doc <- ts(doc_data, start=c(2009, 11), end=c(2013, 7), frequency =12)

#Visualize the data:
par(mfrow=c(1,1), mai=c(0.5,0.8,0.1, 0.1))
plot(doc , typ="l", ylab= "[DOC]", xlab="")

#Remove outliers:
plot(doc , typ="l", ylab= "[DOC]", xlab="")
lines(tsclean(doc) , col="blue")
doc <- tsclean(doc)

#Decompose the time series:
doc.d <- decompose(doc, 'multiplicative')
plot(doc.d)

#Test for stationarity:
#P-value < 0.05 indicates the time series is stationary
adf.test(doc)
#If null hypothesis is not rejected, try differencing the time series w/ diff()
adf.test(diff(doc))
#Still not stationary --> take second difference
adf.test(diff(diff(doc)))

#Look for significant lags:
ccf(diff(diff(doc)),rich, na.action = na.pass, lag.max=40, plot=TRUE)

#Remove Nov 2009 datapoint for richness so it can be modeled using diff(diff(doc)) --> creates same number of observations for x & xreg
rich_data <- Richness[c(48:91), c(2)]
rich <- ts(rich_data, start=c(2009, 12), end=c(2013, 7), frequency =12)

#auto.arima() helps generate optimal p, d, & q
arima.rich2 <-auto.arima(rich, xreg=c(diff(diff(doc)), 0), trace=TRUE)

#Compare to current model
#AIC score indicates which model is better, lower AIC = better model
AIC(arima.rich1, arima.rich2 )

tsdisplay(residuals(arima.rich2), lag.max=45)

par(mfrow=c(1,1))
plot(rich , typ="l", ylab="Richness"); lines(fitted(arima.rich2),col="red")

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.rich2, lag=36)

#Forcast richness for 6 months (prediction) --> does not project correct amount of time when using xreg
plot(forecast(arima.rich2, xreg=diff(diff(doc))))

