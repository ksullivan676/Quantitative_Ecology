
Richness <- read.csv("~/Quantitative_Ecology/Richness.csv")
  View(Richness)

library(zoo)
library(tseries)
library(forecast)
library(xts)

#Specify periodicity of the data (# observations = 12 per year)
rich <- ts(Richness$Richness, start=c(2006, 1), end=c(2018, 10), frequency=12)

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

#Test for stationarity
#Use Dickey-Fuller (ADF), p-value < 0.05 indicates the TS is stationary
adf.test(rich)

#Detecting autocorrelation:
    #Auto correlation function (ACF) plots display correlation between a series & its lags
    #Series correlated with lags indicates seasonal trends/components
acf(rich, lag.max=45)
pacf(rich, lag.max=45)

#Fitting ARIMA Model:
#auto.arima() helps generate optimal p, d, & q
arima.rich1 <-auto.arima(rich, trace=TRUE)

#Test out the auto-generated ARIMA model values:
    #Residuals should be normally distributed
tsdisplay(residuals(arima.rich1), lag.max=45)
#Resulting ACF & PACF plots show repeating pattern at lag 5
    #--> Try out p = 18 or q = 18
arima.rich2 <-arima(rich , order=c(5,0,1), seasonal= list(order=c(1,0,1)))
tsdisplay(residuals(arima.rich2), lag.max= 30)

#Compare observed & predicted values, try to minimize AIC
AIC(arima.rich1, arima.rich2)
par(mfrow=c(1,1))
plot(rich , typ="l"); lines(fitted(arima.rich1),col="red")

#Test for independence
    #Measuring for significant difference from white noise
    #Need p-value greater than 0.05
checkresiduals(arima.rich1, lag=36)

#Forcast richness for 6 months (prediction)
plot(forecast(arima.rich1, h=6))

############### Evaluating the potential drivers of richness ################

#Create a time series object:
doc <- ts(Richness$DOC, start=c(2006, 1), end=c(2018, 10), frequency=12)

#Visualize the data:
par(mfrow=c(1,1), mai=c(0.5,0.8,0.1, 0.1))
plot(doc , typ="l", ylab= "[DOC]", xlab="")

#Remove outliers:
plot(toc , typ="l", ylab= "[DOC]", xlab="")
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

#Look for significant lags:
ccf(diff(doc),rich, na.action = na.pass, lag.max=40, plot=TRUE)

#Explore models of NEE:
#auto.arima() helps generate optimal p, d, & q
arima.rich3 <-auto.arima(rich, xreg=c(diff(doc), 0), trace=TRUE)

#Compare to current model
#AIC score indicates which model is better, lower AIC = better model
AIC(arima.rich1, arima.rich3)

#Plot ARIMA model
par(mfrow=c(1,1))
plot(rich , typ="l"); lines(fitted(arima.rich3),col="red")

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.rich4, lag=36)

#Try to make a better model
#Remove Jan 2006 datapoint for richness so it can be modeled using diff(doc) --> creates same number of observations for x & xreg
rich <--ts(Richness$Richness, start=c(2006, 2), end=c(2018, 10), frequency=12)

arima.rich4 <-arima(-rich, order=c(1,1,1), seasonal= list(order=c(0,1,1), period=12),  xreg=diff(doc), method="ML")
tsdisplay(residuals(arima.rich3), lag.max= 20)

#Compare observed & predicted values, try to minimize AIC
AIC(arima.rich1, arima.rich4)


#Reset richness to original number of observations:
#rich <--ts(Richness$Richness, start=c(2006, 1), end=c(2018, 10), frequency=12)
par(mfrow=c(1,1))

plot(-rich , typ="l", ylab="Richness"); lines(fitted(arima.rich3),col="red")

#Test for independence
#Measuring for significant difference from white noise
#Need p-value greater than 0.05
checkresiduals(arima.rich3, lag=36)

#Forcast richness (projection created based on number of datapoints in xreg, not value of horizon(h))
plot(forecast(arima.rich3, xreg=diff(doc)))

