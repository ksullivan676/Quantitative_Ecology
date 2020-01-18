
load("~/Quantitative_Ecology/NLM_Workshop.RData")

library(nlstools)

#Visualize data
par(ma1=c(1,1,0.1,0.1))
plot(harv$TIMESTAMP, harv$NEE, ylab=expression(paste("NEE(",mu, "mol m"^{-2} ~ s^{-1} ~ ")" )), xlab="")

#Fit light response curve
plot( NEE ~ PAR, data = day)

y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, data=day[which(day$MONTH == 07),],
         start=list(a1= -1 , ax= -1, r= 1),
         na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
summary(y)

#Starting values for nonlinear models
#1. Create a function of the model:
lrcModel <- function(PAR, a1, ax, r) {
  NEE <- (a1 * PAR * ax)/(a1 * PAR + ax) + r
  return(NEE)
}

#2. Create a function that calculates the initial values from the data (initial)
lrc.int <- function(mCall, LHS, data){
  x <- data$PAR
  y <-data$NEE
  
  r <- max(na.omit(y), na.rm=T) # Maximum NEE
  ax <- min(na.omit(y), na.rm=T) # Minimum NEE
  a1 <- (r + ax)/2 # Midway between r and a1
  
  # Create limits for the parameters:
  a1[a1 > 0]<- -0.1
  r[r > 50] <- ax*-1
  r[r < 0] <- 1
  
  value = list(a1, ax, r) # Must include this for the selfStart function
  names(value) <- mCall[c("a1", "ax", "r")] # Must include this for the selfStart function
  return(value)
}

#Selfstart function
SS.lrc <- selfStart(model=lrcModel,initial= lrc.int)

# 3. Find initial values:
iv <- getInitial(NEE ~ SS.lrc('PAR', "a1", "ax", "r"),
                 data = day[which(day$MONTH == 07),])
iv

y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, day[which(day$MONTH == 07),],
         start=list(a1= iv$a1 , ax= iv$ax, r= iv$r),
         na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
summary(y)

res.lrc <- nlsResiduals(y)
par(mfrow=c(2,2))
plot(res.lrc, which=1)# Residulas vs fitted values (Constant Variance)
plot(res.lrc, which=3) # Standardized residuals
plot(res.lrc, which=4) # Autocorrelation
plot(res.lrc, which=5) # Histogram (Normality)

results <- nlsBoot(y, niter=100 )
summary(results)

plot(results, type = "boxplot")