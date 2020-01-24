
load("~/Quantitative_Ecology/NLM_Workshop.RData")

library(nlstools)

#Set parameters
par(ma1=c(1,1,0.1,0.1))

#Visualize data (Net ecosystem exchange rate (NEE) Oct 1990-2016)
plot(harv$TIMESTAMP, harv$NEE, ylab=expression(paste("NEE(",mu, "mol m"^{-2} ~ s^{-1} ~ ")" )), xlab="")

#Fit light response curve (relationship between NEE & PAR)
plot( NEE ~ PAR, data = day)

#Determine nonlinear least squares estimates of the parameters of the nonlinear model
    #a1 = apparent quantum efficiency
    #ax = max ecosystem CO2 uptake rate
    #r = ecosystem respiration
y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, data=day[which(day$MONTH == 07),],
         start=list(a1= -1 , ax= -1, r= 1),
         na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
summary(y)

#Current starting values for a1, ax, and r are crappy, find better starting values
#1. Create a function of the model:
lrcModel <- function(PAR, a1, ax, r) {
  NEE <- (a1 * PAR * ax)/(a1 * PAR + ax) + r
  return(NEE)
}

#2. Create a function that calculates the initial values from the data:
lrc.int <- function(mCall, LHS, data){
  x <- data$PAR
  y <-data$NEE
  
  r <- max(na.omit(y), na.rm=T) # Maximum NEE
  ax <- min(na.omit(y), na.rm=T) # Minimum NEE
  a1 <- (r + ax)/2 # Midway between r and ax
  
  # Create limits for the parameters:
  a1[a1 > 0]<- -0.1
  r[r > 50] <- ax*-1
  r[r < 0] <- 1
  
  value = list(a1, ax, r) # Must include this for the selfStart function
  names(value) <- mCall[c("a1", "ax", "r")] # Must include this for the selfStart function
  return(value)
}

#selfStart function to calculate initial values
SS.lrc <- selfStart(model=lrcModel,initial= lrc.int)

# 3. Find initial values:
iv <- getInitial(NEE ~ SS.lrc('PAR', "a1", "ax", "r"),
                 data = day[which(day$MONTH == 07),])
iv

#Use the initial values in the model
y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, day[which(day$MONTH == 07),],
         start=list(a1= iv$a1 , ax= iv$ax, r= iv$r),
         na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
summary(y)

#Check the assumptions of the model
res.lrc <- nlsResiduals(y)
par(mfrow=c(2,2))
plot(res.lrc, which=1) # Residulas vs fitted values (Constant Variance)
plot(res.lrc, which=3) # Standardized residuals
plot(res.lrc, which=4) # Autocorrelation
plot(res.lrc, which=5) # Histogram (Normality)

#Bootstrap to estimate errors for the parameters via resampling data
results <- nlsBoot(y, niter=100 )
summary(results)

plot(results, type = "boxplot")

########################Exercise#################################

#1. Dataframe to store month parameter values & se
parms.Month <- data.frame(
  MONTH=numeric(),
  a1=numeric(),
  ax=numeric(),
  r=numeric(),
  a1.pvalue=numeric(),
  ax.pvalue=numeric(),
  r.pvalue=numeric(), stringsAsFactors=FALSE, row.names=NULL)
parms.Month[1:12, 1] <- seq(1,12,1) # Adds months to the file

#2. Function to fit the model and extract parameters
nee.day <- function(dataframe){ y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, dataframe,
                                         start=list(a1= iv$a1 , ax= iv$ax, r= iv$r),
                                         na.action=na.exclude, trace=F,
                                         control=nls.control(warnOnly=T))
y.df <- as.data.frame(cbind(t(coef(summary(y)) [1:3, 1]), t(coef(summary(y)) [1:3, 4])))
names(y.df) <-c("a1","ax", "r", "a1.pvalue", "ax.pvalue", "r.pvalue")
return (y.df )}

#3. Loop to fit monthly curves and add parameters to dataframe (parms.Month)
try(for(j in unique(day$MONTH)) {
  # Determines starting values:
  iv <- getInitial(NEE ~ SS.lrc('PAR', "a1", "ax", "r"), data = day[which(day$MONTH == j),])
  # Fits light response curve:
  y3 <- try(nee.day(day[which(day$MONTH == j),]), silent=T)
  # Extracts data and saves it in the dataframe
  try(parms.Month[c(parms.Month$MONTH == j ), 2:7 ] <- cbind(y3), silent=T)
  rm(y3)
}, silent=T)
parms.Month

#4. Bootstrapping: Create file to store parms and se
boot.NEE <- data.frame(parms.Month[, c("MONTH")]); names (boot.NEE) <- "MONTH"
boot.NEE$a1.est <- 0
boot.NEE$ax.est<- 0
boot.NEE$r.est<- 0
boot.NEE$a1.se<- 0
boot.NEE$ax.se<- 0
boot.NEE$r.se<- 0
for ( j in unique(boot.NEE$Month)){
  y1 <-day[which(day$MONTH == j),] # Subsets data
  # Determines the starting values:
  iv <- getInitial(NEE ~ SS.lrc('PAR', "a1", "ax", "r"), data = y1)
  # Fit curve:
  day.fit <- nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, data=y1,
                  start=list(a1= iv$a1 , ax= iv$ax, r= iv$r),
                  na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
  # Bootstrap and extract values:f
      #"try" ensures loop will continue to move through the months even if a model does not converge
  try(results <- nlsBoot(day.fit, niter=100 ), silent=T)
  try(a <- t(results$estiboot)[1, 1:3], silent=T)
  try(names(a) <- c('a1.est', 'ax.est', 'r.est'), silent=T)
  try( b <- t(results$estiboot)[2, 1:3], silent=T)
  try(names(b) <- c('a1.se', 'ax.se', 'r.se'), silent=T)
  try(c <- t(data.frame(c(a,b))), silent=T)
  # Add bootstrap data to dataframe:
  try(boot.NEE[c(boot.NEE$MONTH == j), 2:7] <- c[1, 1:6], silent=T)
  try(rm(day.fit, a, b, c, results, y1), silent=T)
}
lrc <- merge( parms.Month, boot.NEE, by.x="MONTH", by.y="MONTH") # Merge dataframes
lrc

#############################CHALLENGE###############################################

library(nlstools)
rm(list=ls())

load("~/Quantitative_Ecology/NLM_Workshop.RData")

#Visualize data
plot(NEE ~ TA, data = night)
plot(TA ~ MONTH, data = night)

# Selfstart for the trc:
trcModel <- function(TA, a, b) {
  y=a * exp(b*TA) #arrhenius equation
  return(y)
}

# Create a function to find initial values for the selfstart function:
trc.int <- function (mCall, LHS, data){
  x <- data$TA
  y <- data$NEE
  
  a <-1.00703982 + -0.08089044* (min(na.omit(y)))
  b <- 0.051654 + 0.001400 * (min(na.omit(y))) 
  
  value = list(a, b)
  names(value) <- mCall[c("a", "b")]
  return(value)
}

# Selfstart Function
SS.trc <- selfStart(model=trcModel,initial= trc.int)

#__________________________________________________________________________________________
#__________________________________________________________________________________________
#__________________________________________________________________________________________

# Fitting monthly models:

# Create Dataframe to store month parameter values
parms.Month <- data.frame(
  MONTH=numeric(),
  a=numeric(),
  b=numeric(), 
  a.pvalue=numeric(),
  b.pvalue=numeric(), stringsAsFactors=FALSE, row.names=NULL)


parms.Month[1:12, 1] <- seq(1,12,1) # Creates time file to merge with parm file:

#Functions: (to fit model)
nee.night <- function(dataframe){y.df = nls(NEE ~ a * exp(b*TA), 
                                            dataframe, start=list(a= iv$a , b=iv$b ),
                                            na.action=na.exclude, trace=F,
                                            control=nls.control(warnOnly=T))

y.df <- as.data.frame(cbind(t(coef(summary(y.df))[1:2, 1]), t(coef(summary(y.df)) [1:2, 4])))

names(y.df) <- c("a", "b", "a.pvalue", "b.pvalue")                      
return(y.df)}

# This loop fits monthly models (1:12):
try(for(j in unique(night$MONTH)){
  print(j)
  
  iv <- getInitial(NEE ~ SS.trc('TA', "a", "b"), data = night[which(night$MONTH == j),]) 
  
  #fit temp response curve
  y4 <- try(nee.night(night[which(night$MONTH == j),]), silent=T) # Fit night model
  
  try(parms.Month[c(parms.Month$MONTH == j ), 2:5 ] <- cbind(y4), silent=T)
  print(y4)
  rm(y4)
}, silent=T)

parms.Month

# Create file to store parms and se
boot.NEE <- data.frame(parms.Month[, c("MONTH")]); names (boot.NEE) <- "MONTH"
boot.NEE$a.est<- 0
boot.NEE$b.est<- 0
boot.NEE$a.se<- 0
boot.NEE$b.se<- 0

# Night Model & bootstrap to find how accurate the curve is
for ( j in unique(boot.NEE$MONTH)){
  print(j)
  y1 <-night[which(night$MONTH == j),]
  
  #Find intitial values
  iv <- getInitial(NEE ~ SS.trc('TA',"a", "b"), data = y1) 
  
  #Use initial values in the model
  night.fit <- nls(NEE ~ a * exp(b*TA), 
                   data=y1, start=list(a= iv$a , b=iv$b ),
                   na.action=na.exclude, trace=F,
                   control=nls.control(warnOnly=T))
  
  results <- nlsBoot(night.fit, niter=100 )
  a <- t(results$estiboot)[1, 1:2]
  names(a) <- c('a.est', 'b.est')
  b <- t(results$estiboot)[2, 1:2]
  names(b) <- c('a.se', 'b.se')
  c <- t(data.frame(c(a,b)))
  boot.NEE[c(boot.NEE$MONTH == j), 2:5] <- c[1, 1:4]
  rm(day.fit, a, b, c, results, y1)
}

trc <- merge( parms.Month, boot.NEE)

