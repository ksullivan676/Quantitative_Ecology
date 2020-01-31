#Species distribution model for saguaro

library("sp")
library("raster")
library("maptools")
library("rgdal")
library("dismo")

bioclim.data <- getData(name = "worldclim",
                       var = "bio",
                       res = 2.5,
                       path = "~/Quantitative_Ecology/")

#Read in saguaro observations
obs.data <- read.csv(file = "~/Quantitative_Ecology/SDM_Data.csv")
obs.data <- na.omit(obs.data)
summary(obs.data)

#Determine geographic extent of the data
max.lat <- ceiling(max(obs.data$latitude))
min.lat <- floor(min(obs.data$latitude))
max.lon <- ceiling(max(obs.data$longitude))
min.lon <- floor(min(obs.data$longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

#Load the data to use for our base map
data(wrld_simpl)

#Plot the base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

#Add the points for individual observation
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)
#And draw a little box around the graph
box()

#Crop bioclim data to geographic extent of saguaro
bioclim.data <- crop(x = bioclim.data, y = geographic.extent)

#Build species distribution model
bc.model <- bioclim(x = bioclim.data, p = obs.data)

#Drop unused column
obs.data <- obs.data[, c("latitude", "longitude")]

#Build species distribution model
bc.model <- bioclim(x = bioclim.data, p = obs.data)

#Reverse order of columns
obs.data <- obs.data[, c("longitude", "latitude")]

#Build species distribution model
bc.model <- bioclim(x = bioclim.data, p = obs.data)

#Predict presence from model
predict.presence <- dismo::predict(object = bc.model, x = bioclim.data, ext = geographic.extent)

#Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

#Add model probabilities
plot(predict.presence, add = TRUE)

#Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")

#Add original observations
points(obs.data$longitude, obs.data$latitude, col = "olivedrab", pch = 20, cex = 0.75)
box()

#Use the bioclim data files for sampling resolution
bil.files <- list.files(path = "~/Quantitative_Ecology/wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)

#We only need one file, so use the first one in the list of .bil files
mask <- raster(bil.files[1])

#Randomly sample points (same number as our observed points)
background <- randomPoints(mask = mask,     #Provides resolution of sampling points
                           n = nrow(obs.data),      #Number of random points
                           ext = geographic.extent, #Spatially restricts sampling
                           extf = 1.25)             #Expands sampling a little bit
head(background)

#Plot the base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95",
     main = "Presence and pseudo-absence points")

#Add the background points
points(background, col = "grey30", pch = 1, cex = 0.75)

#Add the observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)

box()

#Arbitrarily assign group 1 as the testing data group
testing.group <- 1

#Create vector of group memberships
group.presence <- kfold(x = obs.data, k = 5) #kfold is in dismo package

head(group.presence)

table(group.presence)

#Separate observations into training and testing groups
presence.train <- obs.data[group.presence != testing.group, ]
presence.test <- obs.data[group.presence == testing.group, ]

#Repeat the process for pseudo-absence points
group.background <- kfold(x = background, k = 5)
background.train <- background[group.background != testing.group, ]
background.test <- background[group.background == testing.group, ]

#Build a model using training data
bc.model <- bioclim(x = bioclim.data, p = presence.train)

#Predict presence from model (same as previously, but with the update model)
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data, 
                                   ext = geographic.extent)
#Use testing data for model evaluation
bc.eval <- evaluate(p = presence.test,   #The presence testing data
                    a = background.test, #The absence testing data
                    model = bc.model,    #The model we are evaluating
                    x = bioclim.data)    #Climatic variables for use by model

#Determine minimum threshold for "presence"
bc.threshold <- threshold(x = bc.eval, stat = "spec_sens")

#Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

#Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = "olivedrab")

#And add those observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "black",
       pch = "+", 
       cex = 0.75)

#Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

predict.presence > bc.threshold

#Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

#Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "olivedrab"))

#And add those observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "black",
       pch = "+", 
       cex = 0.75)

#Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()


