# clear out the past 
if(!is.null(dev.list())) dev.off()  
rm(list = ls())
cat("\014")

#################PACKAGES#####################

library(dplyr)
library(tidyverse)
library(corrplot)
library(gridExtra)
library(grid)
library(tidyr) # for splitting on the period see below
library(moments) # for calculating moments for skewness etc.
library(reshape2)
par(mfrow=c(1, 1)) 

#################DATA FETCH###################

setwd("C:/Files/UoE/Modules/Spring/MA334-7-SP Data analysis and statistics with R/Assignment")
proj_data <- read.csv("proportional_species_richness_V2.csv") #reading the V2 data
proj_data_nona <- read.csv("proportional_species_richness_V3.csv") #reading the V3 data
#proj_data <- proj_data[, -c(3, 8, 9, 10)] #selecting only the allocated 7 taxonomic groups
names(proj_data)
eco_selected <- c(2,4,5,6,7,11,12) #selecting only the allocated 7 taxonomic groups


#uni-variate exploration 

summary(proj_data_nona$Bryophytes)

g1 <- ggplot(proj_data_nona, aes(x = Bees)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")

g2 <- ggplot(proj_data_nona, aes(x = Bryophytes)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")

g3 <- ggplot(proj_data_nona, aes(x = Butterflies)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")

g4 <- ggplot(proj_data_nona, aes(x = Carabids)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")

g5 <- ggplot(proj_data_nona, aes(x = Hoverflies)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")

g6 <- ggplot(proj_data_nona, aes(x = Grasshoppers_._Crickets)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")

g7 <- ggplot(proj_data_nona, aes(x = Vascular_plants)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")

grid.arrange(g1, g2, g3, g4, g5, g6, g7, nrow=4, top=textGrob("Distribution of the allocated taxonomic group"))

#correlation matrix of the selected taxonomic groups
cor_BD7 <- cor(proj_data_nona[ ,eco_selected])
corrplot(cor_BD7, method = "square", title = "Correlation of 7 taxonomic groups", 
         tl.col="black", tl.srt=45,  mar=c(0,0,4,3))

##################################################

# Select columns of interest & Calculating mean and SD 
subset_BD7 <- proj_data_nona %>%
  select(Location, Bees, Bryophytes, Butterflies, Carabids, Hoverflies, Grasshoppers_._Crickets, Vascular_plants)
var_stats <- subset_BD7 %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))

# Melting data for visualization
melt_data <- subset_BD7 %>%
  pivot_longer(cols = -Location, names_to = "groups", values_to = "Prop_species_richness")

ggplot(melt_data, aes(y = Prop_species_richness, x = groups)) +
  geom_boxplot() +
  ggtitle("Proportional species richness for different taxonomic groups") +
  ylab("Proportional species richness") +
  xlab("Allocated taxonomic groups")

##################################################
# count incidents (both periods) for each land classification for later selection
n1 <- proj_data%>%group_by(dominantLandClass)%>%count()%>%
  arrange(dominantLandClass)

ggplot(n1, aes(dominantLandClass, n)) +
  geom_histogram(stat = "identity", fill = "lightblue", color = "black", bins = 10) +
  labs(title = "Occurences of Species in different dominantLandClass ", x = "dominantLandClass", y = "count") +
  coord_flip()

# FILTERING VALUES BASED ON THE LOCATION, DOMINANTLANDCLASS, ETC
#proj_data%>%filter(grepl("TM",Location))
#proj_data%>%filter(grepl("w",dominantLandClass))
#proj_data%>%filter(grepl("TM",Location)|grepl("TG",Location))
#proj_data%>%filter(dominantLandClass=="3e")

# selecting allocated  7 chosen predictors to form the trial eco_stat
all <- c(2:12)
eco_not_selected <- all[!(all%in%eco_selected)]
eco_names <- names(proj_data[,2:12])
eco_selected_names <- names(proj_data)[eco_selected]
eco_selected_names

# calculate the bio-div measure over 7 taxonomic groups
mean_selected <- rowMeans(proj_data[,eco_selected],na.rm=TRUE) # mean the 7 columns 
sum(is.na(mean_selected)) # checking that there are no NAs in mean_selected

# adding the biodiversity measure which is the mean over 7 taxonomic groups
proj_data_MA334 <- proj_data%>%mutate(eco_status_7=mean_selected)
names(proj_data_MA334)

##################DATA EXPLORATION###################

# splitting the data by period and comparing the stats before and after 
table <- data.frame()
for(i in eco_selected){
  table <- rbind(table,
                 c(eco_names[i-1],
                   round(mean(proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(sd(proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(skewness(proj_data_MA334[,i],na.rm = TRUE),digits = 2)
                 ))}
colnames(table) <- c("taxi_group","mean","sd","skewness")
table%>%arrange(sd,skewness) # something more could be done here


# extending data exploration; with correlations between continuous variables
names(proj_data_MA334)
cont_vars <- proj_data_MA334%>%select(c(eco_selected,16,17)) # includes easting and northing 
names(cont_vars)
cormat <- round(x = cor(cont_vars,use="pairwise.complete.obs"), digits = 2)
# melt the correlation matrix
meltcor <- melt(cormat)%>%mutate(R2 = value^2)%>%arrange(value)
meltcor %>%mutate(R2 = value^2)%>%arrange(Var1,value) #sorting with the columns

# plotting a map
plot(cont_vars$Northing~cont_vars$Easting, xlab = "Easting", ylab = "Northing", main = "The UK Map")

# now using the eastings and northings (these may be better used as predictors )
plot(proj_data_MA334$eco_status_7~proj_data_MA334$Easting) 
cor(proj_data_MA334$eco_status_7,proj_data_MA334$Easting)
plot(proj_data_MA334$eco_status_7~proj_data_MA334$Northing)
cor(proj_data_MA334$eco_status_7,proj_data_MA334$Northing)

# doing a linear regression with only Northing as a predictor 
lin_mod <- lm(proj_data_MA334$eco_status_7~proj_data$Northing)
summary(lin_mod)
abline(lin_mod,col="green")
plot(jitter(fitted(lin_mod)),residuals(lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0)
qqnorm(lin_mod$residuals)
qqline(lin_mod$residuals,col="red")

# following code splits between the two periods to find the BD7 change
# however it may be better to use period as a predictor

# box plot comparisons for the two periods ignoring all other variables 
eco_status <- proj_data_MA334%>%pull(eco_status_7)
eco_period <- proj_data_MA334%>%pull(period)
ggplot(proj_data_MA334, aes(period, eco_status_7)) + 
  geom_boxplot()

names(proj_data_MA334)
proj_data_MA334_period <- proj_data_MA334%>%select(Location,period,eco_status_7)
proj_data_MA334_split <- proj_data_MA334_period%>%pivot_wider(names_from =period,values_from=eco_status_7)
proj_data_MA334_split <- proj_data_MA334_split%>%mutate(BD7_change=Y00-Y70)
head(proj_data_MA334_split)
hist(proj_data_MA334_split$BD7_change, xlab = "BD7_change", main = "Distribution of the BD7 change") # the distribution of the BD7 change 


BD7_change <- proj_data_MA334_split%>%pull(BD7_change)
# t test with H0: mu=0
t.test(BD7_change,mu=0)  


# comparing the two distributions of bio div based on 7 and 11 taxonomic groups 
par(mfrow=c(1, 1))  # divide graph area in 1 columns
qqplot(proj_data_MA334$eco_status_7,proj_data_MA334$ecologicalStatus)
abline(0,1,col="red")
# both cdfs together  and do a kolmogorov test H0: distributions are the same
BD7_cdf <- ecdf(proj_data_MA334$eco_status_7)
BD11_cdf <- ecdf(proj_data_MA334$ecologicalStatus)
plot(BD11_cdf,col="red")
lines(BD7_cdf,col="green")
ks.test(proj_data_MA334$eco_status_7,proj_data_MA334$ecologicalStatus)


# Simple linear regression
# regressions of eco_status_7 against ecologicalstatus based on all 11
plot(proj_data_MA334$eco_status_7~proj_data_MA334$ecologicalStatus)
abline(0,1,col="red")
lin_mod <- lm(proj_data_MA334$eco_status_7~proj_data_MA334$ecologicalStatus)
abline(lin_mod,col="green")
plot(jitter(fitted(lin_mod)),residuals(lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0,col="blue")
qqnorm(residuals(lin_mod))
qqline(residuals(lin_mod),col="red")
# do the same for each period report and differences 
proj_data_MA334_Y70 <- proj_data_MA334%>%filter(period=="Y70")
lin_mod <- lm(proj_data_MA334_Y70$eco_status_7~proj_data_MA334_Y70$ecologicalStatus)
lin_mod$coefficients
# for later period 
proj_data_MA334_Y00 <- proj_data_MA334%>%filter(period=="Y00")
lin_mod <- lm(proj_data_MA334_Y00$eco_status_7~proj_data_MA334_Y00$ecologicalStatus)
lin_mod$coefficients


# linear regression of BD4 on BD7 
mean_selected <- rowMeans(proj_data[,eco_not_selected ],na.rm=TRUE) # mean the rem 4 columns 
sum(is.na(mean_selected)) # checking that there are no NAs in mean_selected
# add in the biodiversity measure which is the mean over 7 taxonomic groups
proj_data_MA334 <- proj_data_MA334%>%mutate(eco_status_4=mean_selected)
names(proj_data_MA334)

# regressions of means: eco_status_4 against others not inc eco_status_4 data
plot(proj_data_MA334$eco_status_4~proj_data_MA334$eco_status_7)
abline(0,1,col="red")
lin_mod <- lm(proj_data_MA334$eco_status_4~proj_data_MA334$eco_status_7)
summary(lin_mod)
abline(lin_mod,col="green")
plot(jitter(fitted(lin_mod)),residuals(lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0,col="blue")
qqnorm(residuals(lin_mod))
qqline(residuals(lin_mod),col="red")

# multiple linear regression BD4 against the selected 7 

# Create Training and Test data 
trainingRowIndex <- sample(1:nrow(proj_data_MA334), 0.8*nrow(proj_data_MA334))  # row indices for 80% training data
trainingData <- proj_data_MA334[trainingRowIndex, ]  # model training data
testData  <- proj_data_MA334[-trainingRowIndex, ]%>%na.omit # for test data remove NAs 

# Build the model on training data
lmMod_train <- lm(eco_status_4~.,
                  data=trainingData[c(eco_selected_names,"eco_status_4")],
                  na.action=na.omit,y=TRUE)
summary (lmMod_train)  # model summary
cor(lmMod_train$fitted.values,lmMod_train$y) # cor training data 
Eco_4_Pred <- predict(lmMod_train, testData) # predict to check model on test Data
cor(Eco_4_Pred,testData$eco_status_4)
plot(Eco_4_Pred~testData$eco_status_4)
abline(0,1,col="red")
# mis_fit_to_testData are the residuals for the train model fit to the test data 
mis_fit_to_testData <- testData$eco_status_4-Eco_4_Pred
plot(mis_fit_to_testData~Eco_4_Pred) # look for unwanted pattern in residuals
abline(0,0,col="red")
qqnorm(mis_fit_to_testData) # check for normality of residuals in prediction
qqline(mis_fit_to_testData,col="red")

step(lmMod_train, direction = "both")


# multiple linear regression BD7 against period, easting and northing 
mult_lin_mod <- lm(eco_status_7~.,
                   data=proj_data_MA334[c("eco_status_7",
                                          "period","Easting","Northing")],
                   na.action = na.omit,y=TRUE)
summary(mult_lin_mod)
plot(mult_lin_mod$fitted.values~mult_lin_mod$y)
abline(0,1,col="red")
plot(jitter(fitted(mult_lin_mod)),residuals(mult_lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0,col="blue")
qqnorm(residuals(mult_lin_mod))
qqline(residuals(mult_lin_mod),col="red")


# comparing the effect of each significant coefficient to that of period
mult_lin_mod$coefficients
as.numeric(mult_lin_mod$coefficients[3])*mean(proj_data_MA334$Easting)
as.numeric(mult_lin_mod$coefficients[4])*mean(proj_data_MA334$Northing)

# The following PCA method is an extension to the set book 
# PCA for visualizing the multi-dimensional spread of biodiversity values #######################

table(proj_data_MA334_Y70$period); table(proj_data_MA334_Y00$period) # check that these separate periods 
table(proj_data_MA334_Y00$Location==proj_data_MA334_Y70$Location) # check that Locations correspond between the two periods

eco_difference <- proj_data_MA334_Y00[,eco_selected ]-proj_data_MA334_Y70[,eco_selected ] # general differences between the two periods 
head(eco_difference)

################################################

BD7_3e_y70 <- proj_data_MA334 %>%
  filter(dominantLandClass == "3e", period == "Y70") %>%
  group_by(eco_status_7, dominantLandClass, period) %>%
  summarize()

BD7_3e_y00 <- proj_data_MA334 %>%
  filter(dominantLandClass == "3e", period == "Y00") %>%
  group_by(eco_status_7, dominantLandClass, period) %>%
  summarize()


BD7_3e <- rbind(data.frame(period = "Y00", eco_status_7 = BD7_3e_y00$eco_status_7),
            data.frame(period = "Y70", eco_status_7 = BD7_3e_y70$eco_status_7))

ggplot(BD7_3e, aes(x=eco_status_7, fill=period)) +
  geom_histogram(binwidth = 0.01, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  labs(title = "BD7 in E Anglia/S England for different periods",
       x = "BD7", y = "Occurrences")

###################YET TO ADD & INTERPRET###################################


# see ?prcomp the default here is the mean correct but not to scale 
pr.out=prcomp(na.omit(eco_difference)) # Principal Components 
pr.out$center  # gives the mean corrections the "centers"
pr.out$scale  # not scaled
pr.out$rotation[,1:2] # print out first two principal axes
screeplot(pr.out, type="lines") # plot the variances in decreasing order
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], proj_data_MA334_Y00$dominantLandClass, cex=0.5, pos=4, col="red") # location labels

# label by location 
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], proj_data_MA334_Y00$Location, cex=0.4, pos=4, col="red") # location labels

# label by eco increase 
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
BD_inc  <- proj_data_MA334_Y00$eco_status_7-proj_data_MA334_Y70$eco_status_7 # BD differences between the two periods 
text(pr.out$x[,1],pr.out$x[,2], round(BD_inc,2), cex=0.4, pos=4, col="red") # location labels

# label by a particular taxi group (if any) dominant in the first PC 
eco_selected1 <- c(12,2,5,6,3,7,4) # an example Bees !
eco_difference <- proj_data_MA334_Y00[,eco_selected1 ]-proj_data_MA334_Y70[,eco_selected1 ] # general differences between the two periods 
pr.out=prcomp(na.omit(eco_difference)) # Principal Components 
pr.out$rotation[,1:2] # print out first two principal axes
screeplot(pr.out, type="lines") # plot the variances in decreasing order
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], round(eco_difference$Bees,2), cex=0.5, pos=4, col="red") # location labels

