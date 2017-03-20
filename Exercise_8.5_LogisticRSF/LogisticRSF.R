library(lme4)
library(AICcmodavg)
library(adehabitatHR)

#Load text files for each season
data_1 <- read.csv("MD_winter12.csv",header=T)
str(data_1)

############################### MODELING ###############################
data_1$crop=as.factor(data_1$crop)
#dataset$SEASON=factor(dataset$SEASON)
data_1[,2:3]=scale(data_1[,2:3],scale=TRUE)#standardize data to mean of zero
str(data_1)

#Models considered
fit1 = glmer(use ~ relevel(crop,"1")+(1|animal_id), data=data_1, family=binomial(link="logit"),nAGQ = 0)#Sunflower and cover model
fit2 = glmer(use ~ d_cover+(1|animal_id), data=data_1, family=binomial(link="logit"),nAGQ = 0)#Distance to cover only model
fit3 = glmer(use ~ d_roads+(1|animal_id), data=data_1, family=binomial(link="logit"),nAGQ = 0)#Distance to roads only model
fit4 = glmer(use ~ d_cover+d_roads+(1|animal_id), data=data_1, family=binomial(link="logit"),nAGQ = 0)#Distance to cover and roads model
fit5 = glmer(use ~ 1|animal_id, data=data_1, family=binomial(link="logit"),nAGQ = 0)#Intercept model

fit1
fit2
fit3
fit4
fit5

AIC(fit1,fit2,fit3,fit4,fit5)

mynames <- paste("fit", as.character(1:5), sep = "")
myaicc <- aictab(list(fit1,fit2,fit3,fit4,fit5), modnames = mynames)
print(myaicc, LL = FALSE)

#Get confidence intervals from the top model to interpret results 
per1_se <- sqrt(diag(vcov(fit1)))
# table of estimates with 95% CI
tab_per1 <- cbind(Est = fixef(fit1), LL = fixef(fit1) - 1.96 * per1_se, UL = fixef(fit1) + 1.96 * per1_se)

################################################################################
#############################################################################
#CODE BELOW IS TO CREATE PREDICTIVE SURFACES FOR THE TOP MODEL FOR SUMMER 2012
#i.e, fit1).
################################################################################
################################################################################
layer1 <- read.table("layer1.txt",sep=",")
str(layer1)
names(layer1) = c("crop", "d_cover", "d_roads","x", "y")
str(layer1)
head(layer1)
#Need to standardize the raw distance rasters first to match what we modeled
layer1[,2:3]=scale(layer1[,2:3],scale=TRUE)
head(layer1)
layer1$crop <- as.factor(layer1$crop)

# predictions based on best model 
predictions = predict(fit1, newdata=layer1, re.form=NA, type="link")	# based on the scale of the linear predictors
predictions = exp(predictions)
range(predictions)

#-----------------------------------------------------------------------------------
# create Ascii grid of raw predictions if needed
layer1$predictions = predictions
#preds = layer1
#preds = SpatialPixelsDataFrame(points=preds[c("x", "y")], data=preds)
#preds = as(preds, "SpatialGridDataFrame")
#names(preds)
#writeAsciiGrid(preds, "predictions.asc", attr=13) # attr should be column number for 'predictions'

#-----------------------------------------------------------------------------------
# assign each cell or habitat unit to a 'prediction class'.
# classes have (nearly) equal area, if the cells or habitat units have equal areas.
# output is a vector of class assignments (higher is better).
F.prediction.classes <- function(raw.prediction, n.classes){
  # raw.prediction = vector of raw (or scaled) RSF predictions
  # n.classes = number of prediction classes.
  pred.quantiles = quantile(raw.prediction, probs=seq(1/n.classes, 1-1/n.classes, by=1/n.classes))
  ans = rep(n.classes, length(raw.prediction))
  for(i in (n.classes-1):1){
    ans[raw.prediction < pred.quantiles[i]] = i
  }
  return(ans)
}

str(layer1)
layer1$prediction.class = F.prediction.classes(layer1$predictions, 5)

table(layer1$prediction.class)

##############################################
# create map of RSF prediction classes in R
m = SpatialPixelsDataFrame(points = layer1[c("x", "y")], data=layer1)
names(m)
par(mar=c(0,0,0,0))
image(m, attr=7, col=c("grey90", "grey70", "grey50", "grey30", "grey10"))
par(lend=1)
legend("bottomright", col=rev(c("grey90", "grey70", "grey50", "grey30", "grey10")),
       legend=c("High", "Medium-high", "Medium", "Medium-low", "Low"),
       title="Prediction Class", pch=15, cex=1.0,bty != "n", bg="white")

# create Ascii grid of prediction classes
#m = as(m, "SpatialGridDataFrame")
#names(m)
#writeAsciiGrid(m, "PredictionClassess.asc", attr=7)


