# Calculating BMI & Body height 
# (c) Aimee Freiberg (HU Berlin), 2018 - 2024

# load data 
setwd("/Users/aimeefreiberg/Documents/Universtiy Bachelor/AG Brockmann/Bachelorarbeit/R_Data")
mdata <- read.table("Organs.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Genotype <- read.table("Genotype1.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Factors <- read.table("factors_AIL.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))

# calculate BMI
surface <- c()
  #use dubois equation to calculate body suface(m^2)
surface <- (0.007184*((mdata[, "Schlachtgewicht"]/1000)^0.425)*(mdata[, "Körperlänge"]^0.725))*1000
  # calculate BMI with g/m^2
BMI <- cbind(surface, "BMI" = unlist(mdata[, "Schlachtgewicht"] / surface))

# cbind Genotype and factor data
BMI <- as.data.frame.matrix(BMI) 
BMI <- cbind (BMI, Genotype = as.character(unlist(Genotype[,1])), Mother = as.character(unlist(Factors[,1])), WG = as.character(unlist(Factors[,2])), Sex = as.character(unlist(Factors[,3])))

# adjust data for factors
  # null hypothesis
BH0 <- lm(BMI ~ Sex + as.numeric(WG) + Mother, data=BMI)
AIC(BH0)
  
  #adjust for dominance effect
genotypes_num = as.numeric(factor(BMI[,"Genotype"], levels=c("CC", "CT", "TT")))
names(genotypes_num) = rownames(BMI)
genotypes_dom = as.numeric(factor(BMI[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2 

 BHgen <- lm(BMI ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data = BMI)

 #test which model is better
 AIC(BH0,BHgen)
 #BHgen is better than BH0 <- chosen
 
 # make model for plot
 
 Bplot <- lm(BMI ~ Sex + as.numeric(WG) + Mother, data = BMI)
 
 #check for significant differences
 anova(BHgen)
 # significant differences between genotypes
 
# make plot
 
 plot(x = c(0, 4), y=c(2, 6), xaxt="n", main= "Body Mass Index", ylab=" BMI [g/cm^2]", xlab="Genotype", t='n', las=2, cex.lab=1.2, cex.axis=1.2, cex.main=2)
 boxplot(Bplot$residuals + mean(BMI[, "BMI"], na.rm=TRUE) ~ BMI[names(Bplot$residuals), "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n') 
 lines(x=c(1,2.5), y = c(5.5, 5.5))
 lines(x=c(1,1), y = c(5, 5.5))
 lines(x=c(2.5,2.5), y = c(5.5, 5.2))
 lines(x=c(2,3), y = c(5.2, 5.2))
 lines(x=c(3,3), y = c(5, 5.2))
 lines(x=c(2,2), y= c(5, 5.2))
 text(1.7,5.7, paste0("***"), cex = 3)
 # BMI of CC significantly higher 
 
 # extract statistics from data
 
 bmiADj = Bplot$residuals + mean(BMI[,"BMI"], na.rm=T)
 BMI = cbind(BMI, AdjBMI = NA)
 BMI[names(bmiADj), "AdjBMI"] = bmiADj
 
 BMICC <- BMI[which(BMI[,"Genotype"] == "CC"), 7]
 BMICT <- BMI[which(BMI[,"Genotype"] == "CT"), 7]
 BMITT <- BMI[which(BMI[,"Genotype"] == "TT"), 7]
 
 columns <- c("BMICC","BMICT","BMITT") 
 meanBMI <- c()
 sdBMI <- c()
 for(column in columns){
 meanBMI <- c(meanBMI, mean(get(column), na.rm = T))
 sdBMI <- c(sdBMI, sd(get(column), na.rm = T))
 }

 
 
 
###### Body height ##########
 
 # adjust data for factors
 # null hypothesis
 HH0 <- lm(Körperlänge ~ Sex + as.numeric(WG) + Mother, data=mdata)
 
 #bind factor data
 
 mdata <- cbind (mdata, Genotype = as.character(unlist(Genotype[,1])), Mother = as.character(unlist(Factors[,1])), WG = as.character(unlist(Factors[,2])), Sex = as.character(unlist(Factors[,3])))
 
  #adjust for dominance effect
genotypes_num = as.numeric(factor(BMI[,"Genotype"], levels=c("CC", "CT", "TT")))
names(genotypes_num) = rownames(BMI)
genotypes_dom = as.numeric(factor(BMI[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2 
 
Hgen <- lm(Körperlänge ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data = mdata)


 # second alternative hypothesis 

HNoM <- lm(Körperlänge ~ Sex + genotypes_dom, data = mdata)
# check for the best model
AIC(HH0, Hgen, HNoM)
#No model is better, but HNoM uses the least degrees of freedom

anova(HNoM)
#all factors are significant 
# genotypes differ significantly from each other

# make model for graph

Hgenplot <- lm(Körperlänge ~ Sex, data = mdata)

#make plot 

plot(x = c(0, 4), y=c(8, 14), xaxt="n", ylab=" height [cm]", xlab="Genotype", t='n', las=2, cex.lab=1.2, cex.axis=1.2, cex.main=2)
boxplot(Hgenplot$residuals + mean(mdata[, "Körperlänge"], na.rm=TRUE) ~ mdata[names(Hgenplot$residuals), "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n') 
lines(x=c(1,2.5), y = c(13, 13))
lines(x=c(1,1), y = c(12,13))
lines(x=c(2.5,2.5), y = c(13, 12.5))
lines(x=c(2,3), y = c(12.5, 12.5))
lines(x=c(3,3), y = c(12, 12.5))
lines(x=c(2,2), y= c(12, 12.5))
text(1.7,13.2, paste0("*"), cex = 3)

# get statistical data

KörperlängeADj = Hgenplot$residuals + mean(mdata[,"Körperlänge"], na.rm=T)
mdata = cbind(mdata, KörperlängeADj = NA)
mdata[names(KörperlängeADj), "Körperlängeadj"] = KörperlängeADj

KLCC <- mdata[which(mdata[,"Genotype"] == "CC"), 28]
KLCT <- mdata[which(mdata[,"Genotype"] == "CT"), 28]
KLTT <- mdata[which(mdata[,"Genotype"] == "TT"), 28]

columns <- c("KLCC","KLCT","KLTT") 
meanKL <- c()
sdKL <- c()
for(column in columns){
  meanKL <- c(meanKL, mean(get(column), na.rm = T))
  sdKL <- c(sdKL, sd(get(column), na.rm = T))
}

#END