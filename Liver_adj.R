# Analysis of different phenotypes of liver of AIL B6xS1 
# (c) Dany Arends, Aimee Freiberg (HU Berlin), 2018 - 2024

# load data
setwd("/Users/aimeefreiberg/Documents/Universtiy Bachelor/AG Brockmann/Bachelorarbeit/R_Data")
Triglycerides <- read.table("Liv_B6xS1.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Factors <- read.table("factors_AIL.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Genotype <- read.table("Genotype1.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
weights <- read.table("Organs.txt",  sep="\t", row.names=1, header=TRUE, na.strings = c("","NA") )

## calculate liver triglyceried contents 

#calculate tryglceride per protein content 
colnames(Triglycerides) <- c("Genotype","Trigly","Prot")

Trigly.Prot <- c()
Trigly.Prot <- Triglycerides$Trigly/Triglycerides$Prot
Triglycerides[, "Trigly/Prot"] <- Trigly.Prot

#bind factor data into one table

Triglycerides <- cbind(Triglycerides, Mother = NA, WG = NA, Sex = NA)
Triglycerides[rownames(Factors), "Mother"] <- as.character(unlist(Factors[,1]))
Triglycerides[rownames(Factors), "WG"] <- as.character(unlist(Factors[,2]))
Triglycerides[rownames(Factors), "Sex"] <- as.character(unlist(Factors[,3]))

  # set up model for Triglycerides
    # null hypothesis
TH0 <- lm(Trigly/Prot ~ Sex + as.numeric(WG) + Mother, data=Triglycerides)
AIC(TH0)

    # alternative hypothesis
genotypes_num = as.numeric(factor(Triglycerides[,"Genotype"], levels=c("CC", "CT", "TT")))
names(genotypes_num) = rownames(Triglycerides)
genotypes_dom = as.numeric(factor(Triglycerides[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2  
THAdom <- lm(Trigly/Prot ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data=Triglycerides)
AIC(THAdom)

# null hypothesis rejected, when AIC(TH0) > AIC (THAdom)+20 <- 
AIC(TH0, THAdom)
# NULL HYPTHESIS REJECTED

    # set up second alternative hypothesis (remove mothers) 
THAdomNoM <- lm(Trigly/Prot ~ Sex + as.numeric(WG) + genotypes_dom, data=Triglycerides)

    # set up third alternative hypothesis only genotypes
THAdomO <- lm(Trigly/Prot ~ genotypes_dom, data=Triglycerides)

    # Testing of models
AIC(TH0, THAdom, THAdomNoM, THAdomO)
    # THAdom0 chosen as best model due to least degrees of freedom used

  # check normed model HAdomNoM for normal distribution
qqnorm(THAdomO$residuals)
qqline(THAdomO$residuals, col="red")
shapiro.test(THAdomO$residuals)
  # data is not normally distributed    

#test for significance
  anova(THAdomO)
  #no significant difference between the genotypes
  

# test absolute liver weights
  
  # bind factor data into the tabplotp
weights <- cbind(weights, Mother = NA, WG = NA, Sex = NA, Genotypes = NA)
weights[rownames(Factors), "Mother"] <- as.character(unlist(Factors[,1]))
weights[rownames(Factors), "WG"] <- as.character(unlist(Factors[,2]))
weights[rownames(Factors), "Sex"] <- as.character(unlist(Factors[,3]))
weights[rownames(Genotype), "Genotype"] <- as.character(unlist(Genotype[,1]))

  # set up model 
LWH0 <- lm(Leber ~ Sex + as.numeric(WG) + Mother, data=weights)

  
  # set up alternative hypothesis 
  # Genotypes coded as dominant/recessive due to knowledge from literature 
  
  genotypes_dom = as.numeric(factor(weights[,"Genotype"], levels=c("CC", "CT", "TT")))
  genotypes_dom[genotypes_dom == 3] <- 2  
  LWHAdom <- lm(Leber ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data=weights)
  
  AIC(LWHAdom, LWH0)
  # Null hypothesis rejected, LWHAdom model chosen
  
  # check normed model HAdomNoM for normal distribution
  qqnorm(LWHAdom$residuals)
  qqline(LWHAdom$residuals, col="red")
  shapiro.test(LWHAdom$residuals)
  # a few outliers but generally normally distriubted 
  
  # see which factors are significantly different 
  anova(LWHAdom)
  # Liver weight significantly different between genotypes 
#model for plot 
  PLWHAdom <- lm(Leber ~ Sex + as.numeric(WG) + Mother, data=weights)  
  
# test relative liver weights (with lean body weight)
 
   #calculate relative liver weights 
weights <- cbind(weights, "leanBW" = unlist(weights[, "Schlachtgewicht"] - weights[, "WATgon"] - weights["WATsc"])) 
weights <- cbind(weights, "rLiver" = unlist(weights[, "Leber"] / weights["leanBW"])*100)  
  
# set up model for relative liver weight
  
  # null hypothesis

rLH0 <- lm(rLiver ~ Sex + as.numeric(WG) + Mother, data=weights)
AIC(rLH0)

  # alternative hypothesis
  genotypes_num = as.numeric(factor(weights[,"Genotype"], levels=c("CC", "CT", "TT")))
  names(genotypes_num) = rownames(weights)
  genotypes_dom = as.numeric(factor(weights[,"Genotype"], levels=c("CC", "CT", "TT")))
  genotypes_dom[genotypes_dom == 3] <- 2  
  rLdom <- lm(rLiver ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data=weights)
  AIC(rLdom)
  
  #set up second alternative hypothesis
  
  rLNoWG <- lm(rLiver ~ Sex + Mother + genotypes_dom, data=weights)
  
  # set up third alternative hypothesis
  
  rLS <- lm(rLiver ~ Sex + genotypes_dom, data=weights)

  #test for best model
  AIC(rLH0,rLdom, rLNoWG, rLS)
  # no model is significantly better (although lrS is worse)
  #model rLNoWG chosen

  # check normed model rLNoWG for normal distribution
  qqnorm(rLNoWG$residuals)
  qqline(rLNoWG$residuals, col="red")
  shapiro.test(rLNoWG$residuals)
  # a few outliers but generally normally distriubted 
  
  #test for significant differences  
  anova(rLNoWG)  
  # so significant differences between genotypes in realtive liver weights
  
  # set up model for graph 
  rLplot <- lm(rLiver ~ Sex + Mother, data=weights)
  
# Make liver graph

par(mfrow=c(1,3))
par (mar= c(7,5,7,5))  

# absolute liver weight
plot(x = c(0, 4), y=c(0, 4), xaxt="n", main= "absolute weight", ylab=" Weight [g]", xlab="Genotype", t='n', las=2,  cex.axis=1.5, cex.lab=2, cex.main=2)
boxplot(PLWHAdom$residuals + mean(weights[, "Leber"], na.rm=TRUE) ~ weights[names(PLWHAdom$residuals), "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')
lines(x=c(1,2.5), y = c(3.50, 3.50))
lines(x=c(1,1), y = c(3, 3.50))
lines(x=c(2.5,2.5), y = c(3.50, 3.3))
lines(x=c(2,3), y = c(3.3, 3.3))
lines(x=c(3,3), y = c(3, 3.3))
lines(x=c(2,2), y= c(3, 3.3))
text(1.7, 3.7, paste0("***"), cex = 3)


# relative liver weight 

plot(x = c(0, 4), y=c(2, 8.5), xaxt="n", main= "relative weight", ylab=" Weight [%]", xlab="Genotype", t='n', las=2,  cex.axis=1.5, cex.lab=2, cex.main=2)
boxplot(rLplot$residuals + mean(weights[, "rLiver"], na.rm=TRUE) ~ weights[names(rLplot$residuals), "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')

#trigly/prot   
plot(x = c(0, 4), y=c(0, 460), xaxt="n", main="triglyceride concentration ", ylab="Triglycerides (µg) /Protein (µg)", xlab="Genotype", t='n', las=2, cex.lab=2, cex.axis=1.2, cex.main=2)
boxplot(Triglycerides[, "Trigly/Prot"] ~ Triglycerides[, "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')
lines(x=c(1,2.5), y = c(450, 450))
lines(x=c(1,1), y = c(430, 450))
lines(x=c(2.5,2.5), y = c(450, 440))
lines(x=c(2,3), y = c(440, 440))
lines(x=c(3,3), y = c(430, 440))
lines(x=c(2,2), y= c(430, 440))
text(1.7, 460, paste0("."), cex = 8)

# extract statistical data

KörperlängeADj = THAdomO$residuals + mean(Triglycerides[,"Trigly/Prot"], na.rm=T)
weights = cbind(weights, KörperlängeADj = NA)
weights[names(weights), "Körperlängeadj"] = KörperlängeADj

KörperlängeADj = THAdomO$residuals + mean(Triglycerides[,"Trigly/Prot"], na.rm=T)
weights = cbind(weights, KörperlängeADj = NA)
weights[names(KörperlängeADj), "Körperlängeadj"] = KörperlängeADj

KLCC <- Triglycerides[which(Triglycerides[,"Genotype"] == "CC"), 4]
KLCT <- Triglycerides[which(Triglycerides[,"Genotype"] == "CT"), 4]
KLTT <- Triglycerides[which(Triglycerides[,"Genotype"] == "TT"), 4]

columns <- c("KLCC","KLCT","KLTT") 
meanKL <- c()
sdKL <- c()
for(column in columns){
  meanKL <- c(meanKL, mean(get(column), na.rm = T))
  sdKL <- c(sdKL, sd(get(column), na.rm = T))
}

#END