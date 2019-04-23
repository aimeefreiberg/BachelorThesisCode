# Analysis of weight of WATgon, WATsc and BAT
# (c) Danny Arends and Aimee Freiberg (HU Berlin), 2018 - 2024


# load data 
setwd("/Users/aimeefreiberg/Documents/Universtiy Bachelor/AG Brockmann/Bachelorarbeit/R_Data")
mdata <- read.table("Organs.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Genotype <- read.table("Genotype1.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Factors <- read.table("factors_AIL.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))

# bind data into one table 
mdata <- cbind(mdata, Genotype = NA)
mdata[rownames(Genotype), "Genotype"] <- as.character(unlist(Genotype[,1]))

mdata <- cbind(mdata, Mother = NA, WG = NA, Sex = NA)
mdata[rownames(Factors), "Mother"] <- as.character(unlist(Factors[,1]))
mdata[rownames(Factors), "WG"] <- as.character(unlist(Factors[,2]))
mdata[rownames(Factors), "Sex"] <- as.character(unlist(Factors[,3]))

# calculate relative fat storage weight (in %)
mdata <- cbind(mdata, "Rel.WATgon" = unlist(mdata[, "WATgon"] / mdata["Schlachtgewicht"])*100)
mdata <- cbind(mdata, "Rel.WATsc" = unlist(mdata[, "WATsc"] / mdata["Schlachtgewicht"])*100)
mdata <- cbind(mdata, "Rel.BAT" = unlist(mdata[, "BAT"] / mdata["Schlachtgewicht"])*100)

# test for impact of the factors sex, mother and litter size on data 
  # test for Rel.WATgon
    # set up null hypothesis
GH0 <- lm(Rel.WATgon ~ Sex + as.numeric(WG) + Mother, data=mdata)
AIC(GH0)

    # set up alternative hypothesis 
      # Genotypes coded as dominant/recessive due to knowledge from literature 

genotypes_num = as.numeric(factor(mdata[,"Genotype"], levels=c("CC", "CT", "TT")))
names(genotypes_num) = rownames(mdata)
genotypes_dom = as.numeric(factor(mdata[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2  
GHAdom <- lm(Rel.WATgon ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data=mdata)
AIC(GHAdom)

    # null hypothesis rejected, when AIC(H0) > AIC (HAdom)+20 <- 
AIC(GH0, GHAdom)
    # NULL HYPTHESIS REJECTED

    # set up second alternative hypothesis (remove mothers) 
GHAdomNoM <- lm(Rel.WATgon ~ Sex + as.numeric(WG) + genotypes_dom, data=mdata)

    # HAdom rejected, when AIC(HAdom) > AIC (HAdomNoM)+20 <- 
AIC(GHAdomNoM, GHAdom)
    # HAdom rejected in favor of HAdomNoM

    # check normed model HAdomNoM for normal distribution
qqnorm(GHAdomNoM$residuals)
qqline(GHAdomNoM$residuals, col="red")
shapiro.test(GHAdomNoM$residuals)
    # DATA IS NORMALLY DISTRIBUTED
    
    # see which factors are significantly different 
anova(GHAdomNoM)
    # ALL FACTORS SIGNIFICANTLY INFLUENCE REL.WATGON
    # WATgon significantly different between genotypes

    # set up model for plotting WATgon
Gplot <- lm(Rel.WATgon ~ Sex + as.numeric(WG), data=mdata)


  # test for Rel.WATsc
SH0 <- lm(Rel.WATsc ~ Sex + as.numeric(WG) + Mother, data=mdata)
AIC(SH0)

    # set up alternative hypothesis 
    # Genotypes coded as dominant/recessive due to knowledge from literature 

genotypes_dom = as.numeric(factor(mdata[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2  
SHAdom <- lm(Rel.WATsc ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data=mdata)
AIC(SHAdom)

    # null hypothesis rejected, when AIC(H0) > AIC (HAdom)+20 
AIC(SH0, SHAdom)
    # NULL HYPTHESIS REJECTED

    # set up second alternative hypothesis (remove mothers) 
SHAdomNoM <- lm(Rel.WATsc ~ Sex + as.numeric(WG) + genotypes_dom, data=mdata)

    # set uo thrids alternative hypothesis (remove littersize, because so significant influence)
SHAdomSex <- lm(Rel.WATsc ~ Sex + genotypes_dom, data=mdata)

    # check for best model again 
AIC(SHAdomNoM, SHAdom, SH0, SHAdomSex)
  # SHAdomNoMWG, because significant difference to H0 and least degrees of freedoms used

    # check normed model HAdomNoM for normal distribution
qqnorm(SHAdomSex$residuals)
qqline(SHAdomSex$residuals, col="red")
shapiro.test(SHAdomSex$residuals)
    # DATA IS NORMALLY DISTRIBUTED (with a few outliers)

    # see which factors are significantly different 
anova(SHAdomSex)
    # ALL FACTORS (BESIDES LITTER SIZE) SIGNIFICANTLY INFLUENCE REL.WATSC
    # WATsc significantly different between genotypes 

    # set up model for plotting WATgon
Splot <- lm(Rel.WATsc ~ Sex, data=mdata)


  # test for BAT
BH0 <- lm(Rel.BAT ~ Sex + as.numeric(WG) + Mother, data=mdata)
AIC(BH0)

# set up alternative hypothesis 
# Genotypes coded as dominant/recessive due to knowledge from literature 

genotypes_dom = as.numeric(factor(mdata[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2  
BHAdom <- lm(Rel.BAT ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data=mdata)
AIC(BHAdom)

    # null hypothesis rejected, when AIC(H0) > AIC (HAdom)+20 
AIC(BH0, BHAdom)
    # NULL HYPTHESIS ACCEPTED

    # set up second alternative hypothesis (remove mothers) 
BHAdomNoM <- lm(Rel.BAT ~ Sex + as.numeric(WG) + genotypes_dom, data=mdata)

    # set up third alternative hypothesis (without sex and WG)
BHAdomM <- lm(Rel.BAT ~ Mother + genotypes_dom, data=mdata)

    #set up fourth alternative hypothesis 
BHAdom0 <- lm(Rel.BAT ~ genotypes_dom, data=mdata)

    # null hypothesis rejected, when AIC(H0) > AIC (HAdom)+20 
AIC(BH0, BHAdom, BHAdomNoM, BHAdomM, BHAdom0)
    # no significant differences, choose BHAdom0 because of least degrees of freedom used 

    # check normed model HAdomNoM for normal distribution
qqnorm(BHAdom0$residuals)
qqline(BHAdom0$residuals, col="red")
shapiro.test(BHAdomNoM$residuals)
    # DATA IS NORMALLY DISTRIBUTED

    # see which factors are significantly different 
anova(BHAdom0)
    # BAT significantly different between genotypes 


# plot data with normed models 
  
par(mfrow=c(1,3))
  #WATgon
plot(x = c(0, 4), y=c(-1, 12), xaxt="n", main= "WATgon", ylab=" weight [%]", xlab="Genotype", t='n', las=2, cex.lab=2, cex.axis=1.5, cex.main=2)
boxplot(Gplot$residuals + mean(mdata[, "Rel.WATgon"], na.rm=TRUE) ~ mdata[names(Gplot$residuals), "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')
lines(x=c(1,2.5), y = c(11, 11))
lines(x=c(1,1), y = c(10, 11))
lines(x=c(2.5,2.5), y = c(11, 10.5))
lines(x=c(2,3), y = c(10.5, 10.5))
lines(x=c(3,3), y = c(10, 10.5))
lines(x=c(2,2), y= c(10, 10.5))
text(1.7, 11.3, paste0("***"), cex = 3)

  #WATsc
plot(x = c(0, 4), y=c(-1, 12), xaxt="n", main= "WATsc", ylab=" weight [%]", xlab="Genotype", t='n', las=2, cex.lab=2, cex.axis=1.5, cex.main=2)
boxplot(Splot$residuals + mean(mdata[, "Rel.WATsc"], na.rm=TRUE) ~ mdata[names(Splot$residuals), "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n') 
lines(x=c(1,2.5), y = c(6, 6))
lines(x=c(1,1), y = c(5, 6))
lines(x=c(2.5,2.5), y = c(6, 5.5))
lines(x=c(2,3), y = c(5.5, 5.5))
lines(x=c(3,3), y = c(5, 5.5))
lines(x=c(2,2), y= c(5, 5.5))
text(1.7, 6.3, paste0("***"), cex = 3)

  #BAT
plot(x = c(0, 4), y=c(0, 1), xaxt="n", main= "BAT", ylab=" weight [%]", xlab="Genotype", t='n', las=2, cex.lab=2, cex.axis=1.5, cex.main=2)
boxplot(mdata[, "Rel.BAT"] ~ mdata[, "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n') 
lines(x=c(1,2.5), y = c(0.8, 0.8))
lines(x=c(1,1), y = c(0.7, 0.8))
lines(x=c(2.5,2.5), y = c(0.8, 0.75))
lines(x=c(2,3), y = c(0.75, 0.75))
lines(x=c(3,3), y = c(0.7, 0.75))
lines(x=c(2,2), y= c(0.7, 0.75))
text(1.7, 0.83, paste0("***"), cex = 3)

# END
