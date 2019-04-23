# Analysis of weight of body weight and body weight change
# (c) Danny Arends and Aimee Freiberg (HU Berlin), 2018 - 2024

# load data 
setwd("/Users/aimeefreiberg/Documents/Universtiy Bachelor/AG Brockmann/Bachelorarbeit/R_Data")
WW <- read.table("WWeight.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Genotype <- read.table("Genotype1.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Factors <- read.table("factors_AIL.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))

# add genotype data
WW <- cbind(WW, Genotype = NA, Sex = NA, WG = NA, Mother = NA)
WW[rownames(Genotype),"Genotype"] <- as.character(unlist(Genotype[,1]))
WW[rownames(Factors), "Mother"] <- as.character(unlist(Factors[,1]))
WW[rownames(Factors), "WG"] <- as.character(unlist(Factors[,2]))
WW[rownames(Factors), "Sex"] <- as.character(unlist(Factors[,3]))

# factor in dominance effect 
genotypes_num = as.numeric(factor(WW[,"Genotype"], levels=c("CC", "CT", "TT")))
names(genotypes_num) = rownames(WW)
genotypes_dom = as.numeric(factor(WW[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2  

# test for significant differences
columns <- c("X21","X28","X35","X42","X49","X56","X63","X70","X77","X84","X91","X98","X105","X112","X119","X126")
anova.pvals.absolute <- c()
shapiro.pvals.absolute <- c()
WWadjusted = matrix(NA, nrow(WW), length(columns), dimnames=list(rownames(WW), columns))
for(column in columns){
  histo.model <- lm(as.numeric(WW[, column]) ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data = WW)
  histo.model2 <- lm(as.numeric(WW[, column]) ~ Sex + as.numeric(WG) + Mother, data = WW)
  adjusted = histo.model2$residuals + mean(as.numeric(WW[, column]), na.rm = T)
  WWadjusted[names(adjusted),column] = adjusted 
  #Test for differences in S1 vs S2
  anova.pvals.absolute <- c(anova.pvals.absolute, anova(histo.model)[[5]][4])
  # Test if we can believe the anova test !
  shapiro.pvals.absolute <- c(shapiro.pvals.absolute, shapiro.test(histo.model$residuals)$"p.value")
}

names(anova.pvals.absolute) <- columns
names(shapiro.pvals.absolute) <- columns
round(anova.pvals.absolute,4)
#body weight significantly different from 7 weeks  

# check for siginificant impact of factors
anova(histo.model) # All have significant impact 

#bind genotype data
WWadjusted <- as.data.frame(WWadjusted)
WWadjusted <- cbind(WWadjusted, Genotype = NA)
WWadjusted[rownames(Genotype),"Genotype"] <- as.character(unlist(Genotype[,1]))

# seperate for genotypes 

indCC <- WWadjusted[which(WWadjusted[,"Genotype"] == "CC"),1:16]
indCT <- WWadjusted[which(WWadjusted[,"Genotype"] == "CT"),1:16]
indTT <- WWadjusted[which(WWadjusted[,"Genotype"] == "TT"),1:16]

#plot weekly bodyweight change 

plot(x= c(21, 126), y = c(10, 55), t = 'n')
apply(indCC, 1,function(x){ points(x = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126), y = as.numeric(x), col="green", t = 'l') })
apply(indCT, 1,function(x){ points(x = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126), y = as.numeric(x), col="blue", t = 'l') })
apply(indTT, 1,function(x){ points(x = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126), y = as.numeric(x), col="orange", t = 'l') })

statsCC <- apply(indCC, 2, function(x){ c(mean(x, na.rm = TRUE), sqrt(sd(x, na.rm = TRUE)), sd(x, na.rm = TRUE))  })
statsCT <- apply(indCT, 2, function(x){ c(mean(x), sqrt(sd(x)), sd(x)) })
statsTT <- apply(indTT, 2, function(x){ c(mean(x), sqrt(sd(x)), sd(x)) })

plot(x= c(21,126),y= c(10,50), ylab = "body weight (g)", xlab= "time (weeks)", t= 'n', xaxt="n")
axis(1, at = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126), c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126) / 7)
polygon(x = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126, 126, 119, 112, 105, 98, 91, 84, 77, 70, 63, 56, 49, 42, 35, 28, 21), c(statsCC[1,] + statsCC[2,], (statsCC[1,] - statsCC[2,])[16:1]), col = rgb(255, 165, 0, 125, maxColorValue=255), border=NA)
points(x = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126), statsCC[1,], t = 'l', col="orange")

polygon(x = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126, 126, 119, 112, 105, 98, 91, 84, 77, 70, 63, 56, 49, 42, 35, 28, 21), c(statsCT[1,] + statsCT[2,], (statsCT[1,] - statsCT[2,])[16:1]), col = rgb(0, 165, 255, 125, maxColorValue=255), border=NA)
points(x =c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126), statsCT[1,], t = 'l', col="white")

polygon(x = c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126, 126, 119, 112, 105, 98, 91, 84, 77, 70, 63, 56, 49, 42, 35, 28, 21), c(statsTT[1,] + statsTT[2,], (statsTT[1,] - statsTT[2,])[16:1]), col = rgb(0, 255, 0, 125, maxColorValue=255), border=NA)
points(x =c(21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126), statsTT[1,], t = 'l', col="lightgreen")

legend("topleft", c("CC", "CT", "TT"), fill=c("orange", rgb(0, 165, 255, 125, maxColorValue=255), "lightgreen"), border=NA)
text(28,24, paste0("*"), cex = 3)
text(49,32, paste0("*"), cex = 3)
text(56,34, paste0("*"), cex = 3)
text(63,36, paste0("*"), cex = 3)
text(70,37, paste0("*"), cex = 3)
text(77,38, paste0("*"), cex = 3)
text(84,39, paste0("*"), cex = 3)
text(91,40.5, paste0("*"), cex = 3)
text(98,41.5, paste0("*"), cex = 3)
text(105,42.5, paste0("*"), cex = 3)
text(112,43.5, paste0("*"), cex = 3)
text(119,44, paste0("*"), cex = 3)
text(126,45, paste0("*"), cex = 3)





######### calculate weekly weight gain ##############
weightgain <- c()
weightgain <- c(W1 = (WW[,2] - WW[,1]))
weightgain <- cbind(weightgain, W2 = (WW[,3] - WW[,2]))
weightgain <- cbind(weightgain, W3 = (WW[,4] - WW[,3]))
weightgain <- cbind(weightgain, W4 = (WW[,5] - WW[,4]))
weightgain <- cbind(weightgain, W5 = (WW[,6] - WW[,5]))
weightgain <- cbind(weightgain, W6 = (WW[,7] - WW[,6]))
weightgain <- cbind(weightgain, W7 = (WW[,8] - WW[,7]))
weightgain <- cbind(weightgain, W8 = (WW[,9] - WW[,8]))
weightgain <- cbind(weightgain, W9 = (WW[,10] - WW[,9]))
weightgain <- cbind(weightgain, W10 = (WW[,11] - WW[,10]))
weightgain <- cbind(weightgain, W11 = (WW[,12] - WW[,11]))
weightgain <- cbind(weightgain, W12 = (WW[,13] - WW[,12]))
weightgain <- cbind(weightgain, W13 = (WW[,14] - WW[,13]))
weightgain <- cbind(weightgain, W14 = (WW[,15] - WW[,14]))
weightgain <- cbind(weightgain, W15 = (WW[,16] - WW[,15]))
colnames(weightgain)[colnames(weightgain)=="weightgain"] <- "W1"

weightgain <- as.data.frame.matrix(weightgain) 

# add genotype data + factors

weightgain <- cbind(weightgain, Genotype = NA, Sex = NA, WG = NA, Mother = NA)
rownames(weightgain) <- rownames(Genotype)
weightgain[rownames(Genotype),"Genotype"] <- as.character(unlist(Genotype[,1]))
weightgain[rownames(Factors), "Mother"] <- as.character(unlist(Factors[,1]))
weightgain[rownames(Factors), "WG"] <- as.character(unlist(Factors[,2]))
weightgain[rownames(Factors), "Sex"] <- as.character(unlist(Factors[,3]))

# adjust for dominance effect

genotypes_num = as.numeric(factor(weightgain[,"Genotype"], levels=c("CC", "CT", "TT")))
names(genotypes_num) = rownames(weightgain)
genotypes_dom = as.numeric(factor(weightgain[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2 

# test significant differences 

columns <- c("W1","W2","W3","W4","W5","W6","W7","W8", "W9", "W10", "W11", "W12", "W13", "W14","W15")
anova.pvals.absolute <- c()
shapiro.pvals.absolute <- c()
WGadjusted = matrix(NA, nrow(weightgain), length(columns), dimnames=list(rownames(weightgain), columns))
for(column in columns){
  histo.model <- lm(as.numeric(weightgain[, column]) ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data = weightgain)
  histo.model2 <- lm(as.numeric(weightgain[, column]) ~ Sex + as.numeric(WG) + Mother, data = weightgain)
  adjusted = histo.model2$residuals + mean(as.numeric(weightgain[, column]), na.rm = T)
  WGadjusted[names(adjusted),column] = adjusted 
  #Test for differences in S1 vs S2
  anova.pvals.absolute <- c(anova.pvals.absolute, anova(histo.model)[[5]][4])
  # Test if we can believe the anova test !
  shapiro.pvals.absolute <- c(shapiro.pvals.absolute, shapiro.test(histo.model$residuals)$"p.value")
}

names(anova.pvals.absolute) <- columns
names(shapiro.pvals.absolute) <- columns
round(anova.pvals.absolute,4)

# Results W3 (Week 6) - W6 (Week 9) , W10 (Week 13), W11 (Week 14), W13 (Week 16), W14 (17) -  growth differs sigificantly different from each other  
# point W7 (Week 10) - W8 (Week 11) - trend visible

# bind genotype data

WGadjusted <- as.data.frame(WGadjusted)
WGadjusted <- cbind(WGadjusted, Genotype = NA)
WGadjusted[rownames(Genotype),"Genotype"] <- as.character(unlist(Genotype[,1]))

# plot weekly weight gain

indCC <- WGadjusted[which(WGadjusted[,"Genotype"] == "CC"),1:16]
indCT <- WGadjusted[which(WGadjusted[,"Genotype"] == "CT"),1:16]
indTT <- WGadjusted[which(WGadjusted[,"Genotype"] == "TT"),1:16]

indCC <- apply(indCC[,-16],2,as.numeric)
indCT <- apply(indCT[,-16],2,as.numeric)
indTT <- apply(indTT[,-16],2,as.numeric)

plot(x= c(1, 15), y = c(-4, 10), t = 'n')
apply(indCC[, -16], 1,function(x){ points(x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), y = as.numeric(x), col="green", t = 'l') })
apply(indCT[, -16], 1,function(x){ points(x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), y = as.numeric(x), col="blue", t = 'l') })
apply(indTT[, -16], 1,function(x){ points(x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), y = as.numeric(x), col="orange", t = 'l') })

statsCC <- apply(indCC, 2, function(x){ c(mean(x, na.rm = TRUE), sqrt(sd(x, na.rm = TRUE)), sd(x, na.rm = TRUE))})
statsCT <- apply(indCT, 2, function(x){ c(mean(x), sqrt(sd(x)), sd(x)) })
statsTT <- apply(indTT, 2, function(x){ c(mean(x), sqrt(sd(x)), sd(x)) }) 


plot(x= c(0,16),y= c(-3,10), ylab = "weight gain (g)", xlab= "time(weeks)", t= 'n', xaxt="n")
axis(1, at = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), c(28,35,42,49,56,63,70,77,84,91,98,105,112,119,126) / 7)

polygon(x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1), c(statsCC[1,] + statsCC[2,], (statsCC[1,] - statsCC[2,])[15:1]), col = rgb(255, 165, 0, 125, maxColorValue=255), border=NA)
points(x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), statsCC[1,], t = 'l', col="orange")

polygon(x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1), c(statsCT[1,] + statsCT[2,], (statsCT[1,] - statsCT[2,])[15:1]), col = rgb(0, 165, 255, 125, maxColorValue=255), border=NA)
points(x =c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), statsCT[1,], t = 'l', col="white")

polygon(x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1), c(statsTT[1,] + statsTT[2,], (statsTT[1,] - statsTT[2,])[15:1]), col = rgb(0, 255, 0, 125, maxColorValue=255), border=NA)
points(x =c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), statsTT[1,], t = 'l', col="lightgreen")

legend("topleft", c("CC", "CT", "TT"), fill=c("orange", rgb(0, 165, 255, 125, maxColorValue=255), "lightgreen"), border=NA) # weight measurements look spiky 
text(3,5, paste0("*"), cex = 3)
text(4,4, paste0("*"), cex = 3)
text(5,3.8, paste0("*"), cex = 3)
text(6,3.6, paste0("*"), cex = 3)
text(7,3.5, paste0("*"), cex = 3)
text(8,3.3, paste0("°"), cex = 3)
text(10,3.3, paste0("*"), cex = 3)
text(11,3.3, paste0("*"), cex = 3)
text(13,3.3, paste0("*"), cex = 3)
text(14,3.5, paste0("*"), cex = 3)

#END