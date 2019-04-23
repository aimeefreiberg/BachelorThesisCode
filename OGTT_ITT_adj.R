# analyzing OGTT and ITT Test
# (c) Danny Arends and Aimee Freiberg (HU Berlin), 2018 - 2024

setwd("/Users/aimeefreiberg/Documents/Universtiy Bachelor/AG Brockmann/Bachelorarbeit/R_Data")
OGTT <- read.table("OGTT.txt",sep='\t',na.strings=c("", "NA", "-"), header=TRUE, row.names=1,colClasses="character")
Genotype <- read.table("Genotype1.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Factors <- read.table("factors_AIL.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
colnames(OGTT) <- c("BW","Glucose","T0","T15","T30","T60","T120")
OGTT[OGTT == "600+"] <- "600"

#load needed library

library(DescTools)

#calculate area under the curve (AUC) for every mouse -> 30 double intervall, 60 four time interval

aucs <- c()
for(x in 1:nrow(OGTT)){
  aucs <- c(aucs, AUC(x = c(0, 15, 30, 60,120), as.numeric(OGTT[x, 3:7])))
}

OGTT <- cbind(OGTT, auc = aucs)

#calculate adjusted auc 
# auc - base line auc (area from data at point 0 )

aucs.adj <- cbind(apply(OGTT[, 3:7],2, as.numeric) - mean(as.numeric(OGTT[, "T0"])))

aucs.new <- c()
for(x in 1:nrow(aucs.adj)){
  aucs.new <- c(aucs.new, AUC(x = c(0, 15, 30, 60,120), as.numeric(aucs.adj[x, 1:5])))
}

aucs.adj <- cbind(aucs.adj, auc = aucs.new)

#bind Genotypes to Dataframe

aucs.adj <- cbind(aucs.adj, Genotype = NA)
rownames(aucs.adj) <- rownames(Genotype)
aucs.adj[rownames(Genotype),"Genotype"] <- as.character(unlist(Genotype[,1])) 

OGTT <- cbind(OGTT, Genotype = NA)
OGTT[rownames(Genotype),"Genotype"] <- as.character(unlist(Genotype[,1]))



anova(lm(as.numeric(aucs.adj[, "auc"]) ~ aucs.adj[, "Genotype"])) # adjusted data
#genotypes differ significantly

# plot curve adjusted for baseline

#create dataframe

indCC <- OGTT[which(OGTT[,"Genotype"] == "CC"),3:7]
indCT <- OGTT[which(OGTT[,"Genotype"] == "CT"),3:7]
indTT <- OGTT[which(OGTT[,"Genotype"] == "TT"),3:7]
colnames(indCC) <- c("T0","T15","T30","T60","T120")
colnames(indCT) <- c("T0","T15","T30","T60","T120")
colnames(indTT) <- c("T0","T15","T30","T60","T120")

rel.CC = apply(indCC,2,as.numeric) - mean(as.numeric(indCC[, "T0"]))
rel.CT = apply(indCT,2,as.numeric) - mean(as.numeric(indCT[, "T0"]))
rel.TT = apply(indTT,2,as.numeric) - mean(as.numeric(indTT[, "T0"]))



######## ITT #############
ITT <- read.table("ITT.txt",sep='\t',na.strings=c("", "NA", "-"), header=TRUE, row.names=1)
colnames(ITT) <- c("BW","Insulin","T0","T15","T30","T60")

#load needed library

library(DescTools)

#calculate area under the curve (AUC) for every mouse -> 30 double intervall, 60 four time interval

aucs <- c()
for(x in 1:nrow(ITT)){
  aucs <- c(aucs, AUC(x = c(0, 15, 30, 60), as.numeric(ITT[x, 3:6])))
}

ITT <- cbind(ITT, auc = aucs)


#calculate adjusted auc 
# auc - base line auc (area from data at point 0 )

ITTaucs.adj <- cbind(apply(ITT[, 3:6],2, as.numeric) - mean(as.numeric(ITT[, "T0"])))

aucs.new1 <- c()
for(x in 1:nrow(ITTaucs.adj)){
  aucs.new1 <- c(aucs.new1, AUC(x = c(0, 15, 30, 60), as.numeric(ITTaucs.adj[x, 1:4])))
}

ITTaucs.adj <- cbind(ITTaucs.adj, auc = aucs.new1)

#bind Genotypes to Dataframe
ITT <- cbind(ITT, Genotype = NA)
ITT[rownames(Genotype),"Genotype"] <- as.character(unlist(Genotype[,1]))

ITTaucs.adj <- cbind(ITTaucs.adj, Genotype = (as.character(Genotype[, "Genotype"])))

#test for significant difference

anova(lm(as.numeric(ITTaucs.adj[, "auc"]) ~ ITTaucs.adj[, "Genotype"]))

# adjust for baseline

indCC <- ITT[which(ITT[,"Genotype"] == "CC"),3:6]
indCT <- ITT[which(ITT[,"Genotype"] == "CT"),3:6]
indTT <- ITT[which(ITT[,"Genotype"] == "TT"),3:6]

ITTrel.CC = apply(indCC,2,as.numeric) - mean(as.numeric(indCC[, "T0"]))
ITTrel.CT = apply(indCT,2,as.numeric) - mean(as.numeric(indCT[, "T0"]))
ITTrel.TT = apply(indTT,2,as.numeric) - mean(as.numeric(indTT[, "T0"]))

colnames(indCC) <- c("T0","T15","T30","T60")
colnames(indCT) <- c("T0","T15","T30","T60")
colnames(indTT) <- c("T0","T15","T30","T60")

#split up for genotype

aindCC <- ITT[which(ITT[,"Genotype"] == "CC"),3:6]
aindCT <- ITT[which(ITT[,"Genotype"] == "CT"),3:6]
aindTT <- ITT[which(ITT[,"Genotype"] == "TT"),3:6]

# plot 
layout(matrix(c(1, 1, 1, 3, 2, 2, 2, 4), 2, 4, byrow = TRUE))
par(mar = c(4,6,3,1))

#OGTT

statsCC <- apply(rel.CC, 2, function(x){ c(mean(as.numeric(x)), sqrt(sd(x)), sd(x)) })
statsCT <- apply(rel.CT, 2, function(x){ c(mean(as.numeric(x)), sqrt(sd(x)), sd(x)) })
statsTT <- apply(rel.TT, 2, function(x){ c(mean(as.numeric(x)), sqrt(sd(x)), sd(x)) })

plot(x= c(0, 120), y = c(0, 400), t = 'n', ylab="Blood Glucose* [mg/dl] ", xlab="Time [min]", cex.axis=1.7,  cex.lab=1.7, cex.main=1.7)

polygon(x = c(0,15,30,60,120, 120,60, 30, 15, 0), c(statsCC[1,] + statsCC[2,], (statsCC[1,] - statsCC[2,])[5:1]), col = rgb(255, 165, 0, 125, maxColorValue=255), border=NA)
points(x = c(0,15,30,60,120), statsCC[1,], t = 'l', col="orange")

polygon(x = c(0,15,30,60,120, 120,60, 30, 15, 0), c(statsCT[1,] + statsCT[2,], (statsCT[1,] - statsCT[2,])[5:1]), col = rgb(0, 165, 255, 125, maxColorValue=255), border=NA)
points(x = c(0,15,30,60,120), statsCT[1,], t = 'l', col="white")

polygon(x = c(0,15,30,60,120, 120,60, 30, 15, 0), c(statsTT[1,] + statsTT[2,], (statsTT[1,] - statsTT[2,])[5:1]), col = rgb(0, 255, 0, 125, maxColorValue=255), border=NA)
points(x = c(0,15,30,60,120), statsTT[1,], t = 'l', col="lightgreen")

legend("topright", c("CC", "CT", "TT"), fill=c("orange", rgb(0, 165, 255, 125, maxColorValue=255), "lightgreen"), border=NA)

# ITT

#plot

statsCC <- apply(ITTrel.CC, 2, function(x){ c(mean(x), sqrt(sd(x)), sd(x)) })
statsCT <- apply(ITTrel.CT, 2, function(x){ c(mean(x), sqrt(sd(x)), sd(x)) })
statsTT <- apply(ITTrel.TT, 2, function(x){ c(mean(x), sqrt(sd(x)), sd(x)) })

plot(x= c(0, 60), y = c(10,-80), t = 'n', ylab="Blood Glucose* [mg/dl]",xlab = "Time(minutes)", cex.axis=1.7,  cex.lab=1.7, cex.main=1.7)

polygon(x = c(0,15,30,60, 60, 30, 15, 0), c(statsCC[1,] + statsCC[2,], (statsCC[1,] - statsCC[2,])[4:1]), col = rgb(255, 165, 0, 125, maxColorValue=255), border=NA)
points(x = c(0,15,30,60), statsCC[1,], t = 'l', col="orange")

polygon(x = c(0,15,30,60, 60, 30, 15, 0), c(statsCT[1,] + statsCT[2,], (statsCT[1,] - statsCT[2,])[4:1]), col = rgb(0, 165, 255, 125, maxColorValue=255), border=NA)
points(x = c(0,15,30,60), statsCT[1,], t = 'l', col="white")

polygon(x = c(0,15,30,60, 60, 30, 15, 0), c(statsTT[1,] + statsTT[2,], (statsTT[1,] - statsTT[2,])[4:1]), col = rgb(0, 255, 0, 125, maxColorValue=255), border=NA)
points(x = c(0,15,30,60), statsTT[1,], t = 'l', col="lightgreen")

legend("topright", c("CC", "CT", "TT"), fill=c("orange", rgb(0, 165, 255, 125, maxColorValue=255), "lightgreen"), border=NA)

# AUC OGTT

plot(x = c(0, 4), y=c(0, 60000), xaxt="n", ylab="AUC", xlab="Genotype", t='n', las=2, cex.lab=2, cex.axis=1.5, cex.main=2, yaxt="n")
boxplot(as.numeric(as.character(aucs.adj[, "auc"])) ~ aucs.adj[, "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')
lines(x=c(1,2.5), y = c(57000, 57000))
lines(x=c(1,1), y = c(53000., 57000))
lines(x=c(2.5,2.5), y = c(57000, 55000))
lines(x=c(2,3), y = c(55000, 55000))
lines(x=c(3,3), y = c(53000, 55000))
lines(x=c(2,2), y= c(53000, 55000))
text(1.7, 59000, paste0("***"), cex = 3)


#AUC ITT
plot(x = c(0, 4), y=c(-10000, 10000), xaxt="n", ylab="AUC", xlab="Genotype", t='n', las=2, cex.lab=2, cex.axis=1.5, cex.main=2, yaxt="n")
boxplot(-as.numeric(as.character(ITTaucs.adj[, "auc"])) ~ ITTaucs.adj[, "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')
lines(x=c(1,2.5), y = c(9000, 9000))
lines(x=c(1,1), y = c(8000, 9000))
lines(x=c(2.5,2.5), y = c(9000, 8500))
lines(x=c(2,3), y = c(8500, 8500))
lines(x=c(3,3), y = c(8000, 8500))
lines(x=c(2,2), y= c(8000, 8500))
text(1.7, 9500, paste0("***"), cex = 3)

# END