# analyzing fat percentage M. longissimus
# (c) Danny Arends and Aimee Freiberg (HU Berlin), 2018 - 2024

#import data 

setwd("/Users/aimeefreiberg/Documents/Universtiy Bachelor/AG Brockmann/Bachelorarbeit/R_Data")
mdata <- read.table("Musclefat_B6S1.txt", sep="\t", header=TRUE, na.strings = c("","NA"))
Genotype <- read.table("Genotype1.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))
Factors <- read.table("factors_AIL.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))

# calculate mean of dry weight muscle and fat percentage per mouse 
averages <- matrix(NA, length(unique(mdata[,1])), 3, dimnames=list(unique(mdata[,1]), c("MuscleWet", "MuscleDry", "MuscleFat")))
for(mouse in unique(mdata[,1])) {
  # Calculate the mean values per mouse for Dry, Wet and Fat
  averages[mouse, ] <- apply(mdata[which(mdata[,1] == mouse),2:4], 2, mean, na.rm=TRUE)
}

# bind data in one table 
averages <- cbind (averages, Genotype = as.character(unlist(Genotype[,1])), Mother = as.character(unlist(Factors[,1])), WG = as.character(unlist(Factors[,2])), Sex = as.character(unlist(Factors[,3])))

columns <- c("MuscleWet","MuscleDry","MuscleFat")

#adjust muscle fat percentage because weight <1g

averages <- cbind(averages, "MuscleFat.adj" = -((-0.45492)/0.94828)+(as.numeric(averages[, "MuscleFat"])/0.94828)-(0.19858/0.94828*(as.numeric(averages[, "MuscleDry"]))))

# include dominance effect

genotypes_num = as.numeric(factor(averages[,"Genotype"], levels=c("CC", "CT", "TT")))
names(genotypes_num) = rownames(averages)
genotypes_dom = as.numeric(factor(averages[,"Genotype"], levels=c("CC", "CT", "TT")))
genotypes_dom[genotypes_dom == 3] <- 2  

averages <- as.data.frame.matrix(averages)

# Test for significance, and normality

columns <- c("MuscleWet","MuscleDry","MuscleFat.adj")
anova.pvals <- c()
shapiro.pvals <- c()
for(column in columns){
  protein.model <- lm(as.numeric(averages[, column]) ~ Mother + genotypes_dom, data = averages)
  #test for differences in Genotypes
  anova.pvals <- c(anova.pvals, anova(protein.model)[[5]][1])
  #test of we can believe the anova test
  shapiro.pvals <- c(shapiro.pvals, shapiro.test(protein.model$residuals)$"p.value")
}

names(anova.pvals) <- columns
names(shapiro.pvals) <- columns

anova.pvals  ## Results: MuscleFat is significantly different 

# relative muscle weight 
mdata <- read.table("Organs.txt", sep="\t", row.names=1, header=TRUE, na.strings = c("","NA"))

averages <- cbind(averages, "leanweight" = as.numeric(mdata[, "Schlachtgewicht"] - mdata[, "WATgon"] - mdata[, "WATsc"]))
averages <- cbind(averages, "MuscleWet.rel" = as.numeric(as.character(averages[, "MuscleWet"]))/as.numeric(as.character(averages[, "leanweight"]))*100)

columns <- c("MuscleWet.rel")
anova.pvals <- c()
shapiro.pvals <- c()
for(column in columns){
  fat.model <- lm(as.numeric(averages[, column]) ~ Sex + as.numeric(WG) + Mother + genotypes_dom, data = averages)
  #test for differences in Genotypes
  anova.pvals <- c(anova.pvals, anova(fat.model)[[5]][1])
  #test of we can believe the anova test
  shapiro.pvals <- c(shapiro.pvals, shapiro.test(fat.model$residuals)$"p.value")
}

# plot graph for muscle fat and relative weight 


par(mfrow=c(1,3))
par (mar= c(7,5,7,5))

# muscle fat
plot(x = c(0, 4), y=c(0, 47), xaxt="n", main="fat percentage", ylab="Muscle Fat [%]" , xlab="Genotype", t='n', las=2, cex.lab=1.7, cex.axis=1.7, cex.main=1.7)
boxplot(as.numeric(as.character(averages[, "MuscleFat.adj"])) ~ averages[, "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')
lines(x=c(1,2.5), y = c(44, 44))
lines(x=c(1,1), y = c(42, 44))
lines(x=c(2.5,2.5), y = c(44, 43))
lines(x=c(2,3), y = c(43, 43))
lines(x=c(3,3), y = c(42, 43))
lines(x=c(2,2), y= c(42, 43))
text(1.5, 46, paste0("***"), cex = 3)

#muscleweight

plot(x = c(0, 4), y=c(0, 2.8), xaxt="n", main="relative weight", ylab="Muscle Weight [%]" , xlab="Genotype", t='n', las=2, cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
boxplot(as.numeric(as.character(averages[, "MuscleWet.rel"])) ~ averages[, "Genotype"], add=TRUE, col=c(rgb(255, 165, 0, 125, maxColorValue=255), rgb(0, 165, 255, 125, maxColorValue=255), rgb(0, 255, 0, 125, maxColorValue=255)), yaxt='n')
lines(x=c(1,2.5), y = c(2.7, 2.7))
lines(x=c(1,1), y = c(2.5, 2.7))
lines(x=c(2.5,2.5), y = c(2.7, 2.6))
lines(x=c(2,3), y = c(2.6, 2.6))
lines(x=c(3,3), y = c(2.5, 2.6))
lines(x=c(2,2), y= c(2.5, 2.6))
text(1.5, 2.8, paste0("***"), cex = 3)  

#END