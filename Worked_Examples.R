

# Code to perform the worked exmaples in Revisiting and expanding the meta-analysis of variation: The log coefficient of variation ratio, lnCVR by Alistair M. Senior, Wolfgang Viechtbauer, and Shinichi Nakagawa, submitted to Research Synthesis Methods.

# Written by AM Senior @ The University of Sydney

# load metafor package
library(metafor)

# Function to estimate lnCVR and S2

# Returns yi and vi - calculated effects size and sampling variance

calc_lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN, verbose=F, N=EN, cor=0){
	
	# For independent data
	if(cor == 0){
		yi<-log((ESD / EMean) / (CSD / CMean)) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1))) 
		yi<-yi + 0.5 * (((CSD^2) / (CN * CMean^2)) - ((ESD^2) / (EN * EMean^2)))	
		vi<-CSD^2 / (CN * CMean^2) + CSD^4 / (2 * CN^2 * CMean^4) +
			1 / (2 * (CN - 1)) + 1 / (2 * (CN - 1)^2) + 
			ESD^2 / (EN * EMean^2) + ESD^4 / (2 * EN^2 * EMean^4) +
			1 / (2 * (EN - 1)) + 1 / (2 * (EN - 1)^2)
	}

	# For dependent data
	if(cor != 0){	
		yi<-log((ESD / EMean) / (CSD / CMean)) + 0.5 * (CSD^2 / (N * CMean^2) - ESD^2 / (N * EMean^2))
		vi<-CSD^2 / (N * CMean^2) + ESD^2 / (N * EMean^2) - 
			cor * ((2 * CSD * ESD) / (N * CMean * EMean)) +
			1 / (N - 1) - cor^2 * (1 / (N - 1))
	}
		
	# Return as a data frame
	return(as.data.frame(cbind(yi, vi)))
	
}


pdf("Plots.pdf", height=5, width=10)
par(mfrow=c(1,2), cex.lab=1.35)

############ Ecology Example

data<-get(data(dat.curtis1998))
head(data)

# Note group 1 is under elevated CO2 - the treatment, group 2 is under ambient CO2


# We will randomise the order of the data, so we are sure to end up with a random subset for partially dependent analyses
set.seed(123)
data<-data[sample(seq(1, dim(data)[1], 1)),]

# Is there a mean variance relationship - yes

means<-c(data$m1i, data$m2i)
sds<-c(data$sd1i, data$sd2i)
treats<-c(rep(1, dim(data)[1]), rep(16, dim(data)[1]))

plot(log(means), log(sds), pch=treats, xlab=expression(ln~italic(bar(x))), ylab=expression(ln~italic(s)), cex=0.5)
mtext("A", side=2, at=6, line=1, las=2, cex=2)

legend(-2.25, 4.75, c("Treatment", "Control"), pch=c(1, 16))

# What proportions of the data set should we assume are correlated
prop_cor<-seq(0, 1, 0.2)

# Create a dataframe to hold the resuls
results<-data.frame(prop_cor = prop_cor, Est=0, SE=0, LCI=0, UCI=0, Tau2=0, I2=0, Q=0, Q_p=0)

# Add in the correlations
data$cor<-0

# Calculate the effect sizes assuming complete independence
ES<-data.frame(yi=rep(0, dim(data)[1]))
ES$vi<-0
for(i in 1:dim(data)[1]){
	ES[i,]<-calc_lnCVR(CMean = data$m2i[i], CSD = data$sd2i[i], CN = data$n2i[i], EMean = data$m1i[i], ESD = data$sd1i[i], EN = data$n1i[i], cor=data$cor[i])
}


# Fit the model
model<-rma(yi = ES$yi, vi = ES$vi)
results[1,-1]<-c(model$b, model$se, model$ci.lb, model$ci.ub, model$tau2, model$I2, model$QE, model$QEp)

#pdf("funnel_eco.pdf")

# Do the loop for the other proportions
for(c in 2:length(prop_cor)){
	
	# Add in the assumed correlation
	data$cor[c(1:round(dim(data)[1] * prop_cor[c]))]<-0.8
	
	# Recalculate the effect sizes
	ES<-data.frame(yi=rep(0, dim(data)[1]))
	ES$vi<-0
	for(i in 1:dim(data)[1]){
		ES[i,]<-calc_lnCVR(CMean = data$m2i[i], CSD = data$sd2i[i], CN = data$n2i[i], EMean = data$m1i[i], ESD = data$sd1i[i], EN = data$n1i[i], cor=data$cor[i])
	}


	# Fit the model
	model<-rma(yi = ES$yi, vi = ES$vi)
	results[c,-1]<-c(model$b, model$se, model$ci.lb, model$ci.ub, model$tau2, model$I2, model$QE, model$QEp)
	#plot(ES$yi, sqrt(ES$vi)*-1, col=as.factor(data$cor), main=prop_cor[c], ylim=c(-1.2,0))
}

#dev.off()

write.table(signif(results, 4), "res_eco.csv", sep=",", row.names=F, col.names=names(results))

#### Clinical Example

data<-read.csv("JBM_Meta.csv")
head(data)

# Note group 1 is under low GI, which will treat as the 'treatment' grouo, where as group 2 is udner high GI, which we will treat as the control group

# Assume that the n is total unique individuals, so for cross over that is n in each group, but for parrallel design n is halved

data$n1i<-data$n
data$n2i<-data$n
data$n1i[which(data$Design == "P")]<-data$n1i[which(data$Design == "P")] / 2
data$n2i[which(data$Design == "P")]<-data$n2i[which(data$Design == "P")] / 2

# Is there a mean variance relationship?  - yes

means<-c(data$m1i, data$m2i)
sds<-c(data$sd1i, data$sd2i)
units<-as.factor(c(data$Outcome, data$Outcome))
treats<-c(rep(1, dim(data)[1]), rep(16, dim(data)[1]))

plot(log(means), log(sds), col=as.factor(units), pch=treats, xlab=expression(ln~italic(bar(x))), ylab=expression(ln~italic(s)), cex=1)
mtext("B", side=2, at=0.9, line=1, las=2, cex=2)

legend(1.75, -0.75, c("Fructosamine", expression(HBA[1][c])), col=c(1, 2), pch=15)

# JBM says assume cor is 0.34 for carryovers - lets test analysis assuming 0.3, 0.5, 0.8
corrs<-c(0, 0.3, 0.5, 0.8)

# Create a dataframe to hold the resuls
results<-data.frame(assumed_corr = corrs, Est=0, SE=0, LCI=0, UCI=0, Tau2=0, I2=0, Q=0, Q_p=0)

for(c in 1:length(corrs)){

	# Add in the assumed correlation
	data$cor<-corrs[c]
	data$cor[which(data$Design == "P")]<-0
		
	# Calculate the effect sizes assuming complete independence
	ES<-data.frame(yi=rep(0, dim(data)[1]))
	ES$vi<-0
	for(i in 1:dim(data)[1]){
		ES[i,]<-calc_lnCVR(CMean = data$m2i[i], CSD = data$sd2i[i], CN = data$n2i[i], EMean = data$m1i[i], ESD = data$sd1i[i], EN = data$n1i[i], cor=data$cor[i])
	}
	
	# Fit the model
	model<-rma(yi = ES$yi, vi = ES$vi)
	results[c,-1]<-c(model$b, model$se, model$ci.lb, model$ci.ub, model$tau2, model$I2, model$QE, model$QEp)

}

write.table(signif(results, 4), "res_clin.csv", sep=",", row.names=F, col.names=names(results))

100 * (1 - exp(results$Est))

dev.off()
