
# Code to perform the simulations in Revisiting and expanding the meta-analysis of variation: The log coefficient of variation ratio, lnCVR by Alistair M. Senior, Wolfgang Viechtbauer, and Shinichi Nakagawa, submitted to Research Synthesis Methods.

# Code written by AM Senior @ The University of Sydney

# Set the WD
#setwd("/Users/alistairsenior/Dropbox (Sydney Uni)/Work-Home Share/lnCVR_simulation")


# Function to estimate lnCVR and S2

# Takes the mean, sd and n of each group
# Takes the type of lnCVR to calculate: lnCVR_1, or lnCVR_2 - note if the correlation is specified as 0 you get the independent case, but if it is non 0, you get the dependent versions (given in the main text as lnCVR_3 and lnCVR_4) 

# Returns yi and vi - calculated effects size and sampling variance

calc_lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN, EffectSize="lnCVR_1", verbose=F, N=EN, cor=NA){
	
	# For independent data
	if(cor == 0){
		# Estimator 1
		if(EffectSize == "lnCVR_1"){
			# Calcualte the log CVR
			yi<-log((ESD / EMean) / (CSD / CMean))
			vi<-CSD^2 / (CN * CMean^2) + 1 / (2 * (CN - 1)) + 
				ESD^2 / (EN * EMean^2) + 1 / (2 * (EN - 1)) 
		}
	
		# Estimator 2
		if(EffectSize == "lnCVR_2"){
			# Calcualte the log CVR
			yi<-log((ESD / EMean) / (CSD / CMean)) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1))) 
			yi<-yi + 0.5 * (((CSD^2) / (CN * CMean^2)) - ((ESD^2) / (EN * EMean^2)))	
			vi<-CSD^2 / (CN * CMean^2) + CSD^4 / (2 * CN^2 * CMean^4) +
				1 / (2 * (CN - 1)) + 1 / (2 * (CN - 1)^2) + 
				ESD^2 / (EN * EMean^2) + ESD^4 / (2 * EN^2 * EMean^4) +
				1 / (2 * (EN - 1)) + 1 / (2 * (EN - 1)^2)
		}		
	}

	# For dependent data
	if(cor != 0){
		# Estimator dep0
		if(EffectSize == "lnCVR_1"){
			# Calcualte the log CVR
			yi<-log((ESD / EMean) / (CSD / CMean))
			vi<-CSD^2 / (N * CMean^2) + ESD^2 / (N * EMean^2) - 
				cor * ((2 * CSD * ESD) / (N * CMean * EMean)) +
				1 / (N - 1) - cor^2 * (1 / (N - 1))
		}
	
		# Estimator dep1
		if(EffectSize == "lnCVR_2"){
			# Calcualte the log CVR
			yi<-log((ESD / EMean) / (CSD / CMean)) + 0.5 * (CSD^2 / (N * CMean^2) - ESD^2 / (N * EMean^2)) 	
			vi<-CSD^2 / (N * CMean^2) + ESD^2 / (N * EMean^2) - 
				cor * ((2 * CSD * ESD) / (N * CMean * EMean)) +
				(CSD^4 / (2 * N^2 * CMean^4)) + (ESD^4 / (2 * N^2 * EMean^4)) + 
				cor^2 * ((CSD^2 * ESD^2 * (CMean^4 + EMean^4)) / (2 * N^2 * CMean^4 * EMean^4))
			vi<-vi + 1 / (N - 1) - cor^2 * (1 / (N - 1)) + 1 / (N - 1)^2 + cor^4 * ((CSD^8 + ESD^8) / (2 * (N - 1)^2 * CSD^4 * ESD^4))	
		}	
	}
	
	# Return as a data frame
	return(as.data.frame(cbind(yi, vi)))
	
}

# Function to simulate data with specified parameters and calculate a sample effect size

# Takes the mean, sd and n of each group
# Takes the type of lnCVR to calculate: lnCVR_1, lnCVR_2. By Default all are assumed
# Takes the correlation (cor) for dependant data: by default 0 is assumed.
# Verbose states whether warnings are printed in the case that sampling variance B and C are specified, and cor != 0: default is FALSE

# Returns a dataframe [1, n_effects + n_SV] of the estimated effects and sampling variances

sim_effect<-function(CMean, CSD, CN, EMean, ESD, EN, cor=0, EffectSize = c("lnCVR_1", "lnCVR_2"), verbose=F){
	
	# Mass library
	library(MASS)
	
	if((cor != 0) & (CN != EN) & verbose){
		print("Warning - sample sizes are unbalanced and covariance != 0, simulation could be meaningless!")
	}
	
	# Convert the correlation to a covariance
	cov<-cor*ESD*CSD
	
	# Create a covariance matrix
	vcv<-array(cov, c(2, 2))
	diag(vcv)<-c(CSD^2, ESD^2)
	
	# Sample the two groups from multivariate normal distribution - note that: 1) I use the max of the two Ns, but in the case that the two are unbalanced the first 1:N samples in each group are used below, 2) sampled in a while statement to avoid negative sample means
	sample<-cbind(-1, -1)
	while(sum(apply(sample, 2, mean) < 0) > 0){
		sample<-mvrnorm(max(c(CN, EN)), c(CMean, EMean), vcv)
	}
	
	# Calculate the mean and SD of the samples
	dat_NA<-array(NA, c(max(c(CN, EN)) ,2))
	dat_NA[c(1:CN),1]<-sample[c(1:CN),1]
	dat_NA[c(1:EN),2]<-sample[c(1:EN),2]
	C_sample_mean<-mean(dat_NA[,1], na.rm=T)
	E_sample_mean<-mean(dat_NA[,2], na.rm=T)
	C_sample_SD<-sd(dat_NA[,1], na.rm=T)
	E_sample_SD<-sd(dat_NA[,2], na.rm=T)	
	
	# Calculate the correlation between sample values for inclusion in the effect estimation, if necessary, else assume it is 0
	if(cor != 0){
		cor_sample<-cor(dat_NA, use = "pairwise.complete.obs")[2]
	}else{
		cor_sample<-0	
	}
	
	# Create data frame to hold the results
	results<-as.data.frame(array(NA, c(1, length(EffectSize)*2)))
	names(results)<-c(EffectSize, paste0("S2_", EffectSize))
	
	# For each effect size estimator add in the estimated effect
	for(i in 1:length(EffectSize)){
		results[1, i]<-calc_lnCVR(CMean = C_sample_mean, CSD = C_sample_SD, CN = CN, EMean = E_sample_mean, ESD = E_sample_SD, EN = EN, EffectSize = EffectSize[i], verbose=verbose, cor=cor_sample)$yi
		results[1, i+length(EffectSize)]<-calc_lnCVR(CMean = C_sample_mean, CSD = C_sample_SD, CN = CN, EMean = E_sample_mean, ESD = E_sample_SD, EN = EN, EffectSize = EffectSize[i], verbose=verbose, cor=cor_sample)$vi
	}
	
	
	# Return the results	
	return(results)

}


# Function to simulate k effect sizes with specified parameters

# Takes the number of effect sizes to calculate: k
# Takes the mean, sd and n of each group
# Takes the type of lnCVR to calculate: lnCVR_1, lnCVR_2. By Default all are assumed
# Takes the type of estimator for variance lnCVR to calculate: A, A_adjusted, B or C. By default all are assumed
# Takes the correlation (cor) for dependant data: by default 0 is assumed.
# Verbose states whether warnings are printed

# Returns a dataframe containing k estimated effects and sampling variances, the expected values, accuracy (bias), and precision (IQR) of each estimator specified, accuracy and coverage for each effect size * sampling variance combination. 

sim_k_effects<-function(k = 1000, CMean, CSD, CN, EMean, ESD, EN, cor=0, EffectSize = c("lnCVR_1", "lnCVR_2")){
	
	# Simulated the effect sizes
	simulated_effects<-as.data.frame(t(replicate(k, sim_effect(CMean=CMean, CSD=CSD, CN=CN, EMean=EMean, ESD=ESD, EN=EN, cor=cor), simplify = "vector")))
	
	# Note up the sampling variances
	SVs<-paste0("S2_", EffectSize)
	
	# Calculate the expected effects
	expected_effect<-as.data.frame(array(NA, c(1, 1 + length(EffectSize))))
	names(expected_effect)<-paste0("expected_", c("lnCVR", paste0("S2_", EffectSize)))
	# Calculate the expected population lnCVR
	expected_effect$expected_lnCVR<-log((ESD / EMean) / (CSD / CMean))
	# For each SV estimator add in the expected estimated SV based on the population statistics	
	for(i in 1:length(EffectSize)){
		expected_effect[1, i+1]<-calc_lnCVR(CMean=CMean, CSD=CSD, CN=CN, EMean=EMean, ESD=ESD, EN=EN, EffectSize=EffectSize[i], cor=cor)$vi
	}
	
	# The mean deviations for each effect size point estimate, and the IQR of deviations
	deviations_point_est<-as.data.frame(array(NA, c(1, length(EffectSize)*4)))
	names(deviations_point_est)<-c(paste0("accuracy_", EffectSize), paste0("IQR_deviation_", EffectSize), paste0("var_", EffectSize), paste0("mean_S2_", EffectSize))
	# For each effect size estimator add in the mean deviation between the exptected effect from the population stats and each of the k simulated effect sizes, and range of IQR, sd among effect sizes and mean S2 from simulated effects
	for(i in 1:length(EffectSize)){
		
		# Pull out the ith set of simulated data
		sim_effects_i<-unlist(simulated_effects[,i])
		sim_variances_i<-unlist(simulated_effects[,i+length(EffectSize)])
		
		# Ratio of simualted effects to the true effect (mean and IQR)
		deviations_point_est[1, i]<-mean(sim_effects_i - expected_effect[1,1])
		deviations_point_est[1, i+length(EffectSize)]<-IQR(sim_effects_i - expected_effect[1,1])
		
		# SD of effect sizes and mean sampling variance
		deviations_point_est[1, i+length(EffectSize)*2]<-sd(sim_effects_i)^2
		deviations_point_est[1, i+length(EffectSize)*3]<-mean(sim_variances_i)
	}
	
	# The deviation between the SD of the point estimates and each sample variance, and the cover age from confidence intervals
	combinations<-expand.grid(SVs, EffectSize)
	combinations<-paste(combinations$Var2, combinations$Var1, sep="x")
	deviations_samp_var<-as.data.frame(array(NA, c(1, length(combinations)*3)))
	names(deviations_samp_var)<-c(paste0("accuracy_", combinations), paste0("coverageT_", combinations), paste0("coverageZ_", combinations))
	# For each and effect size and SV combination get the accuracy of the SV and coverage
	# Counters
	acc_count<-1
	t_count<-1+length(combinations)
	z_count<-1+length(combinations)*2
	
	# For each effect size estimator
	for(i in 1:length(EffectSize)){
						
		# Pull out the jth effect sizes
		sim_effects_i<-unlist(simulated_effects[,i])	
		
		# For each sampling variance estimator		
		for(j in 1:length(EffectSize)){
			
			# Pull out the jth set of variances
			sim_variances_j<-unlist(simulated_effects[,j+length(EffectSize)])
			
			# Ratio of deviation of jth estimated sampling variance based on population stats and variance in sampled effects - note here the value based on the population stats is closest to the 'true' value
			deviations_samp_var[1,acc_count]<-(expected_effect[,j+1] - sd(sim_effects_i)^2) / sd(sim_effects_i)^2 * 100
			
			# Coverage of 95% CI (t value)
			CI_excludes<-which((abs(sim_effects_i - expected_effect$expected_lnCVR) - qt(0.975, EN+CN-2) * sqrt(sim_variances_j)) > 0)
			deviations_samp_var[1,t_count]<-1 - length(CI_excludes) / length(sim_effects_i)
			
			# Coverage of 95% CI (z value)
			CI_excludes<-which((abs(sim_effects_i - expected_effect$expected_lnCVR) - qnorm(0.975) * sqrt(sim_variances_j)) > 0)
			deviations_samp_var[1,z_count]<-1 - length(CI_excludes) / length(sim_effects_i)
			
			# Iterate the counters
			acc_count<-acc_count+1
			t_count<-t_count+1
			z_count<-z_count+1
		}
	}
		
	# Aggregate the results and return
	output<-cbind(expected_effect, deviations_point_est, deviations_samp_var)
	return(output)
	
}

# Function to run the simulation with specified parameters

# Returns a dataframe with results for all pairwise combinations of parameters values, based on k replicates
# print specifies whether a progress bar should be generated
# core specifes which core on the computer you are running on (useful when running in parrallel)

run_sim<-function(k = k, CMean, CSD, EMean, ESD, n, imbalance, cor, EffectSize = c("lnCVR_1", "lnCVR_2"), print=TRUE, core = 1){
	
	# Load the tcltk package for the progress bar
	library(tcltk)
	
	# Get the combination of parameters
	parameters<-expand.grid(EMean = EMean, ESD = ESD, n=n, imbalance = imbalance, cor = cor)
	parameters$CN<-parameters$n
	parameters$EN<-round(parameters$n * parameters$imbalance)
	
	# Do the first set to get a results table
	results<-cbind(parameters[1,], sim_k_effects(k=k, CMean=CMean, CSD=CSD, CN=parameters$CN[1], EMean=parameters$EMean[1], ESD=parameters$ESD[1], EN=parameters$EN[1], cor=parameters$cor[1]))
	
	# Now do the rest
	for(i in 2:dim(parameters)[1]){
		
		# Print Progress, if we are doing this
		if(print){
			if(i == 2){
				pb<-tkProgressBar(min = 0, max = dim(parameters)[1], title = paste0("Simulation Progress on Core ", core))
			}
			setTkProgressBar(pb, i)
		}
		
		# Do the simulation for the ith set of parameters
		results_i<-cbind(parameters[i,], sim_k_effects(k=k, CMean=CMean, CSD=CSD, CN=parameters$CN[i], EMean=parameters$EMean[i], ESD=parameters$ESD[i], EN=parameters$EN[i], cor=parameters$cor[i]))	
	
		# Bind in to the whole results
		results<-rbind(results, results_i)
	
	}
	close(pb)
	
	# Return the results
	return(results)
	
}

############################################################################################
############################################################################################

# OK lets run this thing

# Load the do Parrallel package - makes life easier
library(doSNOW)
library(MASS)
library(tcltk)

# Set the parameters for the simulation

# Number of replicate effect sizes to simulate
k<-100 * 1000

# The mean and sd in the control group
CMean<-100
CSD<-20

# The mean and SD in the experimental group - values to give nice log ratios between population means and sds
EMean<-CMean*exp(seq(-0.5, 0.5, 0.1))
ESD<-CSD*exp(seq(-0.5, 0.5, 0.1))

# n of the control group
n<-c(8, 16, 42)

############### Start with the independent case

# # # Correlation between treatment and control data - here is set to 0 as we are assuming independent case
cor<-c(0)

# The n of the treatment group is set as round(n_control * imbalance): i.e. 1 is no imbalance
imbalance<-c(1, 1.25)

# I am going to run the simulation across length(EMean) cores - each will handle one value
cl<-makeCluster(length(EMean), outfile="")
registerDoSNOW(cl)

results<-foreach(i = 1:length(EMean)) %dopar% {
	
	# Run the simulation on each of the cores
	output<-run_sim(k = k, CMean = CMean, CSD = CSD, EMean = EMean[i], ESD = ESD, n = n, imbalance = imbalance, cor = cor, core = i)
	return(output)
	
}

# Produces a list of results where each item (n items = n cores it is run over) in the list is a dataframe of results giving the parameters of simulation, expected and observed results.

# Save the results
save(results, file="Sim_res_ind.Rdata")

# Save the parameters
parameters<-list()
parameters$k<-k
parameters$CMean<-CMean
parameters$CSD<-CSD
parameters$EMean<-EMean
parameters$ESD<-ESD
parameters$n<-n
parameters$imbalance<-imbalance
parameters$cor<-cor
save(parameters, file="Sim_parameters_ind.Rdata")


############### Now the dependent case

# Correlation between treatment and control data - here is set to 0 as we are assuming independent case
cor<-0.8

# No imabalance
imbalance<-c(1)

# I am going to run the simulation across length(EMean) cores - each will handle one value
cl<-makeCluster(length(EMean), outfile="")
registerDoSNOW(cl)

results<-foreach(i = 1:length(EMean)) %dopar% {
	
	# Run the simulation on each of the cores
	output<-run_sim(k = k, CMean = CMean, CSD = CSD, EMean = EMean[i], ESD = ESD, n = n, imbalance = imbalance, cor = cor, core = i)
	return(output)
	
}

# Produces a list of results where each item (n items = n cores it is run over) in the list is a dataframe of results giving the parameters of simulation, expected and observed results.

# Save the results
save(results, file="Sim_res_dep.Rdata")


# Save the parameters
parameters<-list()
parameters$k<-k
parameters$CMean<-CMean
parameters$CSD<-CSD
parameters$EMean<-EMean
parameters$ESD<-ESD
parameters$n<-n
parameters$imbalance<-imbalance
parameters$cor<-cor
save(parameters, file="Sim_parameters_dep.Rdata")

