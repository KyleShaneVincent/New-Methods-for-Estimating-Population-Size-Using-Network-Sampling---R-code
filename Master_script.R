#The "Master_script.R" script can be used to run a simulation study using the methods presented in the paper
#"New Methods for Estimating Population Size Using Network Link-Tracing Sampling"

##################################################################################################################################################################
#Load the required packages and set the working directory
##################################################################################################################################################################

#Load the necessary packages (ensure that the pacman package is installed for easy loading)
pacman::p_load(boot, coda, dplyr, igraph, insight, mice, stringr, useful)

#Set the working directory to where the files are saved
setwd("C:/Users/kyles/Desktop/New Methods for Estimating Population Size Using Network Sampling/R Code")


##################################################################################################################################################################
#Load and plot the P90 project data sets (see the corresponding R script for further information)
##################################################################################################################################################################

source("Load_P90_data.R")
CSpop = 5492 #The empirical population (set to either 595 or 5492 as this corresponds to the population size of desired population to use)
num_nominations = 7 #The maximum number of nominations per respondent (set to Inf to keep with original network graph); this induces a directed graph if finite
plot.P90.data = FALSE #Set to TRUE to view a plot of the empirical population
Load_P90_data(CSpop, num_nominations, plot.P90.data)


##################################################################################################################################################################
#Select a sample with the option to plot it
##################################################################################################################################################################

set.seed(4) #Set a seed to recreate the output
source("Sample_selection.R")

alpha=c(0.015,0.0175,0.03,0.0325) #Initial sample selection parameters, corresponds to strata  #Initial sample selection parameters, corresponds to strata
beta_strata=matrix(0,4,4); #The link-tracing parameters, strata to strata
beta_strata[1,1]=0; beta_strata[1,2]=0; beta_strata[1,3]=0; beta_strata[1,4]=1
beta_strata[2,1]=0; beta_strata[2,2]=0; beta_strata[2,3]=0; beta_strata[2,4]=1
beta_strata[3,1]=0; beta_strata[3,2]=0; beta_strata[3,3]=0; beta_strata[3,4]=1 
beta_strata[4,1]=0; beta_strata[4,2]=0; beta_strata[4,3]=0; beta_strata[4,4]=1
beta_sex_worker=1 #The link-tracing parameters, acts on indicator for sex worker (set to NA if using CSpop=595)
beta_net_size=-1 #The link-tracing parameters, acts on out-degree of respondent
num_waves=2 #The number of waves
use_beta_est=TRUE #Must be TRUE in the simulation studies to use the beta parameter estimates in the MCMC procedure
return_beta_est=TRUE #To view the GLM output

#Select the sample
Sample_selection(alpha=alpha, beta_strata=beta_strata, beta_sex_worker=beta_sex_worker,	beta_net_size=beta_net_size, use_beta_est=use_beta_est, return_beta_est=return_beta_est)
for(k in 0:num_waves)
	cat("Size of wave ", k, " is ", length(s[[k+1]]), '\n')
cat("Total sample size is ", length(unlist(s)), '\n')

#Plot the sample
source("Plot_sample.R")
plot.sample=TRUE
Plot_sample(CSpop, plot.sample)

#View the GLM output
if(return_beta_est==TRUE)
	{
	cat("The strata parameter values:", '\n')
	print(beta_strata)	
	if(!is.na(beta_sex_worker))
		{
		cat("The sex worker parameter value:", '\n')
		print(beta_sex_worker)
		}
	cat("The log net size parameter value:", '\n')
	print(beta_net_size)
	cat('\n')
	cat("The estimates for the parameters (if not significant then is set to zero):", '\n')
	print(beta_estimates_list)
	cat('\n', '\n')
	}
	
if(return_beta_est==TRUE)
	{		
	net_size_example=c(1,3,ifelse(num_nominations!=Inf,num_nominations,25))

	if(is.na(beta_sex_worker))
		{
		for(k in 1:length(net_size_example))
			{
			cat("The true selection probabilities for with recruiter network size of", net_size_example[k],'\n')
			print(inv.logit(beta_strata + beta_net_size*log(net_size_example[k])))
			cat('\n')
			cat("The estimated selection probabilities for with recruiter network size of", net_size_example[k],'\n') 
			print(inv.logit(beta_estimates_list$Strata_Parameter_Estimates + beta_estimates_list$Log_Net_Size_Parameter_Estimate*log(net_size_example[k])))	
			cat('\n', '\n')
			}
		}
		
	if(!is.na(beta_sex_worker))
		{
		for(k in 1:length(net_size_example))
			{
			cat("The true selection probabilities for recruiting non-sex workers with recruiter network size of", net_size_example[k],'\n')
			print(inv.logit(beta_strata + beta_sex_worker*0 + beta_net_size*log(net_size_example[k])))
			cat("The estimated selection probabilities for recruiting non-sex workers with recruiter network size of", net_size_example[k],'\n') 
			print(inv.logit(beta_estimates_list$Strata_Parameter_Estimates + beta_estimates_list$Beta_Sex_Worker_Parameter_Estimate*0 +		
			beta_estimates_list$Log_Net_Size_Parameter_Estimate*log(net_size_example[k])))	
			cat("The true selection probabilities for recruiting sex workers with recruiter network size of", net_size_example[k],'\n')
			print(inv.logit(beta_strata + beta_sex_worker*1 + beta_net_size*log(net_size_example[k])))		
			cat('\n')
			cat("The estimated selection probabilities for recruiting sex workers with recruiter network size of", net_size_example[k],'\n') 
			print(inv.logit(beta_estimates_list$Strata_Parameter_Estimates + beta_estimates_list$Beta_Sex_Worker_Parameter_Estimate*1 +		
			beta_estimates_list$Log_Net_Size_Parameter_Estimate*log(net_size_example[k])))			
			cat('\n', '\n')
			}
		}
	}
	

##################################################################################################################################################################
#Calculate the Frank and Snijders (1994) estimate for each stratum based on original approach and the Frank and Snijders style estimate for each 
#stratum based on new approach, then present the point estimates and the standard errors for sample previously selected
##################################################################################################################################################################

var_est = "Seber" #Choose the variance type of estimator (must be "jack" or "Seber")
var_est_new = "Seber" #Choose the variance type of estimator (must be "jack" or "Seber")

#The FS estimators applied to each stratum
source("FS_estimate.R")
FS_estimate(s0=s$Wave_0, var_est=var_est)
cat("Preliminary estimates based on original estimators: linkage and node estimators with estimates for standard errors respectively", '\n')
cat("Point estimates and standard errors for estimators based on linkage counts:", '\n')
N.hat.strata; sqrt(N.hat.strata.var)
cat("Point estimates and standard errors for estimators based on node counts:", '\n')
N.hat.strata.1; sqrt(N.hat.strata.1.var)

#The new estimators applied to each stratum
source("FS_new_estimate.R")
FS_new_estimate(s0=s$Wave_0, var_est_new=var_est_new)
cat("Preliminary estimates based on new estimators: linkage and node estimators with estimates for standard errors respectively", '\n')
cat("Point estimates and standard errors for estimators based on linkage counts:", '\n')
N.hat.strata.new; sqrt(N.hat.strata.new.var)
cat("Point estimates and standard errors for estimators based on node counts:", '\n')
N.hat.strata.1.new; sqrt(N.hat.strata.1.new.var)


##################################################################################################################################################################
#Approximate the Rao-Blackwellized estimates and plot the MCMC chain 
#Then apply a Gelman-Rubin test to determine if convergence has been met based on seeded reorderings that are "overdispersed" and "underdispersed"
##################################################################################################################################################################

#Apply the MCMC procedure to the selected sample
use_beta_est = TRUE #Set to TRUE to use the estimates from the GLM fit
num_mcmc = 500 #The number of MCMC iterations per chain
num_chains = 1 #The number of MCMC chains
search_iterations = 10 #The number of search iterations for seeding the chains within the primary chain
var_est = "Seber" #Choose the variance type of estimator (must be "jack" or "Seber")
var_est_new = "Seber" #Choose the variance type of estimator (must be "jack" or "Seber")
source("MCMC.R")
source("Resample_selection.R")
MCMC(s=s, num_mcmc=num_mcmc, num_chains=num_chains, search_iterations=search_iterations, use_beta_est=use_beta_est, var_est=var_est, var_est_new=var_est_new)

cat("The acceptance rate of the MCMC chain:", '\n')
mean(accept_mcmc)

cat("The size of the strata of the study population", '\n')
table(pop$strata)


#Preliminary and improved estimates based on original estimators
FS_estimate(s0=s$Wave_0, var_est=var_est)
cat("Estimator based on linkage counts, preliminary and improved respectively", '\n')
N.hat.strata
apply(N.hat.strata.mcmc,2,mean)

cat("Estimator based on node counts, preliminary and improved respectively", '\n')
N.hat.strata.1
apply(N.hat.strata.1.mcmc,2,mean)


#Preliminary and improved estimates based on new estimators
FS_new_estimate(s0=s$Wave_0, var_est_new=var_est_new)
cat("New estimator based on linkage counts, preliminary and improved respectively", '\n')
N.hat.strata.new
apply(N.hat.strata.new.mcmc,2,mean)

cat("New estimator based on node counts, preliminary and improved respectively", '\n')
N.hat.strata.1.new
apply(N.hat.strata.1.new.mcmc,2,mean)


#Visual illustrations of the output
#Original estimator based on linkage counts
X11()
par(mfrow=c(num_strata,1))
for(k in 1:num_strata)
	{
	plot(N.hat.strata.mcmc[,k], typ='l', main=paste("MCMC Chain Trace Plot for Strata", k, "Original Population Size Estimator Based on Linkage Counts"), 
		xlab="MCMC Iteration", ylab="Strata Pop Size Est")
	abline(h=mean(N.hat.strata.mcmc[,k]), col="blue")
	}

#Original estimator based on node counts
X11()
par(mfrow=c(num_strata,1))
for(k in 1:num_strata)
	{
	plot(N.hat.strata.1.mcmc[,k], typ='l', main=paste("MCMC Chain Trace Plot for Strata", k, "Original Population Size Estimator Based on Node Counts"), 
		xlab="MCMC Iteration", ylab="Strata Pop Size Est")
	abline(h=mean(N.hat.strata.1.mcmc[,k]), col="blue")
	}

#New estimator based on linkage counts
X11()
par(mfrow=c(num_strata,1))
for(k in 1:num_strata)
	{
	plot(N.hat.strata.new.mcmc[,k], typ='l', main=paste("MCMC Chain Trace Plot for Strata", k, "New Population Size Estimator Based on Linkage Counts"), 
		xlab="MCMC Iteration", ylab="Strata Pop Size Est")
	abline(h=mean(N.hat.strata.new.mcmc[,k]), col="blue")
	}
		
#New estimator based on node counts
X11()
par(mfrow=c(num_strata,1))
for(k in 1:num_strata)
	{
	plot(N.hat.strata.1.new.mcmc[,k], typ='l', main=paste("MCMC Chain Trace Plot for Strata", k, "New Population Size Estimator Based on Node Counts"), 
		xlab="MCMC Iteration", ylab="Strata Pop Size Est")
	abline(h=mean(N.hat.strata.1.new.mcmc[,k]), col="blue")
	}		
		
	
#Search for two sample reorderings which will be used to seed the two primary chains 
#Note: The selection probabilities are larger than those obtained with the original sample ordering.
source("Search_reorderings.R")
search_iterations = 10 #The number of search iterations to seed chain
Search_reorderings(s=s, search_iterations=search_iterations, use_beta_est=use_beta_est, var_est=var_est, var_est_new=var_est_new)
s.reordering_1 = s.reordering
Search_reorderings(s=s, search_iterations=search_iterations, use_beta_est=use_beta_est, var_est=var_est, var_est_new=var_est_new)
s.reordering_2 = s.reordering


#Apply the Gelman-Rubin diagnostics to the two chains seeded through the search for overdispersed reorderings
source("Gelman_Rubin.R")
num_mcmc = 500 #The number of MCMC iterations per chain
num_chains = 1 #The number of chains to run within primary chains
search_iterations = 10 #The number of search iterations to seed chains within primary chains
Gelman_Rubin(s.reordering_1=s.reordering_1, s.reordering_2=s.reordering_2, num_mcmc=num_mcmc, num_chains=num_chains, search_iterations=search_iterations,
				use_beta_est=use_beta_est, var_est=var_est, var_est_new=var_est_new)


##################################################################################################################################################################
#Simulation study
##################################################################################################################################################################

SEEDZ = c(
1015,
1016,
1081, 
1138, 
1139, 
1227, 
1228, 
1229, 
1296, 
1346, 
1413, 
1414, 
1562, 
1567, 
1712, 
1713, 
1772, 
1773, 
1774, 
1837, 
1838, 
1897)

#Source the required files
source("Sample_selection.R")
source("Plot_sample.R")
source("FS_estimate.R")
source("FS_new_estimate.R")
source("MCMC.R")
source("Resample_selection.R")
source("Search_reorderings.R")
source("Gelman_Rubin.R")
source("Simulation_study.R")

#The simulation study parameters
num_sims = 2
alpha=c(0.015,0.0175,0.03,0.0325)#Initial sample selection parameters, corresponds to strata  #Initial sample selection parameters, corresponds to strata
beta_strata=matrix(0,4,4); #The link-tracing parameters, strata to strata
beta_strata[1,1]=0; beta_strata[1,2]=0; beta_strata[1,3]=0; beta_strata[1,4]=1
beta_strata[2,1]=0; beta_strata[2,2]=0; beta_strata[2,3]=0; beta_strata[2,4]=1
beta_strata[3,1]=0; beta_strata[3,2]=0; beta_strata[3,3]=0; beta_strata[3,4]=1 
beta_strata[4,1]=0; beta_strata[4,2]=0; beta_strata[4,3]=0; beta_strata[4,4]=1
beta_sex_worker=1 #The link-tracing parameters, acts on indicator for sex worker (set to NA if using CSpop=595)
beta_net_size=-1 #The link-tracing parameters, acts on out-degree of respondent
num_waves=1 #The number of waves
use_beta_est = TRUE #Set to TRUE to use the estimates from the GLM fit
return_beta_est = FALSE #To view the GLM output
var_est = "Seber" #Choose the variance type of estimator (must be "jack" or "Seber")
var_est_new = "Seber" #Choose the variance type of estimator (must be "jack" or "Seber")
num_mcmc = 100000 #The number of MCMC iterations per chain
num_chains = 3 #The number of MCMC chains
search_iterations = 1000 #The number of search iterations for seeding the second and third chain
seeds =  SEEDZ[23:24] #Enter as "Time" or 1:num_sims
plots = FALSE #Set to TRUE to see a plot of the sample for each simulation run
write.data = TRUE #Set to TRUE to write the data to file

#Run the simulation
Simulation_study(CSpop=CSpop, num_sims=num_sims, use_beta_est=use_beta_est, return_beta_est=return_beta_est, var_est=var_est, var_est_new=var_est_new,
alpha=alpha, beta_strata=beta_strata, beta_sex_worker=beta_sex_worker, beta_net_size=beta_net_size,
num_strata=num_strata, num_waves=num_waves, plots=plots, num_mcmc=num_mcmc, num_chains=num_chains, search_iterations=search_iterations,
write.data=write.data, seeds=seeds)


##################################################################################################################################################################
#View the simulation output
##################################################################################################################################################################

#First, load the necessary packages
pacman::p_load(igraph, stringr, mice)

#Set the working directory to where the files are saved
setwd("C:/Users/kyles/Desktop/New Methods for Estimating Population Size Using Network Sampling/R Code")

#Load and plot the empirical population
source("Load_P90_data.R")
CSpop = 5492 #The empirical population (set to either 595 or 5492 as this corresponds to the population size of desired population to use)
num_nominations = 7 #The maximum number of nominations per respondent (set to Inf to keep with original network graph); this induces a directed graph if finite
plot.P90.data = FALSE #Set to TRUE to view a plot of the empirical population
Load_P90_data(CSpop, num_nominations, plot.P90.data)


#View the results
options(scipen = 999)
source("Simulation_output.R")
read.data=TRUE
#filename="CSpop_5492_Sim_parameters_4_Seber_Seber_0.015_0.0175_0.03_0.0325_TRUE_1_0_0_0_0_1_0_0_0_0_1_0_0_0_0_1_TRUE_-1_1_3_1e+05.csv"
filename="CSpop_5492_Sim_parameters_4_Seber_Seber_0.015_0.0175_0.03_0.0325_TRUE_1_0_0_0_0_1_0_0_0_0_1_0_0_0_0_1_TRUE_-1_2_3_5e+05.csv"
#filename="CSpop_5492_Sim_parameters_4_Seber_Seber_0.015_0.0175_0.03_0.0325_TRUE_0_0_0_0_0_0_0_0_0_0_0_0_1_1_1_1_TRUE_-1_1_3_1e+05.csv"
filename="CSpop_5492_Sim_parameters_4_Seber_Seber_0.015_0.0175_0.03_0.0325_TRUE_0_0_0_0_0_0_0_0_0_0_0_0_1_1_1_1_TRUE_-1_2_3_5e+05.csv"

Simulation_output(read.data, filename)

cat("The simulation study parameters",'\n')
sampling_output


cat("Simulation results based on original estimators",'\n','\n')

cat("Estimator based on linkage counts",'\n')
simulation_output
if(num_mcmc > 0)
	simulation_output.RB
	
cat("Estimator based on node counts",'\n')
simulation_output.1
if(num_mcmc > 0)
	simulation_output.1.RB	
	

cat("Simulation results based on new estimators",'\n','\n')

cat("Estimator based on linkage counts",'\n')	
simulation_output_new
if(num_mcmc > 0)
	simulation_output_new.RB

cat("Estimator based on node counts",'\n')
simulation_output.1_new
if(num_mcmc > 0)
	simulation_output.1_new.RB	

