#The "Simulation_output" function will summarize the simulation study output.

Simulation_output <- function(read.data, filename)
	{
	if(read.data==TRUE)
		{
		sim_output = read.csv(filename,header=FALSE)
		num_strata = mean(sim_output[,2])
		colnames_sim_output = c("seeds", "num_strata")
		use_beta_est = sim_output[1,3]
		colnames_sim_output = c(colnames_sim_output, "use_beta_est")
		var_est = sim_output[1,4]
		colnames_sim_output = c(colnames_sim_output, "var_est")
		var_est_new = sim_output[1,5]
		colnames_sim_output = c(colnames_sim_output, "var_est_new")
		num_waves = mean(sim_output[1,(3+num_strata+num_strata*num_strata+1+1+1+2)])
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("alpha",k,sep="_"))
		for(k in 1:num_strata)
			for(j in 1:num_strata)
				colnames_sim_output = c(colnames_sim_output, paste("beta_strata",k,j,sep="_"))		
		colnames_sim_output = c(colnames_sim_output,"beta_sex_worker","beta_net_size","num_waves")
		for(k in 1:(num_waves+1))
			colnames_sim_output = c(colnames_sim_output, paste("n.sim",k,sep="_"))
		for(k in 1:(num_waves+1))
			colnames_sim_output = c(colnames_sim_output, paste("n.sex.worker.sim",k,sep="_"))
		colnames_sim_output = c(colnames_sim_output,"accept_mcmc.sim","num_mcmc","num_chains","search_iterations")
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.sim",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.var.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.1.sim",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.1.var.sim",k,sep="_"))				
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.strata.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.strata.log",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.strata.1.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.strata.1.log",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.new.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.new.var.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.1.new.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.strata.1.new.var.sim",k,sep="_"))			
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.strata.new.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.strata.new.log",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.strata.1.new.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.strata.1.new.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.sim",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.var.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.var.sim.negative",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.1.sim",k,sep="_"))	
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.1.var.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.1.var.sim.negative",k,sep="_"))			
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.RB.strata.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.RB.strata.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.RB.strata.1.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.RB.strata.1.log",k,sep="_"))			
		for(k in 1:num_strata)		
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.new.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.new.var.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.new.var.sim.negative",k,sep="_"))
		for(k in 1:num_strata)		
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.1.new.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.1.new.var.sim",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("N.hat.RB.strata.1.new.var.sim.negative",k,sep="_"))			
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.RB.strata.new.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.RB.strata.new.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("S.RB.strata.1.new.log",k,sep="_"))
		for(k in 1:num_strata)
			colnames_sim_output = c(colnames_sim_output, paste("C.RB.strata.1.new.log",k,sep="_"))			
			
		colnames(sim_output) = colnames_sim_output
		
		seeds <<- sim_output[,"seeds"]
		num_strata = mean(sim_output[,"num_strata"])
		x = which(str_detect(colnames(sim_output),"use_beta_est"))
		use_beta_est = mean(as.numeric(sim_output[,"use_beta_est"]))
		x = which(str_detect(colnames(sim_output),"var_est"))
		var_est = sim_output[1,"var_est"]
		x = which(str_detect(colnames(sim_output),"var_est_new"))
		var_est_new = sim_output[1,"var_est_new"]
		x = which(str_detect(colnames(sim_output),"alpha"))
		alpha = apply(as.matrix(sim_output[,x]),2,mean); names(alpha)=NULL
		x = which(str_detect(colnames(sim_output),"beta_strata"))
		beta_strata = apply(as.matrix(sim_output[,x]),2,mean); names(beta_strata)=NULL
		x = which(str_detect(colnames(sim_output),"beta_sex_worker"))
		beta_sex_worker = mean(as.numeric(sim_output[,"beta_sex_worker"]))
		x = which(str_detect(colnames(sim_output),"beta_net_size"))
		beta_net_size = mean(as.numeric(sim_output[,"beta_net_size"]))
		num_waves = mean(sim_output[,"num_waves"])
		x = which(str_detect(colnames(sim_output),"n.sim"))
		n.sim = sim_output[,x]; names(n.sim)=NULL
		x = which(str_detect(colnames(sim_output),"n.sex.worker.sim"))
		n.sex.worker.sim = sim_output[,x]; names(n.sex.worker.sim)=NULL
		num_sims = nrow(sim_output)
		num_mcmc <<- mean(sim_output[,"num_mcmc"])
		num_chains <<- mean(sim_output[,"num_chains"])
		search_iterations <<- mean(sim_output[,"search_iterations"])
		accept_mcmc.sim = sim_output[,"accept_mcmc.sim"]		
		x = which(str_detect(colnames(sim_output),"N.hat.strata.sim"))
		N.hat.strata.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.strata.var.sim"))
		N.hat.strata.var.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.strata.1.sim"))
		N.hat.strata.1.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.strata.1.var.sim"))
		N.hat.strata.1.var.sim = sim_output[,x]				
		x = which(str_detect(colnames(sim_output),"S.strata.log"))
		S.strata.log = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"C.strata.log"))
		C.strata.log = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"S.strata.1.log"))
		S.strata.1.log = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"C.strata.1.log"))
		C.strata.1.log = sim_output[,x]		
		x = which(str_detect(colnames(sim_output),"N.hat.strata.new.sim"))
		N.hat.strata.new.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.strata.new.var.sim"))
		N.hat.strata.new.var.sim = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"N.hat.strata.1.new.sim"))
		N.hat.strata.1.new.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.strata.1.new.var.sim"))
		N.hat.strata.1.new.var.sim = sim_output[,x]			
		x = which(str_detect(colnames(sim_output),"S.strata.new.log"))
		S.strata.new.log = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"C.strata.new.log"))
		C.strata.new.log = sim_output[,x]		
		x = which(str_detect(colnames(sim_output),"S.strata.1.new.log"))
		S.strata.1.new.log = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"C.strata.1.new.log"))
		C.strata.1.new.log = sim_output[,x]			
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.sim"))
		N.hat.RB.strata.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.var.sim"))
		N.hat.RB.strata.var.sim = sim_output[,x[1:num_strata]]	
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.var.sim.negative"))
		N.hat.RB.strata.var.sim.negative = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.1.sim"))
		N.hat.RB.strata.1.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.1.var.sim"))
		N.hat.RB.strata.1.var.sim = sim_output[,x[1:num_strata]]	
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.1.var.sim.negative"))
		N.hat.RB.strata.1.var.sim.negative = sim_output[,x]		
		x = which(str_detect(colnames(sim_output),"S.RB.strata.log"))
		S.RB.strata.log = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"C.RB.strata.log"))
		C.RB.strata.log = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"S.RB.strata.1.log"))
		S.RB.strata.1.log = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"C.RB.strata.1.log"))
		C.RB.strata.1.log = sim_output[,x]				
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.new.sim"))
		N.hat.RB.strata.new.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.new.var.sim"))
		N.hat.RB.strata.new.var.sim = sim_output[,x[1:num_strata]]	
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.new.var.sim.negative"))
		N.hat.RB.strata.new.var.sim.negative = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.1.new.sim"))
		N.hat.RB.strata.1.new.sim = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.1.new.var.sim"))
		N.hat.RB.strata.1.new.var.sim = sim_output[,x[1:num_strata]]	
		x = which(str_detect(colnames(sim_output),"N.hat.RB.strata.1.new.var.sim.negative"))
		N.hat.RB.strata.1.new.var.sim.negative = sim_output[,x]					
		x = which(str_detect(colnames(sim_output),"S.RB.strata.new.log"))
		S.RB.strata.new.log = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"C.RB.strata.new.log"))
		C.RB.strata.new.log = sim_output[,x]
		x = which(str_detect(colnames(sim_output),"S.RB.strata.1.new.log"))
		S.RB.strata.1.new.log = sim_output[,x]	
		x = which(str_detect(colnames(sim_output),"C.RB.strata.1.new.log"))
		C.RB.strata.1.new.log = sim_output[,x]		
		sim_output <<- sim_output		
		}

	#Some calculations required for displaying the output
	strata_size = table(pop$strata)

	#The sampling output
	sampling_output = list()
	sampling_output[[1]] = num_strata
	names(sampling_output)[[1]] = "Number of strata"	
	sampling_output[[2]] = alpha
	names(sampling_output)[[2]] = "Bernoulli sampling parameters for strata"
	sampling_output[[3]] = c(beta_strata, beta_sex_worker, beta_net_size)
	names(sampling_output)[[3]] = "Beta logit sampling parameters for tracing links, beta strata, beta_sex_worker, and beta net size"
	sampling_output[[4]] = num_waves
	names(sampling_output)[[4]] = "Number of waves"
	sampling_output[[5]] = apply(n.sim,2,mean)
	names(sampling_output)[[5]] = "Average sample size by wave"	
	sampling_output[[6]] = mean(apply(n.sim,1,sum))
	names(sampling_output)[[6]] = "Average total sample size"		
	sampling_output[[7]] = apply(n.sex.worker.sim,2,mean)
	names(sampling_output)[[7]] = "Average number of sex workers recruited by wave"		
	sampling_output[[8]] = mean(apply(n.sex.worker.sim,1,sum))
	names(sampling_output)[[8]] = "Average total number of sex workers recruited"
	sampling_output[[9]] = ifelse(use_beta_est==TRUE, "Estimated Values Used", "True Values Used")
	names(sampling_output)[[9]] = "Beta parameters"
	sampling_output[[10]] = num_sims
	names(sampling_output)[[10]] = "Number of simulation runs"
	sampling_output[[11]] = var_est
	names(sampling_output)[[11]] = "Variance estimation procedure for original estimators"
	sampling_output[[12]] = var_est_new
	names(sampling_output)[[12]] = "Variance estimation procedure for new estimators"
	k = 12
	if(read.data==TRUE)
		{
		sampling_output[[k+read.data]] = length(unique(seeds))
		names(sampling_output)[[k+read.data]] = "Number of unique seeds"
		}
	if(num_mcmc > 0)
		{
		sampling_output[[k+1+read.data]] = mean(accept_mcmc.sim)
		names(sampling_output)[[k+1+read.data]] = "Mean acceptance rate of MCMC routine"
		sampling_output[[k+2+read.data]] = num_mcmc
		names(sampling_output)[[k+2+read.data]] = "Length of MCMC chain"
		sampling_output[[k+3+read.data]] = num_chains
		names(sampling_output)[[k+3+read.data]] = "Number of MCMC chains per sample"
		sampling_output[[k+4+read.data]] = search_iterations
		names(sampling_output)[[k+4+read.data]] = "Number of search iterations for finding overdispersed reorderings"		
		}

	#Calculations for coverage rates
	N.hat.strata.CI.LB.sim = N.hat.strata.sim - 1.96 * sqrt(N.hat.strata.var.sim)
	N.hat.strata.CI.UB.sim = N.hat.strata.sim + 1.96 * sqrt(N.hat.strata.var.sim)
	N.hat.strata.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.sim)),ncol(as.matrix(N.hat.strata.sim)))
	for(k in 1:num_strata)
		N.hat.strata.CI.coverage[,k] = (as.matrix(N.hat.strata.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.CI.UB.sim)[,k])

	N.hat.strata.log.CI.LB.sim = S.strata.log + (N.hat.strata.sim-S.strata.log)/C.strata.log
	N.hat.strata.log.CI.UB.sim = S.strata.log + (N.hat.strata.sim-S.strata.log)*C.strata.log
	N.hat.strata.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.sim)),ncol(as.matrix(N.hat.strata.sim)))
	for(k in 1:num_strata)
		N.hat.strata.log.CI.coverage[,k] = (as.matrix(N.hat.strata.log.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.log.CI.UB.sim)[,k])
										
										
	N.hat.strata.1.CI.LB.sim = N.hat.strata.1.sim - 1.96 * sqrt(N.hat.strata.1.var.sim)
	N.hat.strata.1.CI.UB.sim = N.hat.strata.1.sim + 1.96 * sqrt(N.hat.strata.1.var.sim)
	N.hat.strata.1.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.1.sim)),ncol(as.matrix(N.hat.strata.1.sim)))
	for(k in 1:num_strata)
		N.hat.strata.1.CI.coverage[,k] = (as.matrix(N.hat.strata.1.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.1.CI.UB.sim)[,k])

	N.hat.strata.1.log.CI.LB.sim = S.strata.1.log + (N.hat.strata.1.sim-S.strata.1.log)/C.strata.1.log
	N.hat.strata.1.log.CI.UB.sim = S.strata.1.log + (N.hat.strata.1.sim-S.strata.1.log)*C.strata.1.log
	N.hat.strata.1.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.1.sim)),ncol(as.matrix(N.hat.strata.1.sim)))
	for(k in 1:num_strata)
		N.hat.strata.1.log.CI.coverage[,k] = (as.matrix(N.hat.strata.1.log.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.1.log.CI.UB.sim)[,k])										
										
										
										
	if(num_mcmc > 0)
		{
		N.hat.RB.strata.CI.LB.sim = N.hat.RB.strata.sim - 1.96 * sqrt(N.hat.RB.strata.var.sim)
		N.hat.RB.strata.CI.UB.sim = N.hat.RB.strata.sim + 1.96 * sqrt(N.hat.RB.strata.var.sim)
		N.hat.RB.strata.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.sim)),ncol(as.matrix(N.hat.RB.strata.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.CI.UB.sim)[,k])	

		N.hat.RB.strata.log.CI.LB.sim = S.RB.strata.log + (N.hat.RB.strata.sim-S.RB.strata.log)/C.RB.strata.log
		N.hat.RB.strata.log.CI.UB.sim = S.RB.strata.log + (N.hat.RB.strata.sim-S.RB.strata.log)*C.RB.strata.log
		N.hat.RB.strata.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.sim)),ncol(as.matrix(N.hat.RB.strata.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.log.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.log.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.log.CI.UB.sim)[,k])	


		N.hat.RB.strata.1.CI.LB.sim = N.hat.RB.strata.1.sim - 1.96 * sqrt(N.hat.RB.strata.1.var.sim)
		N.hat.RB.strata.1.CI.UB.sim = N.hat.RB.strata.1.sim + 1.96 * sqrt(N.hat.RB.strata.1.var.sim)
		N.hat.RB.strata.1.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.1.sim)),ncol(as.matrix(N.hat.RB.strata.1.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.1.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.1.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.1.CI.UB.sim)[,k])	

		N.hat.RB.strata.1.log.CI.LB.sim = S.RB.strata.1.log + (N.hat.RB.strata.1.sim-S.RB.strata.1.log)/C.RB.strata.1.log
		N.hat.RB.strata.1.log.CI.UB.sim = S.RB.strata.1.log + (N.hat.RB.strata.1.sim-S.RB.strata.1.log)*C.RB.strata.1.log
		N.hat.RB.strata.1.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.1.sim)),ncol(as.matrix(N.hat.RB.strata.1.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.1.log.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.1.log.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.1.log.CI.UB.sim)[,k])											
		}								

	N.hat.strata.new.CI.LB.sim = N.hat.strata.new.sim - 1.96 * sqrt(N.hat.strata.new.var.sim)
	N.hat.strata.new.CI.UB.sim = N.hat.strata.new.sim + 1.96 * sqrt(N.hat.strata.new.var.sim)
	N.hat.strata.new.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.new.sim)),ncol(as.matrix(N.hat.strata.new.sim)))
	for(k in 1:num_strata)
		N.hat.strata.new.CI.coverage[,k] = (as.matrix(N.hat.strata.new.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.new.CI.UB.sim)[,k])		

	N.hat.strata.new.log.CI.LB.sim = S.strata.new.log + (N.hat.strata.new.sim-S.strata.new.log)/C.strata.new.log
	N.hat.strata.new.log.CI.UB.sim = S.strata.new.log + (N.hat.strata.new.sim-S.strata.new.log)*C.strata.new.log
	N.hat.strata.new.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.new.sim)),ncol(as.matrix(N.hat.strata.new.sim)))
	for(k in 1:num_strata)
		N.hat.strata.new.log.CI.coverage[,k] = (as.matrix(N.hat.strata.new.log.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.new.log.CI.UB.sim)[,k])
										
	N.hat.strata.1.new.CI.LB.sim = N.hat.strata.1.new.sim - 1.96 * sqrt(N.hat.strata.1.new.var.sim)
	N.hat.strata.1.new.CI.UB.sim = N.hat.strata.1.new.sim + 1.96 * sqrt(N.hat.strata.1.new.var.sim)
	N.hat.strata.1.new.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.1.new.sim)),ncol(as.matrix(N.hat.strata.1.new.sim)))
	for(k in 1:num_strata)
		N.hat.strata.1.new.CI.coverage[,k] = (as.matrix(N.hat.strata.1.new.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.1.new.CI.UB.sim)[,k])		

	N.hat.strata.1.new.log.CI.LB.sim = S.strata.1.new.log + (N.hat.strata.1.new.sim-S.strata.1.new.log)/C.strata.1.new.log
	N.hat.strata.1.new.log.CI.UB.sim = S.strata.1.new.log + (N.hat.strata.1.new.sim-S.strata.1.new.log)*C.strata.1.new.log
	N.hat.strata.1.new.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.strata.1.new.sim)),ncol(as.matrix(N.hat.strata.1.new.sim)))
	for(k in 1:num_strata)
		N.hat.strata.1.new.log.CI.coverage[,k] = (as.matrix(N.hat.strata.1.new.log.CI.LB.sim)[,k] < strata_size[k] & 
										strata_size[k] < as.matrix(N.hat.strata.1.new.log.CI.UB.sim)[,k])										
										
										

	if(num_mcmc > 0)
		{
		N.hat.RB.strata.new.CI.LB.sim = N.hat.RB.strata.new.sim - 1.96 * sqrt(N.hat.RB.strata.new.var.sim)
		N.hat.RB.strata.new.CI.UB.sim = N.hat.RB.strata.new.sim + 1.96 * sqrt(N.hat.RB.strata.new.var.sim)
		N.hat.RB.strata.new.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.new.sim)),ncol(as.matrix(N.hat.RB.strata.new.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.new.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.new.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.new.CI.UB.sim)[,k])							
		
		N.hat.RB.strata.new.log.CI.LB.sim = S.RB.strata.new.log + (N.hat.RB.strata.new.sim-S.RB.strata.new.log)/C.RB.strata.new.log
		N.hat.RB.strata.new.log.CI.UB.sim = S.RB.strata.new.log + (N.hat.RB.strata.new.sim-S.RB.strata.new.log)*C.RB.strata.new.log
		N.hat.RB.strata.new.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.new.sim)),ncol(as.matrix(N.hat.RB.strata.new.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.new.log.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.new.log.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.new.log.CI.UB.sim)[,k])	


		N.hat.RB.strata.1.new.CI.LB.sim = N.hat.RB.strata.1.new.sim - 1.96 * sqrt(N.hat.RB.strata.1.new.var.sim)
		N.hat.RB.strata.1.new.CI.UB.sim = N.hat.RB.strata.1.new.sim + 1.96 * sqrt(N.hat.RB.strata.1.new.var.sim)
		N.hat.RB.strata.1.new.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.1.new.sim)),ncol(as.matrix(N.hat.RB.strata.1.new.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.1.new.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.1.new.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.1.new.CI.UB.sim)[,k])							
		
		N.hat.RB.strata.1.new.log.CI.LB.sim = S.RB.strata.1.new.log + (N.hat.RB.strata.1.new.sim-S.RB.strata.1.new.log)/C.RB.strata.1.new.log
		N.hat.RB.strata.1.new.log.CI.UB.sim = S.RB.strata.1.new.log + (N.hat.RB.strata.1.new.sim-S.RB.strata.1.new.log)*C.RB.strata.1.new.log
		N.hat.RB.strata.1.new.log.CI.coverage = matrix(NA,nrow(as.matrix(N.hat.RB.strata.1.new.sim)),ncol(as.matrix(N.hat.RB.strata.1.new.sim)))
		for(k in 1:num_strata)
			N.hat.RB.strata.1.new.log.CI.coverage[,k] = (as.matrix(N.hat.RB.strata.1.new.log.CI.LB.sim)[,k] < strata_size[k] & 
											strata_size[k] < as.matrix(N.hat.RB.strata.1.new.log.CI.UB.sim)[,k])											
		}									


	#Construct the simulation output for the original Frank and Snijders estimate based on linkage counts
	simulation_output = list()

	simulation_output[[1]] = data.frame(matrix(NA,num_strata,8))
	for(k in 1:num_strata)
		rownames(simulation_output[[1]])[k] = paste("Stratum", k, sep=" ")
	colnames(simulation_output[[1]])[1] = "Strata Sizes"
	simulation_output[[1]][,1] = strata_size
	colnames(simulation_output[[1]])[2] = "Mean"
	simulation_output[[1]][,2] = apply(as.matrix(N.hat.strata.sim),2,mean)
	colnames(simulation_output[[1]])[3] = "SE"
	simulation_output[[1]][,3] = apply(as.matrix(N.hat.strata.sim),2,sd)
	colnames(simulation_output[[1]])[4] = "MSE"
	simulation_output[[1]][,4] = (apply(as.matrix(N.hat.strata.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.strata.sim),2,var)
	colnames(simulation_output[[1]])[5] = "CI Coverage"
	simulation_output[[1]][,5] = apply(N.hat.strata.CI.coverage,2,mean)	
	colnames(simulation_output[[1]])[6] = "Avg. Half-Length CI"
	simulation_output[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.strata.var.sim))*1.96*2,2,mean)
	colnames(simulation_output[[1]])[7] = "Logtrans CI Coverage"
	simulation_output[[1]][,7] = apply(N.hat.strata.log.CI.coverage,2,mean)	
	colnames(simulation_output[[1]])[8] = "Avg. Half-Length Logtrans CI"
	simulation_output[[1]][,8] = 1/2*(apply(as.matrix(N.hat.strata.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.strata.log.CI.LB.sim),2,mean))

	names(simulation_output)[1] = "Performance of Original Frank and Snijders Estimator Per Stratum"

	if(num_mcmc > 0)
		{
		simulation_output.RB = list()

		simulation_output.RB[[1]] = data.frame(matrix(NA,num_strata,10))
		for(k in 1:num_strata)
			rownames(simulation_output.RB[[1]])[k] = paste("Stratum", k, sep=" ")
		colnames(simulation_output.RB[[1]])[1] = "Strata Sizes"
		simulation_output.RB[[1]][,1] = strata_size
		colnames(simulation_output.RB[[1]])[2] = "Mean"
		simulation_output.RB[[1]][,2] = apply(as.matrix(N.hat.RB.strata.sim),2,mean)
		colnames(simulation_output.RB[[1]])[3] = "SE"
		simulation_output.RB[[1]][,3] = apply(as.matrix(N.hat.RB.strata.sim),2,sd)
		colnames(simulation_output.RB[[1]])[4] = "MSE"
		simulation_output.RB[[1]][,4] = (apply(as.matrix(N.hat.RB.strata.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.RB.strata.sim),2,var)
		colnames(simulation_output.RB[[1]])[5] = "CI Coverage"
		simulation_output.RB[[1]][,5] = apply(N.hat.RB.strata.CI.coverage,2,mean)		
		colnames(simulation_output.RB[[1]])[6] = "Avg. Half-Length CI"
		simulation_output.RB[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.RB.strata.var.sim))*1.96*2,2,mean)
		colnames(simulation_output.RB[[1]])[7] = "Log Trans CI Coverage"
		simulation_output.RB[[1]][,7] = apply(N.hat.RB.strata.log.CI.coverage,2,mean)		
		colnames(simulation_output.RB[[1]])[8] = "Avg. Half-Length Logtrans CI"
		simulation_output.RB[[1]][,8] = 1/2*(apply(as.matrix(N.hat.RB.strata.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.RB.strata.log.CI.LB.sim),2,mean))
		colnames(simulation_output.RB[[1]])[9] = "Percent Neg Var Ests"
		simulation_output.RB[[1]][,9] = apply(as.matrix(N.hat.RB.strata.var.sim.negative),2,mean)	
		colnames(simulation_output.RB[[1]])[10] = "Ratio of Vars of RB and Base Ests"
		simulation_output.RB[[1]][,10] = apply(as.matrix(N.hat.RB.strata.sim),2,var)/apply(as.matrix(N.hat.strata.sim),2,var)
		names(simulation_output.RB)[1] = "Performance of Improved Original Frank and Snijders Estimator Per Stratum"
		}



	#Construct the simulation output for the original Frank and Snijders estimate based on node counts
	simulation_output.1 = list()

	simulation_output.1[[1]] = data.frame(matrix(NA,num_strata,8))
	for(k in 1:num_strata)
		rownames(simulation_output.1[[1]])[k] = paste("Stratum", k, sep=" ")
	colnames(simulation_output.1[[1]])[1] = "Strata Sizes"
	simulation_output.1[[1]][,1] = strata_size
	colnames(simulation_output.1[[1]])[2] = "Mean"
	simulation_output.1[[1]][,2] = apply(as.matrix(N.hat.strata.1.sim),2,mean)
	colnames(simulation_output.1[[1]])[3] = "SE"
	simulation_output.1[[1]][,3] = apply(as.matrix(N.hat.strata.1.sim),2,sd)
	colnames(simulation_output.1[[1]])[4] = "MSE"
	simulation_output.1[[1]][,4] = (apply(as.matrix(N.hat.strata.1.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.strata.1.sim),2,var)	
	colnames(simulation_output.1[[1]])[5] = "CI Coverage"
	simulation_output.1[[1]][,5] = apply(N.hat.strata.1.CI.coverage,2,mean)			
	colnames(simulation_output.1[[1]])[6] = "Avg. Half-Length CI"
	simulation_output.1[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.strata.1.var.sim))*1.96*2,2,mean)
	colnames(simulation_output.1[[1]])[7] = "Logtrans CI Coverage"
	simulation_output.1[[1]][,7] = apply(N.hat.strata.1.log.CI.coverage,2,mean)	
	colnames(simulation_output.1[[1]])[8] = "Avg. Half-Length Logtrans CI"
	simulation_output.1[[1]][,8] = 1/2*(apply(as.matrix(N.hat.strata.1.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.strata.1.log.CI.LB.sim),2,mean))
	names(simulation_output.1)[1] = "Performance of Original Frank and Snijders Estimator Per Stratum"

	if(num_mcmc > 0)
		{
		simulation_output.1.RB = list()

		simulation_output.1.RB[[1]] = data.frame(matrix(NA,num_strata,10))
		for(k in 1:num_strata)
			rownames(simulation_output.1.RB[[1]])[k] = paste("Stratum", k, sep=" ")
		colnames(simulation_output.1.RB[[1]])[1] = "Strata Sizes"
		simulation_output.1.RB[[1]][,1] = strata_size
		colnames(simulation_output.1.RB[[1]])[2] = "Mean"
		simulation_output.1.RB[[1]][,2] = apply(as.matrix(N.hat.RB.strata.1.sim),2,mean)
		colnames(simulation_output.1.RB[[1]])[3] = "SE"
		simulation_output.1.RB[[1]][,3] = apply(as.matrix(N.hat.RB.strata.1.sim),2,sd)
		colnames(simulation_output.1.RB[[1]])[4] = "MSE"
		simulation_output.1.RB[[1]][,4] = (apply(as.matrix(N.hat.RB.strata.1.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.RB.strata.1.sim),2,var)		
		colnames(simulation_output.1.RB[[1]])[5] = "CI Coverage"
		simulation_output.1.RB[[1]][,5] = apply(N.hat.RB.strata.1.CI.coverage,2,mean)				
		colnames(simulation_output.1.RB[[1]])[6] = "Avg. Half-Length CI"
		simulation_output.1.RB[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.RB.strata.1.var.sim))*1.96*2,2,mean)		
		colnames(simulation_output.1.RB[[1]])[7] = "Log Trans CI Coverage"
		simulation_output.1.RB[[1]][,7] = apply(N.hat.RB.strata.1.log.CI.coverage,2,mean)			
		colnames(simulation_output.1.RB[[1]])[8] = "Avg. Half-Length Logtrans CI"
		simulation_output.1.RB[[1]][,8] = 1/2*(apply(as.matrix(N.hat.RB.strata.1.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.RB.strata.1.log.CI.LB.sim),2,mean))		
		colnames(simulation_output.1.RB[[1]])[9] = "Percent Neg Var Ests"
		simulation_output.1.RB[[1]][,9] = apply(as.matrix(N.hat.RB.strata.1.var.sim.negative),2,mean)	
		colnames(simulation_output.1.RB[[1]])[10] = "Ratio of Vars of RB and Base Ests"
		simulation_output.1.RB[[1]][,10] = apply(as.matrix(N.hat.RB.strata.1.sim),2,var)/apply(as.matrix(N.hat.strata.1.sim),2,var)
		names(simulation_output.1.RB)[1] = "Performance of Improved Original Frank and Snijders Estimator Per Stratum"
		}


	#Construct the simulation output for the new estimate based on linkage counts
	simulation_output_new = list()

	simulation_output_new[[1]] = data.frame(matrix(NA,num_strata,8))
	for(k in 1:num_strata)
		rownames(simulation_output_new[[1]])[k] = paste("Stratum", k, sep=" ")
	colnames(simulation_output_new[[1]])[1] = "Strata Sizes"
	simulation_output_new[[1]][,1] = strata_size
	colnames(simulation_output_new[[1]])[2] = "Mean"
	simulation_output_new[[1]][,2] = apply(as.matrix(N.hat.strata.new.sim),2,mean)
	colnames(simulation_output_new[[1]])[3] = "SE"
	simulation_output_new[[1]][,3] =  apply(as.matrix(N.hat.strata.new.sim),2,sd)
	colnames(simulation_output_new[[1]])[4] = "MSE"
	simulation_output_new[[1]][,4] = (apply(as.matrix(N.hat.strata.new.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.strata.new.sim),2,var)	
	colnames(simulation_output_new[[1]])[5] = "CI Coverage"
	simulation_output_new[[1]][,5] = apply(N.hat.strata.new.CI.coverage,2,mean)		
	colnames(simulation_output_new[[1]])[6] = "Avg. Half-Length CI"
	simulation_output_new[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.strata.new.var.sim))*1.96*2,2,mean)	
	colnames(simulation_output_new[[1]])[7] = "Logtrans CI Coverage"
	simulation_output_new[[1]][,7] = apply(N.hat.strata.new.log.CI.coverage,2,mean)		
	colnames(simulation_output_new[[1]])[8] = "Avg. Half-Length Logtrans CI"
	simulation_output_new[[1]][,8] = 1/2*(apply(as.matrix(N.hat.strata.new.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.strata.new.log.CI.LB.sim),2,mean))	
	names(simulation_output_new)[1] = "Performance of New Estimator Per Stratum"
	
	if(num_mcmc > 0)
		{
		simulation_output_new.RB = list()

		simulation_output_new.RB[[1]] = data.frame(matrix(NA,num_strata,10))
		for(k in 1:num_strata)
			rownames(simulation_output_new.RB[[1]])[k] = paste("Stratum", k, sep=" ")
		colnames(simulation_output_new.RB[[1]])[1] = "Strata Sizes"
		simulation_output_new.RB[[1]][,1] = strata_size
		colnames(simulation_output_new.RB[[1]])[2] = "Mean"
		simulation_output_new.RB[[1]][,2] = apply(as.matrix(N.hat.RB.strata.new.sim),2,mean)
		colnames(simulation_output_new.RB[[1]])[3] = "SE"
		simulation_output_new.RB[[1]][,3] =  apply(as.matrix(N.hat.RB.strata.new.sim),2,sd)
		colnames(simulation_output_new.RB[[1]])[4] = "MSE"
		simulation_output_new.RB[[1]][,4] = (apply(as.matrix(N.hat.RB.strata.new.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.RB.strata.new.sim),2,var)	
		colnames(simulation_output_new.RB[[1]])[5] = "CI Coverage"
		simulation_output_new.RB[[1]][,5] = apply(N.hat.RB.strata.new.CI.coverage,2,mean)		
		colnames(simulation_output_new.RB[[1]])[6] = "Avg. Half-Length CI"
		simulation_output_new.RB[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.RB.strata.new.var.sim))*1.96*2,2,mean)	
		colnames(simulation_output_new.RB[[1]])[7] = "Logtrans CI Coverage"
		simulation_output_new.RB[[1]][,7] = apply(N.hat.RB.strata.new.log.CI.coverage,2,mean)			
		colnames(simulation_output_new.RB[[1]])[8] = "Avg. Half-Length Logtrans CI"
		simulation_output_new.RB[[1]][,8] = 1/2*(apply(as.matrix(N.hat.RB.strata.new.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.RB.strata.new.log.CI.LB.sim),2,mean))
			colnames(simulation_output_new.RB[[1]])[9] = "Percent Neg Var Ests"
		simulation_output_new.RB[[1]][,9] = apply(as.matrix(N.hat.RB.strata.new.var.sim.negative),2,mean)	
		colnames(simulation_output_new.RB[[1]])[10] = "Ratio of Vars of RB and Base Ests"
		simulation_output_new.RB[[1]][,10] = apply(as.matrix(N.hat.RB.strata.new.sim),2,var)/apply(as.matrix(N.hat.strata.new.sim),2,var)
		names(simulation_output_new.RB)[1] = "Performance of Improved New Estimator Per Stratum"
		}
	

	#Construct the simulation output for the new estimate based on linkage counts
	simulation_output.1_new = list()

	simulation_output.1_new[[1]] = data.frame(matrix(NA,num_strata,8))
	for(k in 1:num_strata)
		rownames(simulation_output.1_new[[1]])[k] = paste("Stratum", k, sep=" ")
	colnames(simulation_output.1_new[[1]])[1] = "Strata Sizes"
	simulation_output.1_new[[1]][,1] = strata_size
	colnames(simulation_output.1_new[[1]])[2] = "Mean"
	simulation_output.1_new[[1]][,2] = apply(as.matrix(N.hat.strata.1.new.sim),2,mean)
	colnames(simulation_output.1_new[[1]])[3] = "SE"
	simulation_output.1_new[[1]][,3] =  apply(as.matrix(N.hat.strata.1.new.sim),2,sd)
	colnames(simulation_output.1_new[[1]])[4] = "MSE"
	simulation_output.1_new[[1]][,4] = (apply(as.matrix(N.hat.strata.1.new.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.strata.1.new.sim),2,var)
	colnames(simulation_output.1_new[[1]])[5] = "CI Coverage"
	simulation_output.1_new[[1]][,5] = apply(N.hat.strata.1.new.CI.coverage,2,mean)		
	colnames(simulation_output.1_new[[1]])[6] = "Avg. Half-Length CI"
	simulation_output.1_new[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.strata.1.new.var.sim))*1.96*2,2,mean)
	colnames(simulation_output.1_new[[1]])[7] = "Logtrans CI Coverage"
	simulation_output.1_new[[1]][,7] = apply(N.hat.strata.1.new.log.CI.coverage,2,mean)	
	colnames(simulation_output.1_new[[1]])[8] = "Avg. Half-Length Logtrans CI"
	simulation_output.1_new[[1]][,8] = 1/2*(apply(as.matrix(N.hat.strata.1.new.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.strata.1.new.log.CI.LB.sim),2,mean))
	
	names(simulation_output.1_new)[1] = "Performance of New Estimator Per Stratum"
	
	if(num_mcmc > 0)
		{
		simulation_output.1_new.RB = list()

		simulation_output.1_new.RB[[1]] = data.frame(matrix(NA,num_strata,10))
		for(k in 1:num_strata)
			rownames(simulation_output.1_new.RB[[1]])[k] = paste("Stratum", k, sep=" ")
		colnames(simulation_output.1_new.RB[[1]])[1] = "Strata Sizes"
		simulation_output.1_new.RB[[1]][,1] = strata_size
		colnames(simulation_output.1_new.RB[[1]])[2] = "Mean"
		simulation_output.1_new.RB[[1]][,2] = apply(as.matrix(N.hat.RB.strata.1.new.sim),2,mean)
		colnames(simulation_output.1_new.RB[[1]])[3] = "SE"
		simulation_output.1_new.RB[[1]][,3] =  apply(as.matrix(N.hat.RB.strata.1.new.sim),2,sd)
		colnames(simulation_output.1_new.RB[[1]])[4] = "MSE"
		simulation_output.1_new.RB[[1]][,4] = (apply(as.matrix(N.hat.RB.strata.1.new.sim),2,mean)-strata_size)^2+apply(as.matrix(N.hat.RB.strata.1.new.sim),2,var)		
		colnames(simulation_output.1_new.RB[[1]])[5] = "CI Coverage"
		simulation_output.1_new.RB[[1]][,5] = apply(N.hat.RB.strata.1.new.CI.coverage,2,mean)		
		colnames(simulation_output.1_new.RB[[1]])[6] = "Avg. Half-Length CI"
		simulation_output.1_new.RB[[1]][,6] = 1/2*apply(sqrt(as.matrix(N.hat.RB.strata.1.new.var.sim))*1.96*2,2,mean)	
		colnames(simulation_output.1_new.RB[[1]])[7] = "Logtrans CI Coverage"
		simulation_output.1_new.RB[[1]][,7] = apply(N.hat.RB.strata.1.new.log.CI.coverage,2,mean)	
		colnames(simulation_output.1_new.RB[[1]])[8] = "Avg. Half-Length Logtrans CI"
		simulation_output.1_new.RB[[1]][,8] = 1/2*(apply(as.matrix(N.hat.RB.strata.1.new.log.CI.UB.sim),2,mean) - apply(as.matrix(N.hat.RB.strata.1.new.log.CI.LB.sim),2,mean))	
		colnames(simulation_output.1_new.RB[[1]])[9] = "Percent Neg Var Ests"
		simulation_output.1_new.RB[[1]][,9] = apply(as.matrix(N.hat.RB.strata.1.new.var.sim.negative),2,mean)	
		colnames(simulation_output.1_new.RB[[1]])[10] = "Ratio of Vars of RB and Base Ests"
		simulation_output.1_new.RB[[1]][,10] = apply(as.matrix(N.hat.RB.strata.1.new.sim),2,var)/apply(as.matrix(N.hat.strata.1.new.sim),2,var)
		names(simulation_output.1_new.RB)[1] = "Performance of Improved New Estimator Per Stratum"
		}
	
	#The output to be displayed
	sampling_output <<- sampling_output
	simulation_output <<- simulation_output
	simulation_output.1 <<- simulation_output.1
	simulation_output_new <<- simulation_output_new
	simulation_output.1_new <<- simulation_output.1_new
	if(num_mcmc > 0)
		{
		simulation_output.RB <<- simulation_output.RB
		simulation_output.1.RB <<- simulation_output.1.RB
		simulation_output_new.RB <<- simulation_output_new.RB
		simulation_output.1_new.RB <<- simulation_output.1_new.RB
		}
	}