#The "Simulation_study" function will run a simulation study.

Simulation_study <- function(CSpop, num_sims, use_beta_est, return_beta_est, var_est, var_est_new, alpha, beta_strata, beta_sex_worker, beta_net_size, num_strata, num_waves, plots, 
								num_mcmc, num_chains, search_iterations, write.data, seeds)
	{	
	if(length(alpha) != num_strata)
		{
		cat("Error: Length of alpha vector must be equal to num_strata '\n'")
		stop()
		}

	n.sim = matrix(NA,num_sims,num_waves+1)
	n.sex.worker.sim = matrix(NA,num_sims,num_waves+1)
	
	N.hat.strata.sim = matrix(NA,num_sims,num_strata)
	N.hat.strata.var.sim = matrix(NA,num_sims,num_strata)
	N.hat.strata.1.sim = matrix(NA,num_sims,num_strata)
	N.hat.strata.1.var.sim = matrix(NA,num_sims,num_strata)
	
	S.strata.log = matrix(NA,num_sims,num_strata)
	C.strata.log = matrix(NA,num_sims,num_strata)
	S.strata.1.log = matrix(NA,num_sims,num_strata)
	C.strata.1.log = matrix(NA,num_sims,num_strata)	
	
	N.hat.strata.new.sim = matrix(NA,num_sims,num_strata)
	N.hat.strata.new.var.sim = matrix(NA,num_sims,num_strata)
	N.hat.strata.1.new.sim = matrix(NA,num_sims,num_strata)
	N.hat.strata.1.new.var.sim = matrix(NA,num_sims,num_strata)
	
	S.strata.new.log = matrix(NA,num_sims,num_strata)
	C.strata.new.log = matrix(NA,num_sims,num_strata)	
	S.strata.1.new.log = matrix(NA,num_sims,num_strata)
	C.strata.1.new.log = matrix(NA,num_sims,num_strata)	
		
	N.hat.RB.strata.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.var.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.var.sim.negative = matrix(0,num_sims,num_strata)
	
	N.hat.RB.strata.1.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.1.var.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.1.var.sim.negative = matrix(0,num_sims,num_strata)		
	
	S.RB.strata.log = matrix(NA,num_sims,num_strata)
	C.RB.strata.log = matrix(NA,num_sims,num_strata)
	S.RB.strata.1.log = matrix(NA,num_sims,num_strata)
	C.RB.strata.1.log = matrix(NA,num_sims,num_strata)		
	
	N.hat.RB.strata.new.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.new.var.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.new.var.sim.negative = matrix(0,num_sims,num_strata)
	
	N.hat.RB.strata.1.new.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.1.new.var.sim = matrix(NA,num_sims,num_strata)
	N.hat.RB.strata.1.new.var.sim.negative = matrix(0,num_sims,num_strata)	
	
	S.RB.strata.new.log = matrix(NA,num_sims,num_strata)
	C.RB.strata.new.log = matrix(NA,num_sims,num_strata)	
	
	S.RB.strata.1.new.log = matrix(NA,num_sims,num_strata)
	C.RB.strata.1.new.log = matrix(NA,num_sims,num_strata)
	
	accept_mcmc.sim = numeric(num_sims)

	cat("Progress of simulation study", '\n')
	progress_sim = txtProgressBar(min = 0, max = num_sims, style = 3)
	for(kk in 1:num_sims)
		{
		#Set the seed
		seed = ifelse(seeds[1] == "Time",as.numeric(Sys.time()), seeds[kk])
		set.seed(seed)
		
		#Select the sample
		Sample_selection(alpha=alpha,beta_strata=beta_strata,beta_sex_worker,beta_net_size=beta_net_size,use_beta_est=use_beta_est,return_beta_est=return_beta_est)
		for(k in 1:(num_waves+1))
			{
			n.sim[kk,k] = length(s[[k]])
			n.sex.worker.sim[kk,k] = length(which(s[[k]] %in% pop.sex.worker[[1]]))
			}
		
		#Plot the sample		
		if(plots==TRUE)
			{	
			if(CSpop==5492)
				{
				s.plot = unlist(s); names(s.plot) = NULL
				n.plot = length(s.plot)
				vertex.color = numeric(n.plot)
				vertex.color[which(pop$strata[s.plot] %in% c(1,2))] = "blue" #White
				vertex.color[which(pop$strata[s.plot] %in% c(3,4))] = "red" #Non-white
				vertex.shape = numeric(n.plot)
				vertex.shape[which(pop$strata[s.plot] %in% c(1,3))]= "square" #Males
				vertex.shape[which(pop$strata[s.plot] %in% c(2,4))]= "circle" #Females
				vertex.size = rep(2,n.plot)
				vertex.size[which(s.plot %in% s$Wave_0)] = 3.5
				G = graph.adjacency(adjacency[c(s.plot),c(s.plot)])
				loc = layout_components(G)
				plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
				vertex.shape=vertex.shape,edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Sample Graph with Seeds Represented as Englarged Nodes",
				xlab="White Males as blue squares, White Females as blue circles, Non-White Males as red squares, Non-White Females as red circles")
				}
			if(CSpop==595)
				{
				#Plot the sample		
				s.plot = unlist(s); names(s.plot) = NULL
				n.plot = length(s.plot)
				vertex.label = s.plot
				vertex.color = numeric(n.plot)
				vertex.color[which(pop$strata[s.plot] == 1)] = "yellow" #non-injection drug user
				vertex.color[which(pop$strata[s.plot] == 2)] = "red"    #injection drug user
				vertex.size = rep(2,n.plot)
				vertex.size[which(s.plot %in% s$Wave_0)] = 3.5
				G = graph.adjacency(adjacency[c(s.plot),c(s.plot)])
				loc = layout_components(G)
				plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
				vertex.shape=vertex.shape,edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Sample Graph with Seeds Represented as Englarged Nodes",
				xlab="Non-injection drug users as yellow circles, Injection drug users as red circles")
				}
			}
				
		#Calculate the Frank and Snijders estimates
		FS_estimate(s0=s$Wave_0, var_est=var_est)

		#Record the estimates
		N.hat.strata.sim[kk,] = N.hat.strata
		N.hat.strata.var.sim[kk,] = N.hat.strata.var	
		N.hat.strata.1.sim[kk,] = N.hat.strata.1
		N.hat.strata.1.var.sim[kk,] = N.hat.strata.1.var			
		for(k in 1:num_strata)
			{
			S.strata.log[kk,k] = S.strata.1.log[kk,k] = length(which(pop$strata[s$Wave_0] == k))
			C.strata.log[kk,k] = exp(1.96*(log(1+N.hat.strata.var.sim[kk,k]/(N.hat.strata.sim[kk,k]-S.strata.log[kk,k])^2))^(0.5))
			C.strata.1.log[kk,k] = exp(1.96*(log(1+N.hat.strata.1.var.sim[kk,k]/(N.hat.strata.1.sim[kk,k]-S.strata.1.log[kk,k])^2))^(0.5))
			}
			
		#Calculate the new style estimates
		FS_new_estimate(s0=s$Wave_0, var_est_new=var_est_new)

		#Record the estimates
		N.hat.strata.new.sim[kk,] = N.hat.strata.new
		N.hat.strata.new.var.sim[kk,] = N.hat.strata.new.var
		N.hat.strata.1.new.sim[kk,] = N.hat.strata.1.new
		N.hat.strata.1.new.var.sim[kk,] = N.hat.strata.1.new.var
		for(k in 1:num_strata)
			{
			S.strata.new.log[kk,k] = S.strata.1.new.log[kk,k] = length(which(pop$strata[s$Wave_0] == k))
			C.strata.new.log[kk,k] = exp(1.96*(log(1+N.hat.strata.new.var.sim[kk,k]/(N.hat.strata.new.sim[kk,k]-S.strata.new.log[kk,k])^2))^(0.5))
			C.strata.1.new.log[kk,k] = exp(1.96*(log(1+N.hat.strata.1.new.var.sim[kk,k]/(N.hat.strata.1.new.sim[kk,k]-S.strata.1.new.log[kk,k])^2))^(0.5))
			}

		#Calculate the Rao-Blackwellized estimates
		if(num_mcmc > 0)
			{
			MCMC(s=s, num_mcmc, num_chains, search_iterations, use_beta_est, var_est, var_est_new)
			N.hat.RB.strata.sim[kk,] = apply(N.hat.strata.mcmc,2,mean)
			N.hat.RB.strata.var.sim[kk,] = apply(N.hat.strata.var.mcmc,2,mean) - apply(N.hat.strata.mcmc,2,var)
			N.hat.RB.strata.var.sim.negative[kk,(N.hat.RB.strata.var.sim[kk,] < 0)] = 1
			N.hat.RB.strata.var.sim[kk,(N.hat.RB.strata.var.sim[kk,] < 0)] = apply(N.hat.strata.var.mcmc,2,mean)[(N.hat.RB.strata.var.sim[kk,] < 0)]				
			N.hat.RB.strata.1.sim[kk,] = apply(N.hat.strata.1.mcmc,2,mean)
			N.hat.RB.strata.1.var.sim[kk,] = apply(N.hat.strata.1.var.mcmc,2,mean) - apply(N.hat.strata.1.mcmc,2,var)
			N.hat.RB.strata.1.var.sim.negative[kk,(N.hat.RB.strata.1.var.sim[kk,] < 0)] = 1
			N.hat.RB.strata.1.var.sim[kk,(N.hat.RB.strata.1.var.sim[kk,] < 0)] = apply(N.hat.strata.1.var.mcmc,2,mean)[(N.hat.RB.strata.1.var.sim[kk,] < 0)]				
			for(k in 1:num_strata)
				{
				S.RB.strata.log[kk,k] = S.RB.strata.1.log[kk,k] = length(which(pop$strata[s$Wave_0] == k))
				C.RB.strata.log[kk,k] = exp(1.96*(log(1+N.hat.RB.strata.var.sim[kk,k]/(N.hat.RB.strata.sim[kk,k]-S.RB.strata.log[kk,k])^2))^(0.5))
				C.RB.strata.1.log[kk,k] = exp(1.96*(log(1+N.hat.RB.strata.1.var.sim[kk,k]/(N.hat.RB.strata.1.sim[kk,k]-S.RB.strata.1.log[kk,k])^2))^(0.5))
				}
			N.hat.RB.strata.new.sim[kk,] = apply(N.hat.strata.new.mcmc, 2, mean)
			N.hat.RB.strata.new.var.sim[kk,] = apply(N.hat.strata.new.var.mcmc,2,mean) - apply(N.hat.strata.new.mcmc,2,var)		
			N.hat.RB.strata.new.var.sim.negative[kk,(N.hat.RB.strata.new.var.sim[kk,] < 0)] = 1
			N.hat.RB.strata.new.var.sim[kk,(N.hat.RB.strata.new.var.sim[kk,] < 0)] = apply(N.hat.strata.new.var.mcmc,2,mean)[(N.hat.RB.strata.new.var.sim[kk,] < 0)]	
			N.hat.RB.strata.1.new.sim[kk,] = apply(N.hat.strata.1.new.mcmc, 2, mean)
			N.hat.RB.strata.1.new.var.sim[kk,] = apply(N.hat.strata.1.new.var.mcmc,2,mean) - apply(N.hat.strata.1.new.mcmc,2,var)		
			N.hat.RB.strata.1.new.var.sim.negative[kk,(N.hat.RB.strata.1.new.var.sim[kk,] < 0)] = 1
			N.hat.RB.strata.1.new.var.sim[kk,(N.hat.RB.strata.1.new.var.sim[kk,] < 0)] = apply(N.hat.strata.1.new.var.mcmc,2,mean)[(N.hat.RB.strata.1.new.var.sim[kk,] < 0)]	
			for(k in 1:num_strata)
				{
				S.RB.strata.new.log[kk,k] = S.RB.strata.1.new.log[kk,k] = length(which(pop$strata[s$Wave_0] == k))
				C.RB.strata.new.log[kk,k] = exp(1.96*(log(1+N.hat.RB.strata.new.var.sim[kk,k]/(N.hat.RB.strata.new.sim[kk,k]-S.RB.strata.new.log[kk,k])^2))^(0.5))
				C.RB.strata.1.new.log[kk,k] = exp(1.96*(log(1+N.hat.RB.strata.1.new.var.sim[kk,k]/(N.hat.RB.strata.1.new.sim[kk,k]-S.RB.strata.1.new.log[kk,k])^2))^(0.5))
				}
			accept_mcmc.sim[kk] = mean(accept_mcmc)
			}	
		if(write.data==TRUE)
			{		
			xx = c(seed,num_strata,use_beta_est,var_est,var_est_new,alpha,as.numeric(beta_strata),beta_sex_worker,beta_net_size,
			num_waves,n.sim[kk,],n.sex.worker.sim[kk,],
			accept_mcmc.sim[kk],num_mcmc,num_chains,search_iterations,
			N.hat.strata.sim[kk,],N.hat.strata.var.sim[kk,],N.hat.strata.1.sim[kk,],N.hat.strata.1.var.sim[kk,],
			S.strata.log[kk,],C.strata.log[kk,],S.strata.1.log[kk,],C.strata.1.log[kk,],
			N.hat.strata.new.sim[kk,],N.hat.strata.new.var.sim[kk,],N.hat.strata.1.new.sim[kk,],N.hat.strata.1.new.var.sim[kk,],
			S.strata.new.log[kk,],C.strata.new.log[kk,],S.strata.1.new.log[kk,],C.strata.1.new.log[kk,],				
			N.hat.RB.strata.sim[kk,],N.hat.RB.strata.var.sim[kk,],N.hat.RB.strata.var.sim.negative[kk,],
			N.hat.RB.strata.1.sim[kk,],N.hat.RB.strata.1.var.sim[kk,],N.hat.RB.strata.1.var.sim.negative[kk,],
			S.RB.strata.log[kk,],C.RB.strata.log[kk,],S.RB.strata.1.log[kk,],C.RB.strata.1.log[kk,],
			N.hat.RB.strata.new.sim[kk,],N.hat.RB.strata.new.var.sim[kk,],N.hat.RB.strata.new.var.sim.negative[kk,],
			N.hat.RB.strata.1.new.sim[kk,],N.hat.RB.strata.1.new.var.sim[kk,],N.hat.RB.strata.1.new.var.sim.negative[kk,],
			S.RB.strata.new.log[kk,],C.RB.strata.new.log[kk,],S.RB.strata.1.new.log[kk,],C.RB.strata.1.new.log[kk,])
			filename=paste("CSpop",CSpop,"Sim_parameters",num_strata,sep="_")	
			filename=paste(filename,var_est,sep="_")
			filename=paste(filename,var_est_new,sep="_")
			for(k in 1:num_strata)
				filename=paste(filename,alpha[k],sep="_")
			filename=paste(filename, use_beta_est,sep="_")
			for(k in 1:(num_strata*num_strata))
				filename=paste(filename,as.numeric(beta_strata)[k],sep="_")
			filename=paste(filename,use_beta_est,beta_net_size,num_waves,num_chains,num_mcmc,sep="_")
			filename=paste(filename,".csv", sep="")
			write(x = xx, file = filename, ncolumns = length(xx), append=TRUE, sep=",")	
			}	
		#Print the progress of the simulation study
		setTxtProgressBar(progress_sim, kk)
		}
	n.sim <<- n.sim
	n.sex.worker.sim <<- n.sex.worker.sim
	
	N.hat.strata.sim <<- N.hat.strata.sim
	N.hat.strata.var.sim <<- N.hat.strata.var.sim
	
	N.hat.strata.1.sim <<- N.hat.strata.1.sim
	N.hat.strata.1.var.sim <<- N.hat.strata.1.var.sim
	
	S.strata.log <<- S.strata.log
	C.strata.log <<- C.strata.log
	
	S.strata.1.log <<- S.strata.1.log
	C.strata.1.log <<- C.strata.1.log
	
	N.hat.strata.new.sim <<- N.hat.strata.new.sim
	N.hat.strata.new.var.sim <<- N.hat.strata.new.var.sim
	
	N.hat.strata.1.new.sim <<- N.hat.strata.1.new.sim
	N.hat.strata.1.new.var.sim <<- N.hat.strata.1.new.var.sim
	
	S.strata.new.log <<- S.strata.new.log
	C.strata.new.log <<- C.strata.new.log
	
	S.strata.1.new.log <<- S.strata.1.new.log
	C.strata.1.new.log <<- C.strata.1.new.log	
		
	N.hat.RB.strata.sim <<- N.hat.RB.strata.sim
	N.hat.RB.strata.var.sim <<- N.hat.RB.strata.var.sim
	N.hat.RB.strata.var.sim.negative <<- N.hat.RB.strata.var.sim.negative
	
	N.hat.RB.strata.1.sim <<- N.hat.RB.strata.1.sim
	N.hat.RB.strata.1.var.sim <<- N.hat.RB.strata.1.var.sim
	N.hat.RB.strata.1.var.sim.negative <<- N.hat.RB.strata.1.var.sim.negative
	
	S.RB.strata.log <<- S.RB.strata.log
	C.RB.strata.log <<- C.RB.strata.log
	
	S.RB.strata.1.log <<- S.RB.strata.1.log
	C.RB.strata.1.log <<- C.RB.strata.1.log
	
	N.hat.RB.strata.new.sim <<- N.hat.RB.strata.new.sim
	N.hat.RB.strata.new.var.sim <<- N.hat.RB.strata.new.var.sim
	N.hat.RB.strata.new.var.sim.negative <<- N.hat.RB.strata.new.var.sim.negative
	
	N.hat.RB.strata.1.new.sim <<- N.hat.RB.strata.1.new.sim
	N.hat.RB.strata.1.new.var.sim <<- N.hat.RB.strata.1.new.var.sim
	N.hat.RB.strata.1.new.var.sim.negative <<- N.hat.RB.strata.1.new.var.sim.negative
	
	S.RB.strata.new.log <<- S.RB.strata.new.log
	C.RB.strata.new.log <<- C.RB.strata.new.log
	
	S.RB.strata.1.new.log <<- S.RB.strata.1.new.log
	C.RB.strata.1.new.log <<- C.RB.strata.1.new.log
	
	accept_mcmc.sim <<- accept_mcmc.sim	
	}
	