#The "MCMC" function will run the MCMC chain to calculate an approximation for the Rao-Blackwellized estimates.

MCMC <- function(s, num_mcmc, num_chains, search_iterations, use_beta_est, var_est, var_est_new)
	{
	#Declare the chain output matrices	
	accept_mcmc = numeric(num_mcmc*num_chains)
	u = runif(num_mcmc*num_chains)
	
	N.hat.strata.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)
	N.hat.strata.var.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)
	N.hat.strata.new.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)
	N.hat.strata.new.var.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)	
	
	N.hat.strata.1.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)
	N.hat.strata.1.var.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)
	N.hat.strata.1.new.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)
	N.hat.strata.1.new.var.mcmc = matrix(NA,num_mcmc*num_chains,num_strata)	
	
	#Extract the disjoint clusters from the networked sample
	s_elements = unlist(s)
	names(s_elements) = NULL
	s_elements = sort(s_elements)
	Graph_s = graph.adjacency(pop$x[s_elements,s_elements])

	components_Graph_s = components(Graph_s)

	Graph_s_clusters <<- components_Graph_s$membership
	Graph_s_clusters_size <<- components_Graph_s$csize
	Graph_s_numclusters <<- components_Graph_s$no
	
	#Use the appropriate values for the beta parameters, i.e., if they are the true or estimated values	
	if(use_beta_est == TRUE)
		{
		beta_strata_calc = beta_estimates_list$Strata_Parameter_Estimates
		if(is.na(beta_sex_worker))
			beta_sex_worker_calc = NA else
			beta_sex_worker_calc = beta_estimates_list$Beta_Sex_Worker_Parameter_Estimate		
		beta_net_size_calc = beta_estimates_list$Log_Net_Size_Parameter_Estimate
		} else
			{
			beta_strata_calc = beta_strata
			if(is.na(beta_sex_worker))
				beta_sex_worker_calc = NA else
				beta_sex_worker_calc = beta_estimates_list$Beta_Sex_Worker_Parameter_Estimate
			beta_net_size_calc = beta_net_size
			}

	#Determine the proportional empirical log probability of obtaining the original sample reordering
	#We have to account for probability of tracing and not tracing links
	#Tracing links to next wave
	P.prob.list = list()
	for(cluster in 1:Graph_s_numclusters)
		{
		P.prob.list[[cluster]] = log(1)
		s.cluster = as.numeric(names(which(Graph_s_clusters == cluster)))
		s.cluster.list = list()
				
		for(k in 1:(num_waves+1))
			s.cluster.list[[k]] = s.cluster[which(s.cluster %in% s[[k]])]			
		
		for(k in 2:(num_waves+1))
			{
			candidates.cluster = setdiff(which(matrix(pop$x[s.cluster.list[[k-1]],],nrow=length(s.cluster.list[[k-1]]))==1,arr.ind=TRUE)[,2], 
									unlist(s.cluster.list[1:(k-1)]))
			candidates.cluster_selected = intersect(candidates.cluster,s.cluster.list[[k]])
			candidates.cluster_not_selected = setdiff(candidates.cluster,s.cluster.list[[k]])		
			for(i in candidates.cluster_selected)
				{
				strata_i = pop$strata[[i]]		
				ngbrs_i = s.cluster.list[[k-1]][which(matrix(pop$x[s.cluster.list[[k-1]],i],nrow=length(s.cluster.list[[k-1]]))==1)]
				strata_ngbrs_i = pop$strata[ngbrs_i]
				sex_worker_ngrs_i = pop$sex.worker[ngbrs_i]
				if(length(ngbrs_i) == 1)
					net_size_ngbrs_i = pop$net_size[ngbrs_i] else
						net_size_ngbrs_i = pop$net_size[ngbrs_i]
				P.prob.list[[cluster]] = P.prob.list[[cluster]] + log(1 - prod(1-inv.logit(beta_strata_calc[strata_ngbrs_i,strata_i] + 
												ifelse(is.na(beta_sex_worker),0,beta_sex_worker_calc*sex_worker_ngrs_i) +
												beta_net_size_calc*log(net_size_ngbrs_i))))
				}
				
			for(i in candidates.cluster_not_selected)
				{
				strata_i = pop$strata[[i]]		
				ngbrs_i = s.cluster.list[[k-1]][which(matrix(pop$x[s.cluster.list[[k-1]],i],nrow=length(s.cluster.list[[k-1]]))==1)]
				strata_ngbrs_i = pop$strata[ngbrs_i]
				sex_worker_ngrs_i = pop$sex.worker[ngbrs_i]
				if(length(ngbrs_i) == 1)
					net_size_ngbrs_i = pop$net_size[ngbrs_i] else
						net_size_ngbrs_i = pop$net_size[ngbrs_i]
				P.prob.list[[cluster]] = P.prob.list[[cluster]]  + log(prod(1-inv.logit(beta_strata_calc[strata_ngbrs_i,strata_i] + 
										ifelse(is.na(beta_sex_worker),0,beta_sex_worker_calc*sex_worker_ngrs_i) +
										beta_net_size_calc*log(net_size_ngbrs_i))))
				}	
			}
		}

	P.prob.list.original = P.prob.list

	#Seeding the chains
	for(m in 1:num_chains)
		{
		if(m == 1)
			{	
			FS_estimate(s0=s$Wave_0, var_est=var_est)
			N.hat.strata.mcmc[1,] = N.hat.strata
			N.hat.strata.var.mcmc[1,] = N.hat.strata.var
			N.hat.strata.1.mcmc[1,] = N.hat.strata.1
			N.hat.strata.1.var.mcmc[1,] = N.hat.strata.1.var
			FS_new_estimate(s0=s$Wave_0, var_est_new=var_est_new)	
			N.hat.strata.new.mcmc[1,] = N.hat.strata.new
			N.hat.strata.new.var.mcmc[1,] = N.hat.strata.new.var	
			N.hat.strata.1.new.mcmc[1,] = N.hat.strata.1.new
			N.hat.strata.1.new.var.mcmc[1,] = N.hat.strata.1.new.var			
			accept_mcmc[1] = 1		
			s.resample = s	
			} else
		if(m > 1)
			{
			P.prob.list = P.prob.list.original
			s.resample = s	
			for(k in 1:search_iterations)
				{
				Resample_selection(s.resample=s.resample, beta_strata_calc=beta_strata_calc, beta_sex_worker_calc=beta_sex_worker_calc, 
									beta_net_size_calc=beta_net_size_calc)
				if(Q.prob > log(0))								
					if(Q.prob < P.prob.list[[cluster_select]])
						{
						P.prob.list[[cluster_select]] = Q.prob
						s.resample = s.resample.proposed
						}
				}				
			FS_estimate(s0=s.resample$Wave_0, var_est=var_est)	
			N.hat.strata.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata
			N.hat.strata.var.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata.var
			N.hat.strata.1.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata.1
			N.hat.strata.1.var.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata.1.var
			FS_new_estimate(s0=s.resample$Wave_0, var_est_new=var_est_new)	
			N.hat.strata.new.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata.new
			N.hat.strata.new.var.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata.new.var	
			N.hat.strata.1.new.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata.1.new
			N.hat.strata.1.new.var.mcmc[(m-1)*num_mcmc+1,] = N.hat.strata.1.new.var				
			accept_mcmc[(m-1)*num_mcmc+1] = 1	
			} 
	
		for(t in ((m-1)*num_mcmc+2):(m*num_mcmc))
			{
			#cat("MCMC iterataion", t, '\n')
			Resample_selection(s.resample=s.resample, beta_strata_calc=beta_strata_calc, beta_sex_worker_calc=beta_sex_worker_calc, 
								beta_net_size_calc=beta_net_size_calc)
			#cat("Q.prob is", Q.prob, '\n')
			#cat("P.prob.list[[cluster_select]] is", P.prob.list[[cluster_select]], '\n')
			#if(Q.prob > log(0))
			#	{
			#	cat("s.resample.proposed is", '\n')
			#	print(s.resample.proposed)
			#	}
			
			if(Q.prob == log(0) || skip.propose == 1)
				{
				N.hat.strata.mcmc[t,] = N.hat.strata.mcmc[t-1,]
				N.hat.strata.var.mcmc[t,] = N.hat.strata.var.mcmc[t-1,]
				N.hat.strata.new.mcmc[t,] = N.hat.strata.new.mcmc[t-1,]
				N.hat.strata.new.var.mcmc[t,] = N.hat.strata.new.var.mcmc[t-1,]
				N.hat.strata.1.mcmc[t,] = N.hat.strata.1.mcmc[t-1,]
				N.hat.strata.1.var.mcmc[t,] = N.hat.strata.1.var.mcmc[t-1,]
				N.hat.strata.1.new.mcmc[t,] = N.hat.strata.1.new.mcmc[t-1,]
				N.hat.strata.1.new.var.mcmc[t,] = N.hat.strata.1.new.var.mcmc[t-1,]				
				if(skip.propose==1)
					accept_mcmc[t] = 1
				} else
					{
					#cat('At the accept/reject step, Q.prob is', Q.prob, 'P.prob is', P.prob, '\n')
					if(log(u[t]) < (Q.prob - P.prob.list[[cluster_select]]))
						{
						#cat("Accepted", '\n')
						s.resample = s.resample.proposed
						#cat("s.resample is ", '\n')
						#print(s.resample)
						P.prob.list[[cluster_select]] = Q.prob
						#cat("Updated P.prob is ",P.prob, '\n')
						FS_estimate(s0=s.resample$Wave_0, var_est=var_est)	
						N.hat.strata.mcmc[t,] = N.hat.strata
						N.hat.strata.var.mcmc[t,] = N.hat.strata.var
						N.hat.strata.1.mcmc[t,] = N.hat.strata.1
						N.hat.strata.1.var.mcmc[t,] = N.hat.strata.1.var
						FS_new_estimate(s0=s.resample$Wave_0, var_est_new=var_est_new)	
						N.hat.strata.new.mcmc[t,] = N.hat.strata.new
						N.hat.strata.new.var.mcmc[t,] = N.hat.strata.new.var
						N.hat.strata.1.new.mcmc[t,] = N.hat.strata.1.new
						N.hat.strata.1.new.var.mcmc[t,] = N.hat.strata.1.new.var						
						accept_mcmc[t] = 1
						} else
						{
						#cat("Rejected", '\n')
						N.hat.strata.mcmc[t,] = N.hat.strata.mcmc[t-1,]
						N.hat.strata.var.mcmc[t,] = N.hat.strata.var.mcmc[t-1,]
						N.hat.strata.new.mcmc[t,] = N.hat.strata.new.mcmc[t-1,]
						N.hat.strata.new.var.mcmc[t,] = N.hat.strata.new.var.mcmc[t-1,]
						N.hat.strata.1.mcmc[t,] = N.hat.strata.1.mcmc[t-1,]
						N.hat.strata.1.var.mcmc[t,] = N.hat.strata.1.var.mcmc[t-1,]
						N.hat.strata.1.new.mcmc[t,] = N.hat.strata.1.new.mcmc[t-1,]
						N.hat.strata.1.new.var.mcmc[t,] = N.hat.strata.1.new.var.mcmc[t-1,]
						}
					}			
			}
		}
	N.hat.strata.mcmc <<- N.hat.strata.mcmc 
	N.hat.strata.var.mcmc <<- N.hat.strata.var.mcmc
	N.hat.strata.new.mcmc <<- N.hat.strata.new.mcmc
	N.hat.strata.new.var.mcmc <<- N.hat.strata.new.var.mcmc
	N.hat.strata.1.mcmc <<- N.hat.strata.1.mcmc 
	N.hat.strata.1.var.mcmc <<- N.hat.strata.1.var.mcmc
	N.hat.strata.1.new.mcmc <<- N.hat.strata.1.new.mcmc
	N.hat.strata.1.new.var.mcmc <<- N.hat.strata.1.new.var.mcmc	
	accept_mcmc <<- accept_mcmc
	}
	