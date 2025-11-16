#The "Sample_selection" function selects an initial sample according to a Bernoulli sampling design applied independently within each stratum and then 
#links are traced independently from each wave with a logit-based probability.

Sample_selection <- function(alpha, beta_strata, beta_sex_worker, beta_net_size, use_beta_est, return_beta_est)
	{
	#The sampled units by wave will be recorded as a list
	s = list()

	#Select the initial sample
	s[[1]] = numeric()
	for(k in 1:num_strata)
		{
		s0_candidate = which(pop$strata == k)
		u = runif(length(s0_candidate))
		s0_k = s0_candidate[which(u<alpha[k])]
		if(length(s0_k) < 2)	#Always make sure that there are at least two individuals per stratum in the initial sample
			s0_k = sample(s0_candidate,2,replace=FALSE)		
		s[[1]] = c(s[[1]], s0_k)
		}
	names(s)[1] = "Wave_0"	

	#adjacency.plot will be the adjacency matrix required for plotting the traced links, it starts with all zeroes 
	R = nrow(adjacency)
	C = ncol(adjacency)
	adjacency.plot = matrix(0,R,C)

	#Set up the logistic regression matrix
	log_reg_matrix = list()
	j = 0
	
	#Select wave k
	for(k in 1:num_waves)
		{	
		s[[k+1]] = numeric()
		if(length(s[[k]]) > 0)
			{
			#The elements selected at the next wave 
			#Note that coupons can still be passed to those already observed up to and including the current wave, such observations are not used
			#in the estimation/MCMC procedure but can be used to improve GLM fit
			selected_units = list()
			for(i in s[[k]])
				{
				j = j+1
				i_candidates = which(pop$x[i,] == 1)
				selected_units[[j]] = numeric()		
				num_candidates = length(i_candidates)
				log_reg_matrix[[j]] = matrix(0,num_candidates,5) 
				if(num_candidates > 0)
					{
					i_net_size = pop$net_size[i]
					u.i_candidates = runif(num_candidates)
					strata_i = pop$strata[i]
					strata_i_candidates = pop$strata[i_candidates]
					sex_worker_i_candidates = pop$sex.worker[i_candidates]
					candidates.select = which(u.i_candidates < inv.logit(beta_strata[strata_i,strata_i_candidates] + 
												ifelse(is.na(beta_sex_worker),0,beta_sex_worker*sex_worker_i_candidates) +
												beta_net_size*log(i_net_size)))
					selected_units[[j]] = i_candidates[candidates.select]	
					adjacency.plot[i,selected_units[[j]]] = 1			
					log_reg_matrix[[j]][candidates.select,1] = 1
					log_reg_matrix[[j]][,2] = paste(strata_i, strata_i_candidates, sep="_")
					log_reg_matrix[[j]][,3] = sex_worker_i_candidates
					log_reg_matrix[[j]][,4] = log(i_net_size)
					log_reg_matrix[[j]][,5] = paste(i,i_candidates,sep="_") 
					}
				}			
			selected_units = unlist(selected_units)			
			s_k = unlist(s[1:k])
			s[[k+1]] = setdiff(unique(selected_units), s_k)
			names(s)[k+1] = paste("Wave", k, sep="_")			
			}
		names(s)[k+1] = paste("Wave", k, sep="_")	
		}

	if(use_beta_est==TRUE)
		{
		log_reg_matrix_setup = numeric()
		for(k in 1:length(log_reg_matrix))
			log_reg_matrix_setup = rbind(log_reg_matrix_setup, log_reg_matrix[[k]])

		colnames(log_reg_matrix_setup) = c("Selection", "Strata_combination", "Receiver_sex_worker", "Sender_log_net_size", "Sender_Receiver")
		log_reg_matrix_setup = data.frame(log_reg_matrix_setup)
		
		#In case there is an absence of observations between two strata, add an untraced observation
		for(k in 1:num_strata)
			for(j in 1:num_strata)
				if(length(which(log_reg_matrix_setup$Strata_combination == paste(k,j,sep="_")))==0)
					log_reg_matrix_setup = rbind(log_reg_matrix_setup,c(0,paste(k,j,sep="_"),0,0,"0_0"))		
		log_reg_matrix_setup[,"Selection"] = as.factor(log_reg_matrix_setup[,"Selection"])
		log_reg_matrix_setup[,"Strata_combination"] = as.factor(log_reg_matrix_setup[,"Strata_combination"])
		log_reg_matrix_setup[,"Receiver_sex_worker"] = as.factor(log_reg_matrix_setup[,"Receiver_sex_worker"])		
		log_reg_matrix_setup[,"Sender_log_net_size"] = as.numeric(log_reg_matrix_setup[,"Sender_log_net_size"])				
				
		glm_formula = ifelse(is.na(beta_sex_worker),"Selection ~ -1 + Strata_combination + Sender_log_net_size",
								"Selection ~ -1 + Strata_combination + Receiver_sex_worker + Sender_log_net_size")
			
		glm_model = glm(glm_formula, family=binomial, data=log_reg_matrix_setup)
		glm_output = summary(glm_model)
		x = which(glm_output$coefficients[,"Pr(>|z|)"] >= 0.05)
		beta_estimates = glm_output$coefficients[,1]
		beta_estimates[x] = 0
		if(return_beta_est==TRUE)
			print(glm_output)
		}

	s <<- s
	adjacency.plot <<- adjacency.plot
	if(use_beta_est==TRUE)
		{
		k = 1
		beta_estimates_list = list()
		beta_estimates_list[[k]] = matrix(beta_estimates[1:(num_strata*num_strata)],nrow=num_strata,ncol=num_strata,byrow=TRUE)
		names(beta_estimates_list)[k] = "Strata_Parameter_Estimates"
		if(!is.na(beta_sex_worker))
			{
			k = k+1
			beta_estimates_list[[k]] = beta_estimates[(num_strata*num_strata) + 1]
			names(beta_estimates_list)[k] = "Beta_Sex_Worker_Parameter_Estimate"	
			}
		k = k+1
		beta_estimates_list[[k]] = beta_estimates[(num_strata*num_strata) + 1 + !is.na(beta_sex_worker)]
		names(beta_estimates_list)[k] = "Log_Net_Size_Parameter_Estimate"
		beta_estimates_list <<- beta_estimates_list
		log_reg_matrix_setup <<- log_reg_matrix_setup
		glm_model <<- glm_model
		glm_output <<- glm_output
		}	
	}