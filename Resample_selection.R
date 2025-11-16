#The "Resample_selection" function selects a sample reordering by first selecting a cluster with probability proportional to size. It then
#selects a random number with probability inversely proportional to size (up to the size of the cluster) and assigns a randomly chosen set of individuals 
#of this size to a wave for membership completely at random while ensuring that the same number show up in the initial sample. 
#The proposal distribution is symmetric between two points and calculations are therefore not required.

Resample_selection <- function(s.resample, beta_strata_calc, beta_sex_worker_calc, beta_net_size_calc)
	{
	#The proportional probability of selecting the sample reordering in the full graph setting will be denoted as Q.prob and will be in the log scale
	Q.prob = log(1)
		
	#Choose a cluster to reorder (selection is done with probability equal to the size of the cluster)
	cluster_select = sample(Graph_s_numclusters,1,replace=FALSE,prob=(Graph_s_clusters_size-1))
	cluster_select_elements = as.numeric(names(which(Graph_s_clusters == cluster_select)))

	cluster_wave_membership = list()
	for(j in 0:num_waves)	
		{	
		cluster_wave_membership[[j+1]] = cluster_select_elements[which(cluster_select_elements %in% s.resample[[j+1]])]
		names(cluster_wave_membership)[j+1] = paste("Wave", j, sep="_")
		}			

	cluster_wave_membership_reordered = cluster_wave_membership
		
	for(k in 1:num_strata)
		{
		cluster_wave_membership_stratum = unlist(cluster_wave_membership)[which(pop$strata[unlist(cluster_wave_membership)]==k)]	
		names(cluster_wave_membership_stratum) = NULL	
		
		size_cluster_stratum = length(cluster_wave_membership_stratum) 
		
		#The approach here is to select an arbitrary number of individuals from the cluster and by stratum to "permute" across waves
		if(size_cluster_stratum == 1 && cluster_wave_membership_stratum %in% s.resample$Wave_0 == FALSE && num_waves > 1)
			{
			elements_permute = cluster_wave_membership_stratum
			elements_permute_wave = sample(2:(num_waves+1),1,replace=FALSE)
			for(kk in 1:(num_waves+1))
				cluster_wave_membership_reordered[[kk]] = setdiff(cluster_wave_membership_reordered[[kk]], elements_permute)
						
			cluster_wave_membership_reordered[[elements_permute_wave]] = c(cluster_wave_membership_reordered[[elements_permute_wave]], elements_permute)					
			}		
		
		#The approach here is to select an arbitrary number of individuals from the cluster and by stratum to "permute" across waves	
		if(size_cluster_stratum > 1)
			{			
			sample_size_cluster_stratum = 
				sample(1:size_cluster_stratum,1,replace=FALSE, prob=dgeom(1:size_cluster_stratum,prob=(num_waves/(num_waves+1))))
			#cat("sample_size_cluster_stratum is", sample_size_cluster_stratum,'\n')
			elements_permute = sample(cluster_wave_membership_stratum, sample_size_cluster_stratum, replace=FALSE)
			num_initial_elements_permute = length(which(elements_permute %in% cluster_wave_membership$Wave_0))
			elements_permute_wave = numeric(sample_size_cluster_stratum)
			if(num_initial_elements_permute > 0)
				{
				if(num_initial_elements_permute==1 && sample_size_cluster_stratum == 1)
					elements_permute_wave = 1
				if(sample_size_cluster_stratum > 1)
					elements_permute_wave[sample(1:sample_size_cluster_stratum,num_initial_elements_permute,replace=FALSE)] = 1
				}
						
			if(num_waves==1)
				elements_permute_wave[which(elements_permute_wave==0)] = 2	
			if(num_waves > 1)
				elements_permute_wave[which(elements_permute_wave==0)] = sample(2:(num_waves+1), sample_size_cluster_stratum-num_initial_elements_permute, replace=TRUE)
						
			for(kk in 1:(num_waves+1))
				cluster_wave_membership_reordered[[kk]] = setdiff(cluster_wave_membership_reordered[[kk]], elements_permute)
						
			for(kk in 1:sample_size_cluster_stratum)
				{
				cluster_wave_membership_reordered[[elements_permute_wave[kk]]] = 
					c(cluster_wave_membership_reordered[[elements_permute_wave[kk]]], elements_permute[kk])
				}
			}
		}

	#Determine the empirical probability of obtaining this reordering
	#First check if stratum sample sizes for the initial wave are the same as the original
	for(k in 1:num_strata)
			{
			strata_wave_count_original = length(which(pop$strata[cluster_wave_membership[[1]]] == k))
			strata_wave_count_permuted = length(which(pop$strata[cluster_wave_membership_reordered[[1]]] == k))
			if(strata_wave_count_original != strata_wave_count_permuted)
				Q.prob = log(0)
			}	

	#If the proposed sample reordering is identical to the most recently accepted reordering then skip ahead to another proposal	
	skip.propose = ifelse(prod(compare.list(cluster_wave_membership, cluster_wave_membership_reordered)==1),1,0)			

	if(Q.prob > log(0) && skip.propose == 0)
		{	
		#Configure the proposed sample reordering
		s.resample.proposed = s.resample
		for(k in 1:(num_waves+1))
			{
			s.resample.proposed[[k]] = unique(c(setdiff(s.resample.proposed[[k]], cluster_wave_membership[[k]]), cluster_wave_membership_reordered[[k]]))
			names(s.resample.proposed[[k]]) = NULL
			}

		#Determine the proportional empirical log probability of obtaining this reordering
		s.cluster = as.numeric(names(which(Graph_s_clusters == cluster_select)))
		s.cluster.list = list()
				
		for(k in 1:(num_waves+1))
			s.cluster.list[[k]] = s.cluster[which(s.cluster %in% s.resample.proposed[[k]])]
				
		for(k in 2:(num_waves+1))
			{
			candidates.cluster = setdiff(which(matrix(pop$x[s.cluster.list[[k-1]],],nrow=length(s.cluster.list[[k-1]]))==1,arr.ind=TRUE)[,2], 
									unlist(s.cluster.list[1:(k-1)]))
			candidates.cluster_selected = intersect(candidates.cluster,s.cluster.list[[k]])
			if(prod(s.cluster.list[[k]] %in% candidates.cluster_selected) == 0)
				Q.prob = log(0)
			candidates.cluster_not_selected = setdiff(candidates.cluster,s.cluster.list[[k]])		
			for(i in candidates.cluster_selected)
				{
				strata_i = pop$strata[[i]]		
				ngbrs_i = s.cluster.list[[k-1]][which(matrix(pop$x[s.cluster.list[[k-1]],i],nrow=length(s.cluster.list[[k-1]]))==1)]
				strata_ngbrs_i = pop$strata[ngbrs_i]
				sex_worker_ngrs_i = pop$sex.worker[ngbrs_i]
				if(length(ngbrs_i) == 0)
					Q.prob = log(0) else
					net_size_ngbrs_i = pop$net_size[ngbrs_i] 
				Q.prob = Q.prob + log(1 - prod(1-inv.logit(beta_strata_calc[strata_ngbrs_i,strata_i] + 
												ifelse(is.na(beta_sex_worker),0,beta_sex_worker_calc*sex_worker_ngrs_i) +
												beta_net_size_calc*log(net_size_ngbrs_i))))
				}
				
			for(i in candidates.cluster_not_selected)
				{
				strata_i = pop$strata[[i]]		
				ngbrs_i = s.cluster.list[[k-1]][which(matrix(pop$x[s.cluster.list[[k-1]],i],nrow=length(s.cluster.list[[k-1]]))==1)]
				strata_ngbrs_i = pop$strata[ngbrs_i]
				sex_worker_ngrs_i = pop$sex.worker[ngbrs_i]
				if(length(ngbrs_i) == 0)
					Q.prob = log(0) else
					net_size_ngbrs_i = pop$net_size[ngbrs_i] 	
				Q.prob = Q.prob + log(prod(1-inv.logit(beta_strata_calc[strata_ngbrs_i,strata_i] + 
										ifelse(is.na(beta_sex_worker),0,beta_sex_worker_calc*sex_worker_ngrs_i) +
										beta_net_size_calc*log(net_size_ngbrs_i))))
				}	
			}
		}

	Q.prob <<- Q.prob
	skip.propose <<- skip.propose
	if(Q.prob > log(0)  && skip.propose == 0)
		{
		cluster_select <<- cluster_select
		s.resample.proposed <<- s.resample.proposed
		}
	}	
