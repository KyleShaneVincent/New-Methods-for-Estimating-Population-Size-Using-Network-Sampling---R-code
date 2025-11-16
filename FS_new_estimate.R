#The "FS_new_estimate" function will calculate the new Frank and Snijders (1994) style estimates for each stratum.

FS_new_estimate <- function(s0, var_est_new)
	{
	#Determine the size of the initial sample
	n0 = length(s0)
		
	#Determine the initial sample size for each stratum 
	s0_strata = list()
	initial.size.strata = numeric()
	for(k in 1:num_strata)
		{
		s0_strata[[k]] = intersect(s0,pop.strata[[k]])
		initial.size.strata[k] = length(s0_strata[[k]])
		}
		
	#The set difference of stratum k with initial sample of stratum k
	strata_setdiff_s0 = list()
	for(k in 1:num_strata)
		strata_setdiff_s0[[k]] = setdiff(pop.strata[[k]],s0_strata[[k]])
			
	#For the stratum size estimator
	r = matrix(0,num_strata,num_strata) #The number of links from stratum k initial sample to stratum j initial sample
	r.1 = list() #The number of individuals within stratum k initial sample linked to at least one other individual within initial sample
	r.1.records = list()
	for(k in 1:num_strata)
		{
		r.1.records = append(r.1.records, list(list()))	
		r.1 = append(r.1, list(list()))
		r.1[[k]] = numeric()
		}
	for(k in 1:num_strata)
		{
		for(j in 1:num_strata)
			{
			r.1.records[[k]][[j]] = apply(pop$x[s0_strata[[k]],s0_strata[[j]]],2,sum)
			r[k,j] = sum(r.1.records[[k]][[j]])
			r.1[[j]] = unlist(c(r.1[[j]], names(which(r.1.records[[k]][[j]] > 0))))
			}
		}
	for(k in 1:num_strata)
		r.1[[k]] = length(unique(r.1[[k]]))	

	l = matrix(0,num_strata,num_strata) #The number of links from stratum k initial sample to outside the initial sample and stratum j
	l.1 = list() #The number of individuals in stratum k and outside initial sample that are linked to at least one other individual within initial sample
	l.1.records = list()
	for(k in 1:num_strata)
		{
		l.1.records = append(l.1.records, list(list()))	
		l.1 = append(l.1, list(list()))		
		l.1[[k]] = numeric()
		}
	for(k in 1:num_strata)
		for(j in 1:num_strata)
			{
			l.1.records[[k]][[j]] = apply(pop$x[s0_strata[[k]],strata_setdiff_s0[[j]]],2,sum)
			l[k,j] = sum(l.1.records[[k]][[j]])
			l.1[[j]] = unlist(c(l.1[[j]], names(which(l.1.records[[k]][[j]] > 0))))
			}
	for(k in 1:num_strata)
		l.1[[k]] = length(unique(l.1[[k]]))		
	
	#Stratum population size estimator
	N.hat.strata.new = numeric()
	N.hat.strata.1.new = numeric()
	for(k in 1:num_strata)
		{
		r.k = sum(r[,k])
		l.k = sum(l[,k])
		N.hat.strata.new[k] = (initial.size.strata[k] + 1) * (r.k + l.k + 1) / (r.k + 1) - 1
		
		r.1.k = r.1[[k]]
		l.1.k = l.1[[k]]
		N.hat.strata.1.new[k] = (initial.size.strata[k] + 1) * (r.1.k + l.1.k + 1) / (r.1.k + 1) - 1
		}	
	
	#Seber-type method for variance estimation	
	if(var_est_new=="Seber")
		{
		#Stratum population size variance estimator
		N.hat.strata.new.var = numeric()
		N.hat.strata.1.new.var = numeric()
		for(k in 1:num_strata)
			{
			#Adapted from Seber for link-count estimator
			r.k = r.1[[k]]
			l.k = sum(l[,k]) * ifelse(r.1[[k]]==0,1,r.1[[k]])/ifelse(sum(r[,k])==0,1,sum(r[,k]))
			N.hat.strata.new.var[k] = (initial.size.strata[k] + 1) * (r.k + l.k + 1) * (initial.size.strata[k] - r.k) *  l.k / ((r.k + 1)^2 * (r.k + 2))	
			
			#An approximation method for node-count estimator, resembles Seber estimator
			r.1.k = r.1[[k]]
			l.1.k = l.1[[k]]			
			N.hat.strata.1.new.var[k] = (initial.size.strata[k] + 1) * (r.1.k + l.1.k + 1) * (initial.size.strata[k] - r.1.k) *  l.1.k / ((r.1.k + 1)^2 * (r.1.k + 2))
			}
		}
	
	#Jackknife method for variance estimation		
	if(var_est_new=="jack")	
		{
		N.hat.strata.new.jack = matrix(NA, num_strata, n0)		
		N.hat.strata.new.var = numeric()	
		N.hat.strata.1.new.jack = matrix(NA, num_strata, n0)		
		N.hat.strata.1.new.var = numeric()	
		jj = 0
		for(i in s0)
			{	
			jj = jj+1
		
			#The stratum membership and index location of this element
			stratum_i = pop$strata[i]	
			i_index = which(names(r.1.records[[stratum_i]][[stratum_i]]) == i)		

			#N.hat.strata.new estimators with this element removed from the initial sample
			r.jack = r 	
			r.jack[stratum_i,stratum_i] = r[stratum_i,stratum_i] - r.1.records[[stratum_i]][[stratum_i]][i_index] - sum(pop$x[i,s0_strata[[stratum_i]]])
			for(k in setdiff(1:num_strata, stratum_i))
					{
					r.jack[k,stratum_i] = r[k,stratum_i] - sum(pop$x[s0_strata[[k]],i])
					r.jack[stratum_i,k] = r[stratum_i,k] - sum(pop$x[i,s0_strata[[k]]])
					}

			l.jack = l
			l.jack[stratum_i,stratum_i] = l[stratum_i,stratum_i] - sum(pop$x[i,strata_setdiff_s0[[stratum_i]]]) + sum(pop$x[setdiff(s0_strata[[stratum_i]],i),i])  
			for(k in setdiff(1:num_strata,stratum_i))
				{
				l.jack[k,stratum_i] = l[k,stratum_i] + sum(pop$x[s0_strata[[k]],i])
				l.jack[stratum_i,k] = l[stratum_i,k] - sum(pop$x[i,strata_setdiff_s0[[k]]])
				}
				
			#For the stratum size estimator
			r.1.records.jack = r.1.records
			r.1.jack = list()
			r.1.records.jack[[stratum_i]][[stratum_i]] = r.1.records.jack[[stratum_i]][[stratum_i]][-i_index]
			r.1.records.jack[[stratum_i]][[stratum_i]] = r.1.records.jack[[stratum_i]][[stratum_i]] - pop$x[i,setdiff(s0_strata[[stratum_i]],i)]	
			names(r.1.records.jack[[stratum_i]][[stratum_i]]) = setdiff(s0_strata[[stratum_i]],i)
				
			r.1.jack[[stratum_i]] = numeric()
			for(k in setdiff(1:num_strata, stratum_i))
				{
				r.1.records.jack[[k]][[stratum_i]] = r.1.records.jack[[k]][[stratum_i]][-i_index]
				names(r.1.records.jack[[k]][[stratum_i]]) = setdiff(s0_strata[[stratum_i]],i)
				r.1.records.jack[[stratum_i]][[k]] = r.1.records.jack[[stratum_i]][[k]] - pop$x[i,s0_strata[[k]]]
				names(r.1.records.jack[[stratum_i]][[k]]) = s0_strata[[k]]
				r.1.jack[[k]] = numeric()
				}			
			for(k in 1:num_strata)
				{
				for(j in 1:num_strata)
					r.1.jack[[k]] = unlist(c(r.1.jack[[k]], names(which(r.1.records.jack[[j]][[k]] > 0))))
				r.1.jack[[k]] = length(unique(r.1.jack[[k]]))							
				}

			#For the stratum size estimator
			l.1.records.jack = l.1.records
			l.1.jack = list()			
			l.1.records.jack[[stratum_i]][[stratum_i]] = l.1.records.jack[[stratum_i]][[stratum_i]] - pop$x[i,strata_setdiff_s0[[stratum_i]]]
			l.1.records.jack[[stratum_i]][[stratum_i]] = c(l.1.records.jack[[stratum_i]][[stratum_i]], sum(pop$x[setdiff(s0_strata[[stratum_i]],i),i]))
			names(l.1.records.jack[[stratum_i]][[stratum_i]]) = c(strata_setdiff_s0[[stratum_i]],i)
			l.1.jack[[stratum_i]] = numeric()
			for(k in setdiff(1:num_strata, stratum_i))
				{
				l.1.records.jack[[k]][[stratum_i]] = c(l.1.records.jack[[k]][[stratum_i]], sum(pop$x[s0_strata[[k]],i]))
				names(l.1.records.jack[[k]][[stratum_i]]) = c(strata_setdiff_s0[[stratum_i]],i)				
				l.1.records.jack[[stratum_i]][[k]] = l.1.records.jack[[stratum_i]][[k]] - pop$x[i,strata_setdiff_s0[[k]]]
				names(l.1.records.jack[[stratum_i]][[k]]) = strata_setdiff_s0[[k]]		
				l.1.jack[[k]] = numeric()
				}				
			for(k in 1:num_strata)
				{
				for(j in 1:num_strata)
					l.1.jack[[k]] = unlist(c(l.1.jack[[k]], names(which(l.1.records.jack[[j]][[k]] > 0))))
				l.1.jack[[k]] = length(unique(l.1.jack[[k]]))	
				}					

			#Stratum population size estimator
			for(k in 1:num_strata)
				{
				r.jack.k = sum(r.jack[,k])
				l.jack.k = sum(l.jack[,k])
				N.hat.strata.new.jack[k,jj] = initial.size.strata[k] * (r.jack.k + l.jack.k + 1) / (r.jack.k + 1) - 1
				
				r.1.jack.k = r.1.jack[[k]]
				l.1.jack.k = l.1.jack[[k]]
				N.hat.strata.1.new.jack[k,jj] = initial.size.strata[k] * (r.1.jack.k + l.1.jack.k + 1) / (r.1.jack.k + 1) - 1
				}			
			}
					
		for(k in 1:num_strata)
			{
			N.hat.strata.new.jack.k = na.omit(N.hat.strata.new.jack[k,])
			N.hat.strata.new.jack.mean.k = mean(N.hat.strata.new.jack.k)
			N.hat.strata.new.var[k] = (n0-2)/(2*n0) * sum((N.hat.strata.new.jack.k - N.hat.strata.new.jack.mean.k)^2)
		
			N.hat.strata.1.new.jack.k = na.omit(N.hat.strata.1.new.jack[k,])
			N.hat.strata.1.new.jack.mean.k = mean(N.hat.strata.1.new.jack.k)
			N.hat.strata.1.new.var[k] = (n0-2)/(2*n0) * sum((N.hat.strata.1.new.jack.k - N.hat.strata.1.new.jack.mean.k)^2)			
			}	
		}
		
	#Make the estimates available outside the function
	N.hat.strata.new <<- N.hat.strata.new
	N.hat.strata.new.var <<- N.hat.strata.new.var	
	N.hat.strata.1.new <<- N.hat.strata.1.new
	N.hat.strata.1.new.var <<- N.hat.strata.1.new.var
	}
