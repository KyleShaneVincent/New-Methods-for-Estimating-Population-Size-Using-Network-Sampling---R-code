#The "FS_estimate" function will calculate the Frank and Snijders (1994) estimates for each stratum.

FS_estimate <- function(s0, var_est)
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
	r = numeric() #The number of links within stratum k initial sample
	r.1 = numeric()	#The number of individuals within stratum k initial sample linked to at least one other individual within stratum k initial sample
	r.1.records = list()
	for(k in 1:num_strata)
		{
		r.1.records[[k]] = apply(pop$x[s0_strata[[k]],s0_strata[[k]]],2,sum)
		r[k] = sum(r.1.records[[k]])
		r.1[k] = length(which(r.1.records[[k]] > 0)) 
		}

	l = numeric() #The number of links from stratum k initial sample to outside the initial sample and stratum k
	l.1 = numeric() #The number of individuals in stratum k and outside initial sample that are linked to at least one other individual within stratum k initial sample
	l.1.records = list()
	for(k in 1:num_strata)
		{
		l.1.records[[k]] = apply(pop$x[s0_strata[[k]],strata_setdiff_s0[[k]]],2,sum)
		l[k] = sum(l.1.records[[k]])
		l.1[k] = length(which(l.1.records[[k]] > 0))
		}
		
	#Stratum population size estimator
	N.hat.strata = numeric()
	N.hat.strata.1 = numeric()
	for(k in 1:num_strata)
		{
		initial.size.strata.k=initial.size.strata[k]; r.k=r[k]; l.k=l[k]; r.1.k=r.1[k]; l.1.k=l.1[k]
		N.hat.strata[k] = (initial.size.strata.k + 1) * (r.k + l.k + 1) / (r.k + 1) - 1
		N.hat.strata.1[k] = (initial.size.strata.k + 1) * (r.1.k + l.1.k + 1) / (r.1.k + 1) - 1
		}		
		
	#Seber-type method for variance estimation	
	if(var_est=="Seber")
		{
		#Stratum population size variance estimator
		N.hat.strata.var = numeric()
		N.hat.strata.1.var = numeric()
		for(k in 1:num_strata)
			{
			#Adapted from Seber for link-count estimator
			r.k = r.1[k]
			l.k = l[k] * ifelse(r.1[k]==0,1,r.1[k])/ifelse(r[k]==0,1,r[k])
			N.hat.strata.var[k] = (initial.size.strata[k] + 1) * (r.k + l.k + 1) * (initial.size.strata[k] - r.k) *  l.k / ((r.k + 1)^2 * (r.k + 2))		
			
			#An approximation method for node-count estimator, resembles Seber estimator
			r.1.k = r.1[k]
			l.1.k = l.1[k]
			N.hat.strata.1.var[k] = (initial.size.strata[k] + 1) * (r.1.k + l.1.k + 1) * (initial.size.strata[k] - r.1.k) *  l.1.k / ((r.1.k + 1)^2 * (r.1.k + 2))
			}		
		}
			
	#Jackknife method for variance estimation		
	if(var_est=="jack")	
		{
		N.hat.strata.jack = matrix(NA, num_strata, n0)		
		N.hat.strata.var = numeric()	
		N.hat.strata.1.jack = matrix(NA, num_strata, n0)		
		N.hat.strata.1.var = numeric()				
		jj = 0
		for(i in s0)
			{	
			jj = jj+1
		
			#The stratum membership and index location of this element
			stratum_i = pop$strata[i]
			i_index = which(names(r.1.records[[stratum_i]]) == i)

			#N.hat.strata estimators with this element removed from the initial sample
			r.jack = r[stratum_i] - r.1.records[[stratum_i]][i_index] - sum(pop$x[i,s0_strata[[stratum_i]]])
			l.jack = l[stratum_i] - sum(pop$x[i,strata_setdiff_s0[[stratum_i]]]) + sum(pop$x[setdiff(s0_strata[[stratum_i]],i),i])
			N.hat.strata.jack[stratum_i,jj] = initial.size.strata[stratum_i] * (r.jack + l.jack + 1) / (r.jack + 1) - 1

			r.1.jack = length(which((r.1.records[[stratum_i]][-i_index] - pop$x[i,setdiff(s0_strata[[stratum_i]],i)]) > 0))
			l.1.jack = length(which((l.1.records[[stratum_i]] - pop$x[i,strata_setdiff_s0[[stratum_i]]]) > 0)) + 
							(sum(pop$x[setdiff(s0_strata[[stratum_i]],i),i]) > 0)
			N.hat.strata.1.jack[stratum_i,jj] = initial.size.strata[stratum_i] * (r.1.jack + l.1.jack + 1) / (r.1.jack + 1) - 1		
			}
					
		for(k in 1:num_strata)
			{
			N.hat.strata.jack.k = na.omit(N.hat.strata.jack[k,])
			N.hat.strata.jack.mean.k = mean(N.hat.strata.jack.k)
			N.hat.strata.var[k] = (initial.size.strata[k]-2)/(2*initial.size.strata[k]) * sum((N.hat.strata.jack.k - N.hat.strata.jack.mean.k)^2)
			
			N.hat.strata.1.jack.k = na.omit(N.hat.strata.1.jack[k,])
			N.hat.strata.1.jack.mean.k = mean(N.hat.strata.1.jack.k)
			N.hat.strata.1.var[k] = (initial.size.strata[k]-2)/(2*initial.size.strata[k]) * sum((N.hat.strata.1.jack.k - N.hat.strata.1.jack.mean.k)^2)			
			}
		}
				
	#Make the estimates available outside the function
	N.hat.strata <<- N.hat.strata
	N.hat.strata.var <<- N.hat.strata.var
	
	N.hat.strata.1 <<- N.hat.strata.1
	N.hat.strata.1.var <<- N.hat.strata.1.var
	}
