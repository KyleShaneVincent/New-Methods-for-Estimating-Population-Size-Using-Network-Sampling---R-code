#The "Load_P90_data" function will load the P90 project data; for more information on this data set, see here: https://oprdata.princeton.edu/Archive/P90/

Load_P90_data <- function(CSpop, num_nominations, plot.P90.data)
	{	
	##The population of size 595	
	if(CSpop == 595)
		{
		#Read in the nodes data set (this gives observations for each respondent's covariates)
		P90_data_nodes = read.table("yvalues")
		P90_data_nodes = P90_data_nodes[,2]
		
		#View the dimensions of the nodes data set followed by the first few lines
		cat("The P90 nodes data set:",'\n')
		print(P90_data_nodes)
		cat('\n')
		
		#Review the summary of the nodes data set
		cat("A summary of the P90 nodes data set:", '\n')
		print(table(P90_data_nodes, useNA="always"))
		cat('\n')
		
		#Read in the edges data set (this gives the information required to construct the population adjacency matrix)
		P90_data_edges = as.matrix(read.table("xdrugcs.txt", header=FALSE))

		#View the dimensions of the edges data set followed by the first few lines
		cat("The size of the P90 edges data set with the first few entries:", '\n')
		print(dim(P90_data_edges))
		print(P90_data_edges[1:5,1:5])	
		cat('\n')
		
		#Construct the population adjacency matrix
		N = length(P90_data_nodes) #The population size
		cat("The population size is:", '\n')
		print(N)
		cat('\n')
		adjacency = P90_data_edges
		rownames(adjacency) = colnames(adjacency) = 1:595
				
		#The total number of links in the networked population 
		cat("The dimensions of the adjacency data set:", '\n')
		print(dim(adjacency))
		cat('\n')
		cat("The total number of links in the population network:", '\n')
		print(sum(adjacency))
		cat('\n')
			
		#The list of information required for the simulation study
		pop = list()
		pop$U = 1:N
		pop$x = adjacency
		rownames(pop$x) = colnames(pop$x) = 1:N

		#Plot the population network with the aid of the 'igraph' package
		if(plot.P90.data==TRUE)
			{
			#This will be based on drug-use status 
			cat("Plotting...", '\n')
			vertex.color = numeric(N)
			vertex.color[which(P90_data_nodes==0)] = "yellow" 	#Non-injection drug-user
			vertex.color[which(P90_data_nodes==1)] = "red" 		#Injection drug-user
			vertex.size = rep(2,N)
			G = graph.adjacency(adjacency)
			loc = layout_components(G)
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Population Graph",
			xlab="Non-injection drug users as yellow circles, Injection drug users as red circles")
			}
			
		pop$strata = NA
		x = which(P90_data_nodes==0) 	#Non-injection drug-user
		pop$strata[x] = 1 
		x = which(P90_data_nodes==1) 	#Injection drug-user
		pop$strata[x] = 2

		num_strata = max(pop$strata)

		cat("The size of each stratum is ", table(pop$strata), '\n')

		pop.strata = list()
		for(k in 1:num_strata)
			pop.strata[[k]] = which(pop$strata==k)	
			
		pop$sex.worker = rep(NA,595)
		pop.sex.worker = rep(NA,595)
		
		#Redefine the adjacency matrix based on num_nominations; priority is given to drug-use status
		if(num_nominations != Inf)
			{
			set.seed(1)
			for(i in 1:595)
				{
				noms_i = which(pop$x[i,] == 1)
				if(length(noms_i) > num_nominations)	
					{
					k = 0
					prob_noms_i = numeric()
					for(j in noms_i)
						{
						k = k+1
						prob_noms_i[k] = sum((P90_data_nodes[k] == P90_data_nodes[j])+1)
						}
					pop$x[i,] = rep(0,595)
					pop$x[i,sample(noms_i, num_nominations, replace=FALSE, prob=prob_noms_i)] = 1
					}
				}
				
			adjacency = pop$x	
			
			cat("Based on a maximum number of nominations equal to", num_nominations, "the total number of links in the population network is:", '\n')
			print(sum(adjacency))
			cat('\n')
						
			#Plot the restricted population network with the aid of the 'igraph' package
			if(plot.P90.data==TRUE)
				{
				X11()
				G = graph.adjacency(adjacency)
				loc = layout_components(G)
				
				plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
				edge.width=0.2,edge.color="black",edge.arrow.mode=0, 
				main=paste("Population Graph Based on Restricted Number of Nominations of ", num_nominations, sep = " "),
							xlab="Non-injection drug users as yellow circles, Injection drug users as red circles")
				}			
			}			
			
		#The network size
		pop$net_size = apply(pop$x,1,sum)
		names(pop$net_size) = NULL	
		
		adjacency <<- adjacency
		pop <<- pop	
		pop.strata <<- pop.strata
		num_strata <<- num_strata
		pop.sex.worker <<- pop.sex.worker
		P90_data_nodes <<- P90_data_nodes
		P90_data_edges <<- P90_data_edges
		}
	
	##The population of size 5492
	if(CSpop == 5492)
		{
		#Read in the nodes data set (this gives observations for each respondent's covariates)
		P90_data_nodes = read.table("nodes.tsv", header=TRUE)

		#View the dimensions of the nodes data set followed by the first few lines
		cat("The size of the P90 nodes data set with the first few lines:",'\n')
		print(dim(P90_data_nodes))
		print(head(P90_data_nodes))
		cat('\n')

		#Notice that there are some missing entries
		cat("A summary of the P90 nodes data set (notice that there are some missing entries that need to be imputed): ", '\n')
		print(summary(P90_data_nodes))
		cat('\n')

		#Impute the missing entries with the aid of the 'mice' package so that they can be used in the simulation study
		cat("Imputing missing entries:")
		P90_data_nodes.impute <- mice(P90_data_nodes, m=1, maxit=1, seed=1)
		cat('\n')
		
		#The "complete" function will return a data frame with all values and we will call this P90_data_nodes
		P90_data_nodes <- complete(P90_data_nodes.impute)

		#Review the summary of the nodes data set
		cat("A summary of the P90 nodes data set after imputing the missing entries:", '\n')
		print(summary(P90_data_nodes))
		cat('\n')
		
		#Read in the edges data set (this gives the information required to construct the population adjacency matrix)
		P90_data_edges = read.table("edges.tsv",header=TRUE)

		#View the dimensions of the edges data set followed by the first few lines
		cat("The size of the P90 edges data set with the first few lines:", '\n')
		print(dim(P90_data_edges))
		print(head(P90_data_edges))	
		cat('\n')

		#Construct the population adjacency matrix
		N = nrow(P90_data_nodes) #The population size
		cat("The population size is:", '\n')
		print(N)
		cat('\n')
		adjacency = matrix(0,N,N)
		adjacency[cbind(P90_data_edges[,1],P90_data_edges[,2])] = 1
			
		#The total number of links in the networked population 
		cat("The dimensions of the adjacency data set:", '\n')
		print(dim(adjacency))
		cat('\n')
		cat("The total number of links in the population network:", '\n')
		print(sum(adjacency))
		cat('\n')
			
		#The list of information requried for the simulation study
		pop = list()
		pop$U = 1:N
		pop$x = adjacency
		rownames(pop$x) = colnames(pop$x) = 1:N

		#Plot the population network with the aid of the 'igraph' package
		if(plot.P90.data==TRUE)
			{
			#This will be based on the race (white and non-white) and gender (male and female) of the individuals
			#Race: 1 = Native American, 2 = Black, 3 = Asian/Pacific Islander, 4 = White, 5 = Other
			#Gender: 0 = Male, 1 = Female
			cat("Plotting full network graph...", '\n', '\n')
			vertex.color = numeric(N)
			vertex.color[which(P90_data_nodes$race==4)] = "yellow" #White
			vertex.color[which(P90_data_nodes$race!=4)] = "red" #Non-white
			vertex.shape = numeric(N)
			vertex.shape[which(P90_data_nodes$gender==0)] = "square" #Males
			vertex.shape[which(P90_data_nodes$gender==1)] = "circle" #Females
			vertex.size = rep(2,N)
			G = graph.adjacency(adjacency)
			loc = layout_components(G)
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			vertex.shape=vertex.shape,edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Population Graph",
			xlab="White Males as yellow squares, White Females as yellow circles, Non-White Males as red squares, Non-White Females as red circles")
			}
			
		pop$strata = NA
		x = which(P90_data_nodes$race==4 & P90_data_nodes$gender==0) #White males
		pop$strata[x] = 1 
		x = which(P90_data_nodes$race==4 & P90_data_nodes$gender==1) #White females
		pop$strata[x] = 2
		x = which(P90_data_nodes$race!=4 & P90_data_nodes$gender==0) #Non-white males
		pop$strata[x] = 3 
		x = which(P90_data_nodes$race!=4 & P90_data_nodes$gender==1) #Non-white females
		pop$strata[x] = 4

		num_strata = max(pop$strata)

		cat("The size of each stratum is ", table(pop$strata), '\n')

		pop.strata = list()
		for(k in 1:num_strata)
			pop.strata[[k]] = which(pop$strata==k)	
			
			
		pop$sex.worker = rep(NA,5492)
		x = which(P90_data_nodes$sex.worker==0) #Non-sex workers
		pop$sex.worker[x] = 0 	
		x = which(P90_data_nodes$sex.worker==1) #Sex workers
		pop$sex.worker[x] = 1 	
		
		cat("The distribution of sex workers is ", table(pop$sex.worker), '\n', '\n', '\n')		
		
		pop.sex.worker = list()
		pop.sex.worker[[1]] = which(pop$sex.worker==1)
		pop.sex.worker[[2]] = which(pop$sex.worker==0)		
		
		#Redefine the adjacency matrix based on num_nominations; priority is given to race and gender
		if(num_nominations != Inf)
			{
			set.seed(1)
			for(i in 1:5492)
				{
				noms_i = which(pop$x[i,] == 1)
				if(length(noms_i) > num_nominations)	
					{
					k = 0
					prob_noms_i = numeric()
					for(j in noms_i)
						{
						k = k+1
						prob_noms_i[k] = sum((P90_data_nodes[k,] == P90_data_nodes[j,])*c(0,rep(1,2),rep(1/2,11)))
						}
					pop$x[i,] = rep(0,5492)
					pop$x[i,sample(noms_i, num_nominations, replace=FALSE, prob=prob_noms_i)] = 1
					}
				}
				
			adjacency = pop$x	
			
			cat("Based on a maximum number of nominations equal to", num_nominations, "the total number of links in the population network is:", '\n')
			print(sum(adjacency))
			cat('\n')
					
			#Plot the restricted population network with the aid of the 'igraph' package
			if(plot.P90.data==TRUE)
				{
				X11()
				G = graph.adjacency(adjacency)
				loc = layout_components(G)
				plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
				vertex.shape=vertex.shape,edge.width=0.2,edge.color="black",edge.arrow.mode=0, 
				main=paste("Population Graph Based on Restricted Number of Nominations of ", num_nominations, sep = " "),
				xlab="White Males as yellow squares, White Females as yellow circles, Non-White Males as red squares, Non-White Females as red circles")
				}			
			}	
			
		#The network size
		pop$net_size = apply(pop$x,1,sum)
		names(pop$net_size) = NULL	
			
		adjacency <<- adjacency
		pop <<- pop	
		pop.strata <<- pop.strata
		num_strata <<- num_strata
		pop.sex.worker <<- pop.sex.worker
		P90_data_nodes <<- P90_data_nodes
		P90_data_edges <<- P90_data_edges
		}		
	}
	