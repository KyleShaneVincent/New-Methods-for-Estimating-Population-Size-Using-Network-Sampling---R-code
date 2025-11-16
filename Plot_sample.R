#The "Plot_sample" function will first plot the network sample based on the traced links to new members and then plot the fully observed network sample.

Plot_sample <- function(CSpop, plot.sample)
	{
	if(plot.sample==TRUE)
		{
		if(CSpop==595)
			{
			par(mfrow=c(2,2))
			#Plot the sample		
			s.plot = unlist(s); names(s.plot) = NULL
			n.plot = length(s.plot)
			vertex.label = s.plot
			vertex.color = numeric(n.plot)
			vertex.color[which(pop$strata[s.plot] == 1)] = "yellow" #Non-injection drug user
			vertex.color[which(pop$strata[s.plot] == 2)] = "red"    #Injection drug user
			vertex.size = rep(3,n.plot)
			vertex.size[which(s.plot %in% s$Wave_0)] = 4.5

			#Plot the graph based on the traced links to new members
			adjacency.plot = adjacency.plot[s.plot,s.plot]
			G = graph.adjacency(adjacency.plot)
			loc = layout_components(G)
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Sample Graph with Seeds Represented as Englarged Nodes\n
			Only Traced Links Depicted in Graph",
			xlab="Non-injection drug users as yellow circles, Injection drug users as red circles")

			#Plot the full sample network graph
			G = graph.adjacency(adjacency[s.plot,s.plot])
			loc = loc
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Full Network Sample Graph with Seeds Represented as Englarged Nodes",
			xlab="Non-injection drug users as yellow circles, Injection drug users as red circles")

			#Replot the full sample network graph with organized coordinates
			loc = layout_components(G)
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Reorganized Full Network Sample Graph with Seeds Represented as Englarged Nodes",
			xlab="Non-injection drug users as yellow circles, Injection drug users as red circles")
			}	
			
		if(CSpop==5492)
			{
			par(mfrow=c(2,2))
			#Plot the sample		
			s.plot = unlist(s); names(s.plot) = NULL
			n.plot = length(s.plot)
			vertex.label = s.plot
			vertex.color = numeric(n.plot)
			vertex.color[which(pop$strata[s.plot] %in% c(1,2))] = "yellow" #White
			vertex.color[which(pop$strata[s.plot] %in% c(3,4))] = "red" #Non-white
			vertex.shape = numeric(n.plot)
			vertex.shape[which(pop$strata[s.plot] %in% c(1,3))]= "square" #Males
			vertex.shape[which(pop$strata[s.plot] %in% c(2,4))]= "circle" #Females
			vertex.size = rep(3,n.plot)
			vertex.size[which(s.plot %in% s$Wave_0)] = 4.5

			#Plot the graph based on the traced links to new members
			adjacency.plot = adjacency.plot[s.plot,s.plot]
			G = graph.adjacency(adjacency.plot)
			loc = layout_components(G)
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			vertex.shape=vertex.shape,edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Sample Graph with Seeds Represented as Englarged Nodes\n
			Only Traced Links Depicted in Graph",
			xlab="White Males as yellow squares, White Females as yellow circles, Non-White Males as red squares, Non-White Females as red circles")

			#Plot the full sample network graph
			G = graph.adjacency(adjacency[s.plot,s.plot])
			loc = loc
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			vertex.shape=vertex.shape,edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Full Network Sample Graph with Seeds Represented as Englarged Nodes",
			xlab="White Males as yellow squares, White Females as yellow circles, Non-White Males as red squares, Non-White Females as red circles")

			#Replot the full sample network graph with organized coordinates
			loc = layout_components(G)
			plot.igraph(G,layout=loc,vertex.size=vertex.size,vertex.label=NA,vertex.color=vertex.color,
			vertex.shape=vertex.shape,edge.width=0.2,edge.color="black",edge.arrow.mode=0, main="Reorganized Full Network Sample Graph with Seeds Represented as Englarged Nodes",
			xlab="White Males as yellow squares, White Females as yellow circles, Non-White Males as red squares, Non-White Females as red circles")
			}	
		}			
	}
