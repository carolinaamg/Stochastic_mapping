rm(list = ls())

library(ape)
library(phytools)
library(tidytree)
library(hash)

dat <- read.table(file="/home/carolina/Documents/stoch_mapp/matrix_example.csv", sep="\t", header=T, check.names=F, row.names = 1)
tree <- read.tree("/home/carolina/Documents/stoch_mapp/phylogeny_example.nwk")
clusters <- colnames(dat); clusters
replicates <- 100 ###CHANGE DEPENDING ON THE NUMBER OF REPLICATES

#####################################################################################
############################# Stochastic mapping ####################################
#####################################################################################
total_spp <- length(tree$tip.label); total_spp
total_nodes <- length(tree$node.label); total_nodes
dat[dat>1] <- 1
final <- colSums(dat) < total_spp & colSums(dat) > 1
dat <- dat[,final]

## Standard stochastic mapping
all_stochastic <- list()
all_simulation <- list()
for(i in 1:dim(dat)[2]) {
    cluster <- clusters[[i]]
    print(cluster)

    ###Getting cluster p/a dataframe###
    orth2 <-setNames(dat[,i],rownames(dat)); orth2

    ###Getting stochastic mapping###
    stoch_map <- data.frame(matrix(nrow = total_nodes, ncol =replicates)) ###CHANGE
    colnames(stoch_map) <- 1:replicates ###CHANGE
    row.names(stoch_map) <- tree$node.label
    smap <- make.simmap(tree, orth2, model="ARD", nsim=replicates) #Stochastic mapping ###CHANGE

    ###Checking how good Q matrices are
    if(smap[[1]]$Q[1,1] == 0){
      print("Changing Q matrix first row...")
      data <- c(-0.00001, 0.00001, as.numeric(smap[[1]]$Q[2,1]), as.numeric(smap[[1]]$Q[2,2]))
      Q = matrix(data, nrow=2, ncol = 2, byrow = TRUE)
      print(Q)
    } else if (smap[[1]]$Q[2,1] == 0){
      print("Changing Q matrix second row...")
      data <- c(as.numeric(smap[[1]]$Q[1,1]), as.numeric(smap[[1]]$Q[1,2]), 0.00001, -0.00001)
      Q = matrix(data, nrow=2, ncol = 2, byrow = TRUE)
      print(Q)
    } else {
      print("Keeping Q matrix as it is...")
      Q = smap[[1]]$Q
      print(Q)
    }

    for(n in 1:replicates){###CHANGE
      summ = summary(smap[[n]])
      stoch_map[,n] <- round(as.numeric(summ$states))
    }

    ###Getting simulations###
    simulations <- data.frame(matrix(nrow = total_nodes, ncol =replicates)) ###CHANGE
    colnames(simulations) <- 1:replicates ###CHANGE
    row.names(simulations) <- tree$node.label
    nulldist <- sim.history(tree, Q, nsim=replicates)#Simulations ###CHANGE

    for(x in 1:replicates){
      summ = summary(nulldist[[x]])
      simulations[,x] <- round(as.numeric(summ$states))
    }
    all_stochastic[[i]] <- stoch_map
    all_simulation[[i]] <- simulations
}

#####################################################
###Getting losses and gains for STOCHASTIC MAPPING###
nodes_clusters_gained <- data.frame(matrix(ncol=dim(stoch_map)[2], nrow=total_nodes-1, "NA"))
nodes_clusters_lost <- data.frame(matrix(ncol=dim(stoch_map)[2], nrow=total_nodes-1, "NA"))
for(n in 1:dim(stoch_map)[2]){
    ###Getting an empty dataframe
    full_anc <- data.frame(matrix(ncol=dim(dat)[2], nrow=total_nodes, "NA"))
    row.names(full_anc) <- tree$node.label
    colnames(full_anc) <- colnames(dat)

    for(i in 1:dim(dat)[2]){
      cluster = clusters[i]
      print(cluster)
      matrix_cluster = all_stochastic[[i]]
      full_anc[,i] <- matrix_cluster[,n]
    }

    ###Finding losses and gains in the tree
    x_tree <- as_tibble(tree)
    x_prob <- full_anc
    f_table <- rbind(x_prob, dat, colnames(x_prob))

    #Get full table with internal nodes and leaves
    lg <- list()
    ll <- list()
    gained_clusters <- hash()
    lost_clusters <- hash()
    nodes <- list()
    all_nodes <- tree$node.label

    ##Looping across the nodes
    for(q in 2:length(all_nodes)){
      node <- all_nodes[q]
      p <- parent(x_tree, node)
      parent <- p$label
      gained <- colnames(x_prob)[f_table[parent,] == 0 & f_table[node,] == 1] # genes gained
      lost   <- colnames(x_prob)[f_table[parent,] == 1 & f_table[node,] == 0] # genes lost
      lg[q] <- length(gained)
      ll[q] <- length(lost)
      nodes[q] <- node
      gained_clusters[node] <- list(gained)
      lost_clusters[node] <- list(lost)
    }
    ###Getting a dataframe for gains
    nodes_clusters_gained[,n] <- as.character(values(gained_clusters))
    nodes_clusters_lost[,n] <- as.character(values(lost_clusters))
}

row.names(nodes_clusters_gained) <- keys(gained_clusters)
row.names(nodes_clusters_lost) <- keys(lost_clusters)
write.table(nodes_clusters_gained, file='gains_stochastic_SM.tsv', sep = "\t", quote=F, row.names=T)
write.table(nodes_clusters_lost, file='losses_stochastic_SM.tsv', sep = "\t", quote=F, row.names=T)

##########################################
###Getting losses and gains SIMULATIONS###
nodes_clusters_gained_sim <- data.frame(matrix(ncol=dim(simulations)[2], nrow=total_nodes-1, "NA"))
nodes_clusters_lost_sim <- data.frame(matrix(ncol=dim(simulations)[2], nrow=total_nodes-1, "NA"))
for(n in 1:dim(simulations)[2]){
    ###Getting an empty dataframe
    full_anc <- data.frame(matrix(ncol=dim(dat)[2], nrow=total_nodes, "NA"))
    row.names(full_anc) <- tree$node.label
    colnames(full_anc) <- colnames(dat)

    for(i in 1:dim(dat)[2]){
      cluster = clusters[i]
      print(cluster)
      matrix_cluster = all_simulation[[i]]
      full_anc[,i] <- matrix_cluster[,n]
    }

    ###Finding losses and gains in the tree
    x_tree <- as_tibble(tree)
    x_prob <- full_anc
    f_table <- rbind(x_prob, dat, colnames(x_prob))

    #Get full table with internal nodes and leaves
    lg <- list()
    ll <- list()
    gained_clusters <- hash()
    lost_clusters <- hash()
    nodes <- list()
    all_nodes <- tree$node.label

    ##Looping across the nodes
    for(q in 2:length(all_nodes)){
      node <- all_nodes[q]
      p <- parent(x_tree, node)
      parent <- p$label
      gained <- colnames(x_prob)[f_table[parent,] == 0 & f_table[node,] == 1] # genes gained
      lost   <- colnames(x_prob)[f_table[parent,] == 1 & f_table[node,] == 0] # genes lost
      lg[q] <- length(gained)
      ll[q] <- length(lost)
      nodes[q] <- node
      gained_clusters[node] <- list(gained)
      lost_clusters[node] <- list(lost)
    }
    ###Getting a dataframe for gains
    nodes_clusters_gained_sim[,n] <- as.character(values(gained_clusters))
    nodes_clusters_lost_sim[,n] <- as.character(values(lost_clusters))
}

row.names(nodes_clusters_gained_sim) <- keys(gained_clusters)
row.names(nodes_clusters_lost_sim) <- keys(lost_clusters)
write.table(nodes_clusters_gained_sim, file='gains_simulation_SM.tsv', sep = "\t", quote=F, row.names=T)
write.table(nodes_clusters_lost_sim, file='losses_simulation_SM.tsv', sep = "\t", quote=F, row.names=T)
