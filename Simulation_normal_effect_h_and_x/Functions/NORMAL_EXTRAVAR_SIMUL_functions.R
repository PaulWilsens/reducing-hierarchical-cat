
tree_aggregation = function(embeddings,tree_levels,embed_dim){
  
  embedding_list = list()
  
  assign(paste0("embeddings_lev",tree_levels),embeddings[,1:embed_dim])
  
  embedding_list[[tree_levels]] = get(paste0("embeddings_lev",tree_levels))
  
  #####
  
  for (i in 1:(tree_levels-1)) {
    
    
    temp = get(paste0("embeddings_lev",tree_levels-i+1))
    temp_sizes = get(paste0("subgroupsizes",tree_levels-i+1))
    
    inter = matrix(nrow=length(temp_sizes),ncol=embed_dim)
    
    k=0
    for (j in 1:(length(temp_sizes))) {
      
      inter[j,] = colMeans(  temp[(k+1):(k+temp_sizes[j]),]   )
      
      k = k + temp_sizes[j]
    }
    
    assign(paste0("embeddings_lev",tree_levels-i),inter)
    
    append(embedding_list,get(paste0("embeddings_lev",tree_levels-i)))
    
    embedding_list[[tree_levels-i]] = get(paste0("embeddings_lev",tree_levels-i))
  }
  
  #####
  
  
  return(embedding_list)
}


all_children = function(dataframe,z) {
  
  # z=0: looking at child nodes of a cluster on level j-1
  # z=1: looking at child nodes of a cluster on level j
  
  child_set = data.frame()
  
  #getting subgroup sizes of level of child nodes
  subgroupsizes_temp = get(paste0("subgroupsizes",j+z))
  
  
  
  for (i in 1:(nrow(dataframe))) {   
    
    #index matching which element it was in original structure on relevant level
    emb_num = which(apply(embedding_list[[j-1+z]], 1, function(x) return(all(x == as.vector(dataframe[i,])))))    
    
    
    if (emb_num == 1) {
      child_set = rbind(child_set, embedding_list[[j+z]][1:(cumsum(subgroupsizes_temp)[emb_num]),])
    } else {
      child_set = rbind(child_set, embedding_list[[j+z]][(cumsum(subgroupsizes_temp)[emb_num-1]+1):(cumsum(subgroupsizes_temp)[emb_num]),])
    }
    
    
  }
  
  
  return(child_set)
}


horizontal_clustering = function(level,cutoff=0.7) {
  
  j = level
  
  if (level != 1) {
    
    
    
    G_list_temp = list()
    
    fake_temp = data.frame()
    if (length(get(paste0("fake_lev",j))) > 0) {
      for (i in 1:(length(get(paste0("fake_lev",j)))) ) {
        fake_temp = rbind(fake_temp, as.data.frame(get(paste0("fake_lev",j))[[i]][2]))
      }
    }
    
    
    #####
    
    
    if (length(get(paste0("clusters_lev",j-1))) > 0) {
      
      for (l in 1:(length(get(paste0("clusters_lev",j-1)))) ) {
        
        temp = all_children(get(paste0("clusters_lev",j-1))[[l]],z=0)
        
        
        temp2 = data.frame() 
        for (i in 1:(nrow(temp))) {
          
          if(sum(apply(fake_temp, 1, function(x) return(all(x == temp[i,] ))),na.rm=TRUE) == 0){
            temp2 = rbind(temp2,temp[i,])
          }
          
        }
        
        
        if(length(temp2) != 0){
          G_list_temp[[length(G_list_temp) + 1]] = temp2
        }
        
      }  
      
    }
    
    
    
    if (length(get(paste0("fake_lev",j-1))) != 0) {
      
      for (l in 1:(length(get(paste0("fake_lev",j-1))))) {
        
        temp = all_children(get(paste0("fake_lev",j-1))[[l]][[2]],z=0)
        
        temp2 = data.frame() 
        for (i in 1:(nrow(temp))) {
          
          if(sum(apply(fake_temp, 1, function(x) return(all(x == temp[i,] ))),na.rm=TRUE) == 0){
            temp2 = rbind(temp2,temp[i,])
          }
          
        }
        
        if(length(temp2) != 0){
          G_list_temp[[length(G_list_temp) + 1]] = temp2
        }
        
      }
      
      
    }
    #####
    
    
    assign(paste0("G_lev",j),G_list_temp)    
    
  } 
  
  
  
  
  #################################
  ##   CLUSTERING 
  #################################
  
  
  temp_list = list()
  
  if (length(get(paste0("G_lev",j))) > 0) {
    
    for (l in 1:(length(get(paste0("G_lev",j))))  ) {
      
      G_temp = get(paste0("G_lev",j))[[l]]
      
      if (nrow(G_temp)  >= 3) {
        
        silhouette_temp = fviz_nbclust(G_temp, clara, method='silhouette', k.max = nrow(G_temp)-1)
        
        #
        if (max(silhouette_temp$data$y) >= cutoff) {
          n_clust = which(silhouette_temp$data$y == max(silhouette_temp$data$y))
          
          clust_object = clara(G_temp,n_clust)
          
          clust_temp = split(as.data.frame(G_temp),clust_object$clustering)
          
          ##
          for (i in 1:(length(clust_temp))  ) {
            
            temp_list[[length(temp_list) + 1]] = as.matrix(clust_temp[[i]])
            
          }
          ##
        } else {
          temp_list[[length(temp_list) + 1]] = as.matrix(G_temp)
        }
        #
        
        
      } else {
        temp_list[[length(temp_list) + 1]] = as.matrix(G_temp)
      }
      
      
      
    }
    
  }
  
  
  
  #################################
  
  return(temp_list)
  
}



vertical_clustering = function(level,cutoff=0.7) {
  
  
  j = level
  
  
  if (j < length(embedding_list)) {
    
    
    temp_list = list()
    
    if (length(get(paste0("clusters_lev",j))) != 0) {
      
      ##########
      
      for (l in 1:(length(get(paste0("clusters_lev",j))))    ) {
        
        temp  = colMeans(get(paste0("clusters_lev",j))[[l]])
        temp2 = all_children(get(paste0("clusters_lev",j))[[l]],z=1) 
        
        temp = rbind(temp,temp2) 
        
        silhouette_temp = fviz_nbclust(temp, clara, method='silhouette', k.max = nrow(temp)-1)
        
        #
        if (max(silhouette_temp$data$y) >= cutoff) {
          n_clust = which(silhouette_temp$data$y == max(silhouette_temp$data$y))
          
          clust_object = clara(temp,n_clust)
          
          id_temp = clust_object$clustering[1]
          
          ##
          
          
          if (sum(clust_object$clustering == id_temp) > 1) {
            
            clust_temp = split(as.data.frame(temp),clust_object$clustering)[[id_temp]]
            
            
            silhouette_temp2 = silhouette(clust_object$clustering,dist(temp))[,3][silhouette(clust_object$clustering,dist(temp))[,1]==id_temp]
            
            
            
            
            ##
            if (min(silhouette_temp2) < silhouette_temp2[1]) {
              
              f_element_temp = list()
              
              f_element_temp[[1]] = as.matrix(clust_temp[1,])   
              f_element_temp[[2]] = as.matrix(clust_temp[2:nrow(clust_temp),])
              
              temp_list[[length(temp_list) + 1]] = f_element_temp
              
            } 
            
            
            ##
            
          }
          ##
          
          
        } else {
          
          
          
          f_element_temp = list()
          
          f_element_temp[[1]] = as.matrix(temp[1,])   
          f_element_temp[[2]] = as.matrix(temp[2:nrow(temp),])
          
          temp_list[[length(temp_list) + 1]] = f_element_temp
          
        }
        #
        
      }
      
      ##########
      
    }
    
    
    ###################
    ###################
    
    
    if (length(get(paste0("fake_lev",j))) != 0) {
      
      
      for (l in 1:(length(get(paste0("fake_lev",j))))    ) {
        
        temp  = as.matrix(get(paste0("fake_lev",j))[[l]][[1]])
        colnames(temp) = c(paste0("X",1:embed_dim))
        temp2 = as.matrix(all_children(get(paste0("fake_lev",j))[[l]][[2]],z=1))
        
        temp = rbind(temp,temp2) 
        
        silhouette_temp = fviz_nbclust(temp, clara, method='silhouette', k.max = nrow(temp)-1)
        
        #
        if (max(silhouette_temp$data$y) >= cutoff) {
          n_clust = which(silhouette_temp$data$y == max(silhouette_temp$data$y))
          
          clust_object = clara(temp,n_clust)
          
          id_temp = clust_object$clustering[1]
          
          ##
          
          if (sum(clust_object$clustering == id_temp) > 1) {
            
            clust_temp = split(as.data.frame(temp),clust_object$clustering)[[id_temp]]
            
            
            silhouette_temp2 = silhouette(clust_object$clustering,dist(temp))[,3][silhouette(clust_object$clustering,dist(temp))[,1]==id_temp]
            
            
            
            
            ##
            if (min(silhouette_temp2) < silhouette_temp2[1]) {
              
              f_element_temp = list()
              
              f_element_temp[[1]] = as.matrix(clust_temp[1,])   
              f_element_temp[[2]] = as.matrix(clust_temp[2:nrow(clust_temp),])
              
              temp_list[[length(temp_list) + 1]] = f_element_temp
              
            } 
            
            
            ##
            
            
            
            
            
          }
          ##
          
          
        } else {
          
          f_element_temp = list()
          
          f_element_temp[[1]] = matrix(temp[1,],ncol=embed_dim)   
          f_element_temp[[2]] = as.matrix(temp[2:nrow(temp),])
          
          temp_list[[length(temp_list) + 1]] = f_element_temp
          
        }
        #
        
      }  
      
    }
    
    
    
    
  } else {
    temp_list = list()
  }
  
  
  
  return(temp_list)
  
}


condense_result = function(listH,listV) {
  
  tempH = list()
  tempV = list()
  
  
  #condense HORIZONTAL results ####
  for (l in 1:length(listH)) {
    
    temp_list = list()
    
    if (length(listH[[l]]) > 0) {
      
      for (i in 1:length(listH[[l]])) {
        
        temp_vec = c()
        
        for (k in 1:(nrow(listH[[l]][[i]]))   ) {
          
          emb_temp = which(apply(embedding_list[[l]], 1, function(x) return(all(x == as.vector(listH[[l]][[i]][k,])))))
          temp_vec = c(temp_vec,emb_temp)
        }
        
        temp_list[[length(temp_list)+1]] = temp_vec
        
      }
      
    }
    
    tempH[[length(tempH)+1]] = temp_list
  }
  ####
  
  #condense VERTICAL results ####
  for (l in 1:length(listV)) {
    
    temp_list = list()
    
    if (length(listV[[l]]) > 0) {
      
      for (i in 1:length(listV[[l]])) {
        
        temp_vec = c()
        
        for (k in 1:(nrow(listV[[l]][[i]][[2]]))   ) {
          
          emb_temp = which(apply(embedding_list[[l]], 1, function(x) return(all(x == as.vector(listV[[l]][[i]][[2]][k,])))))
          temp_vec = c(temp_vec,emb_temp)
        }
        
        temp_list[[length(temp_list)+1]] = temp_vec
        
      }
      
    }
    
    tempV[[length(tempV)+1]] = temp_list
  }
  ####
  
  
  ##
  temp = list(tempH,tempV)
  return(temp)
  
}

