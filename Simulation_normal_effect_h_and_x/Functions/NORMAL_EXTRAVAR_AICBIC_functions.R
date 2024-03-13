designfactor = function(final){
  
  for (r in 0:(tree_levels-2) ) {
    
    k = tree_levels-r
    
    subgroup_temp = get(paste0("subgroupsizes",k))
    temp = c()
    
    for (i in 1:length(subgroup_temp)) {
      
      temp = c(temp,rep(i,subgroup_temp[i]))
      
    }
    
    assign(paste0("parents_",k),temp)
    
  }
  
  
  #################################################
  
  parent_cluster_list = list()
  
  for (r in 1:length(final[[1]])) {
    
    temp_parent_list = list()
    
    if(r>1) {
      
      if(length(final[[1]][[r]]) > 0) {
        for (s in 1:length(final[[1]][[r]])) {
          
          temp_class = final[[1]][[r]][[s]][1]
          temp_class_parent = get(paste0("parents_",r))[temp_class]
          r_parent = r-1
          indic = 0
          
          
          while (r_parent>1 & indic==0) { 
            
            temp_list = final[[2]][[r_parent]]
            
            ##
            if(any(sapply(temp_list, function(vec) temp_class_parent %in% vec))) {
              
              temp_class_parent = get(paste0("parents_",r_parent))[temp_class_parent]
              r_parent = r_parent-1
              
            } else {
              indic = indic + 1
            }
            ##
            
          }
          
          
          
          s_parent = which(sapply(final[[1]][[r_parent]], function(vec) temp_class_parent %in% vec))
          
          temp_parent_list[[length(temp_parent_list)+1]] = c(r_parent,s_parent)
          
        }
      }
      
    }
    
    
    parent_cluster_list[[length(parent_cluster_list)+1]] = temp_parent_list
    
    
  }
  
  #################################################
  
  subgroupsizes_collapsed = list()
  
  for (r in 1:length(final[[1]])) {
    
    if(r==1) {
      
      subgroupsizes_collapsed[[length(subgroupsizes_collapsed)+1]] = length(final[[1]][[r]])
      
    } else {
      
      temp_list = parent_cluster_list[[r]]
      count_list= list()
      
      for (vec in temp_list) {
        vec_str = paste(vec, collapse = ",")
        
        if (vec_str %in% names(count_list)) {
          count_list[[vec_str]] = count_list[[vec_str]] + 1
        } else {
          count_list[[vec_str]] = 1
        }
      }
      
      if(length(count_list) == 0) {
        count_list[[1]] = 0
      }
      
      subgroupsizes_collapsed[[length(subgroupsizes_collapsed)+1]] = as.numeric(unlist(count_list))
      
    }
    
  }
  
  
  
  n_indic = 0
  n_total = 0
  for (r in 1:length(subgroupsizes_collapsed)) {
    
    if(length(subgroupsizes_collapsed[[r]]) > 0) {
      n_total = n_total + sum(subgroupsizes_collapsed[[r]])
      
      for (s in 1:length(subgroupsizes_collapsed[[r]]) ) {
        
        subgroupsizes_collapsed[[r]][s] 
        
        if(subgroupsizes_collapsed[[r]][s] == 1) {
          n_indic = n_indic + 1
        } else {
          n_indic = n_indic + max(subgroupsizes_collapsed[[r]][s]-1,0)
        }
        
      } 
    }
    
    
  }
  n_indic
  n_total
  
  #########################################################
  
  k_2 = 0 
  
  design_list_collapsed = list()
  
  for (r in 1:length(final[[1]])) {
    
    temp_design_list = list()
    
    for (s in 1:length(final[[1]][[r]])) {
      
      k_2 = k_2 + 1
      
      temp_design_list[[length(temp_design_list)+1]] = k_2
      
    }
    
    
    design_list_collapsed[[r]] = temp_design_list 
    
    
  }
  
  
  #########################################################
  
  
  design_matrix_new = data.frame()
  for (i in 1:sum(get(paste0("subgroupsizes",tree_levels))) ) {
    
    r = tree_levels
    indic = 0
    
    temp_class = i
    
    
    while (r>1 & indic==0) { 
      
      temp_list = final[[2]][[r]]
      
      ##
      if(any(sapply(temp_list, function(vec) temp_class %in% vec))) {
        
        temp_class = get(paste0("parents_",r))[temp_class]
        r = r-1
        
      } else {
        indic = indic + 1
      }
      ##
      
    }
    
    ###
    
    s = which(sapply(final[[1]][[r]], function(vec) temp_class %in% vec))
    
    design_matrix_new = rbind(design_matrix_new,design_list_collapsed[[r]][[s]])
    
  }
  
  
  return(design_matrix_new)
}