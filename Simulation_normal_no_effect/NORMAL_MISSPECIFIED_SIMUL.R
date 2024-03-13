#### Working directory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Functions

source("Functions/NORMAL_MISSPECIFIED_SIMUL_functions.R")

#### Load packages

library(keras)
library(ggplot2)
library(cluster)
library(factoextra)
library(tidyr)

#### Hierarchical structure

n_leaf = 86
tree_levels = 3

subgroupsizes1 = 5
subgroupsizes2 = c(4,4,4,4,4)
subgroupsizes3 = c(4,4,4,4,7,7,4,4,4,4,4,4,4,4,4,4,4,4,4,4)

#### Reducing hierarchy

results_summary = list()
NN_seeds = 20 + 50*(1:5)

for (l in 1:100) {
  
  data = read.csv(paste0("Data/Data_3fact_MISSPECIFIED_SIMULATED_",l,".csv"),header = T)
  data$leaf     = as.factor(data$leaf)
  data$leaf_int = as.numeric(data$leaf_int)
  data$y = scale(data$y)
  
  results_temp  = list()
  results_temp2 = list()
  
  for (k in 1:5) {
    
    tensorflow::set_random_seed(NN_seeds[k],disable_gpu = TRUE)
    
    ### Creating NN w embedding layer 
    
    embed_dim = 2
    leaf_dim = length(unique(data$leaf_int))     
    
    
    Leaf = layer_input(shape = c(1),   dtype = 'int32', name = 'Leaf')       
    
    
    model = keras_model_sequential()%>%
      layer_embedding(input_dim = leaf_dim, output_dim = embed_dim, input_length = 1, name = 'LeafEmb') %>% 
      layer_flatten(name='Leaf_flat')  %>%    
      layer_dense(units=2) %>%             
      layer_dense(units=1)
    
    model %>% 
      compile(loss = "mse", optimizer = "adam")
    
    fit = model %>% 
      fit(x = as.matrix(data$leaf_int), y= as.matrix(data$y), epochs = 50, batch_size = 25)
    
    
    ### Extract embedding
    
    layer = get_layer(model, "LeafEmb") 
    embeddings = data.frame(layer$get_weights()[[1]])
    
    ### Aggregating embeddings
    
    embedding_list = tree_aggregation(embeddings,tree_levels,embed_dim)
    
    ### Apply top-down clustering
    
    fake_lev1   = list()                   
    G_lev1      = list()
    G_lev1[[1]] = embedding_list[[1]]
    
    listH = list()
    listV = list()
    listV[[1]] = fake_lev1
    
    for (j in 1:tree_levels) {
      
      assign(paste0("clusters_lev",j),horizontal_clustering(j,cutoff=0.7))
      listH[[length(listH)+1]] = get(paste0("clusters_lev",j))
      
      if (j < length(embedding_list)) {
        assign(paste0("fake_lev",j+1),vertical_clustering(j,cutoff=0.7))
        listV[[length(listV)+1]] = get(paste0("fake_lev",j+1))
      }
      
    }
    
    ### Compare to true structure
    
    final = condense_result(listH,listV)
    
    benchmark = readRDS("true_noeffect")
    
    result = identical(final,benchmark)
    
    ###
    
    results_temp[[length(results_temp)+1]] = result
    results_temp2[[length(results_temp2)+1]] = final
    
    print("beep")
  }
  
  temp = list(results_temp,results_temp2)
  
  results_summary[[length(results_summary)+1]] = temp
  
}

####

