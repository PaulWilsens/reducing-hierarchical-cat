#### Working directory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Functions

source("Functions/NORMAL_EXTRAVAR_AICBIC_functions.R")

#### Load previous result

result = readRDS("NORMAL_EXTRAVAR_result")

#### Hierarchical structure

subgroupsizes1 = 5 #5 classes at level 1
subgroupsizes2 = c(4,4,4,4,4) #20 classes at level 2
subgroupsizes3 = c(4,4,4,4,7,7,4,4,4,4,4,4,4,4,4,4,4,4,4,4) #86 classes at level 3

tree_levels = 3
n_leaf = sum(get(paste0("subgroupsizes",tree_levels)))

#### Computing AIC and BIC

BIC_frame_NORMAL_EXTRAVAR = data.frame()
AIC_frame_NORMAL_EXTRAVAR = data.frame()

structure_completecollapse = readRDS("true_noeffect")

for (i in 1:length(result)) {
  
  data = read.csv(paste0("Data/Data_3fact_EXTRAVAR_SIMULATED_",i,".csv"),header = T)
  data$leaf = as.factor(data$leaf)
  
  for (j in 1:5) {
    
    temp_final = result[[i]][[2]][[j]]
    
    if(identical(temp_final,structure_completecollapse)) {
      
      formula_old = as.formula(paste("data$y ~ data$leaf + data$x2_simul + data$x3_simul + data$x4_simul"))
      
      model_old = lm(formula_old)
      model_new = lm(data$y ~ data$x2_simul + data$x3_simul + data$x4_simul)
      
      ###
      
      temp = c(BIC(model_old),BIC(model_new))
      temp2 = c(AIC(model_old),AIC(model_new))
      
      
      BIC_frame_NORMAL_EXTRAVAR = rbind(BIC_frame_NORMAL_EXTRAVAR,temp)
      AIC_frame_NORMAL_EXTRAVAR = rbind(AIC_frame_NORMAL_EXTRAVAR,temp2)
      
    } else {
      
      ###
      
      new_factor = designfactor(temp_final)
      
      ###
      
      indices = data$leaf_int + 1
      new_data = as.data.frame(new_factor[indices, ])
      new_data[,1] = as.factor(new_data[,1])
      
      ###
      
      formula_old = as.formula(paste("data$y ~ data$leaf + data$x2_simul + data$x3_simul + data$x4_simul"))
      formula_new = as.formula(paste("data$y ~ new_data[,1] + data$x2_simul + data$x3_simul + data$x4_simul"))
      
      ###
      
      model_old = lm(formula_old)
      
      model_new = lm(formula_new)
      
      ###
      
      temp = c(BIC(model_old),BIC(model_new))
      temp2 = c(AIC(model_old),AIC(model_new))
      
      
      BIC_frame_NORMAL_EXTRAVAR = rbind(BIC_frame_NORMAL_EXTRAVAR,temp)
      AIC_frame_NORMAL_EXTRAVAR = rbind(AIC_frame_NORMAL_EXTRAVAR,temp2)
    }
    
    
  }
}

####


