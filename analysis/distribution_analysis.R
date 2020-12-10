rm(list = ls())
library(ggplot2)
library(dplyr)

# Look at a single seed of data
# Returns a trimmed data frame. One row corresponds to a bin of genotypes that WERE present
# Returns a string as an error code on failure
analyze_data = function(in_dir, out_dir){
  # Load the data, return with error code if file / data not found
  filename = paste0(in_dir, 'final_pred_genotypes.csv')
  if(!file.exists(filename)){
    return('file_error')
  }
  data = read.csv(filename)
  if(nrow(data) == 0){
    return('extinction')
  }
  
  # Split y axis into two bins [0, 0.5) and (0.5, 1]
  data$y_bin = 0
  data[data$y > 0.5,]$y_bin = rep(1, nrow(data[data$y > 0.5,]))
  # Split x axis into three bins
  # If y is in bin 0, there are three bins: (-inf, -0.5], (-0.5, 0.5), [0.5, inf]
  data$x_bin = -1
  data[data$y_bin == 0 & data$x > -0.5,]$x_bin = rep(0, nrow(data[data$y_bin == 0 & data$x > -0.5,]))
  data[data$y_bin == 0 & data$x >= 0.5,]$x_bin = rep(1, nrow(data[data$y_bin == 0 & data$x >= 0.5,]))
  # If y is in bin 1, there are three bins: (-inf, y], (-0.5, 0.5), [0.5, inf]
  data[data$y_bin == 1,]$x_bin = rep(0, nrow(data[data$y_bin ==1,]))
  data[data$y_bin == 1 &  data$x <= -0.5 & (-1 * data$x) > data$y,]$x_bin = rep(-1, nrow(data[data$y_bin == 1 & data$x <= -0.5 & (-1 * data$x) > data$y,]))
  data[data$y_bin == 1 &  data$x >= -0.5 & data$x >= data$y,]$x_bin = rep(1, nrow(data[data$y_bin == 1 &  data$x >= -0.5 & data$x >= data$y,]))
  # Combine x and y into a single bin
  data$bin = paste0('(', data$x_bin, ', ', data$y_bin, ')')
  # Some bins are actually just pieces of larger groups. Combine them.
  data$bin_name = 'Specialist: Prey 1'
  data[data$bin == '(-1, 0)',]$bin_name = rep('Specialist: Prey 1', sum(data$bin=='(-1, 0)'))
  data[data$bin == '(-1, 1)',]$bin_name = rep('Specialist: Prey 1', sum(data$bin=='(-1, 1)'))
  data[data$bin == '(1, 0)',]$bin_name =  rep('Specialist: Prey 2', sum(data$bin=='(1, 0)'))
  data[data$bin == '(1, 1)',]$bin_name =  rep('Specialist: Prey 2', sum(data$bin=='(1, 1)'))
  data[data$bin == '(0, 0)',]$bin_name =  rep('Default', sum(data$bin=='(0, 0)'))
  data[data$bin == '(0, 1)',]$bin_name =  rep('Plastic', sum(data$bin=='(0, 1)'))
  data$bin_factor = as.factor(data$bin_name)
  
  # Count genotypes in each bin
  data_grouped = dplyr::group_by(data, bin_name)
  data_summary = dplyr::summarize(data_grouped, count = n())
  # Filter out noise (< 10% of points in bin) 
  total_count = sum(data_summary$count)
  data_summary$pct = data_summary$count / total_count
  data_summary$pass = F
  data_summary[data_summary$pct >= 0.1,]$pass = T
  return(data_summary)
}

# Create a dataframe that will hold all the results (one row per seed)
data_agg = data.frame(data = matrix(nrow = 0, ncol = 11))
colnames(data_agg) = c('experiment', 'k', 'm', 's', 'seed', 'default', 'prey_1', 'prey_2', 'plastic', 'extinction', 'error')

# Load meta data (info on each condition tested)
data_meta = read.csv('~/documents/school/cur_fa20/pop_comm/ibio896_indiv_stochastic_plastic_model/analysis/analysis_data.csv')
for(cond_idx in 1:nrow(data_meta)){
  print(paste0('cond_idx: ', cond_idx))
  print(data_meta[cond_idx,])
  data_cond = data_meta[cond_idx,]
  exp = data_cond$exp
  k = data_cond$k
  m = data_cond$m
  s = data_cond$s
  seed_min = data_cond$seed_min
  seed_max = data_cond$seed_max
  for(seed in (seed_min):(seed_max)){
    in_dir = paste0('~/documents/school/cur_fa20/pop_comm/ibio896_indiv_stochastic_plastic_model/data/',seed_min,'_stopla__K_',k,'__M_',m,'__S_',s,'/', seed, '/')
    out_dir = '~/'
    data_res = analyze_data(in_dir, out_dir)
    # Did we have an error?
    if(!is.data.frame(data_res)){
      #cat(paste0('Seed: ', seed, '\n'))
      if(data_res == 'extinction'){
        #cat(paste0('\t','Extinction', '\n'))
        data_agg[nrow(data_agg) + 1,] = c(exp, k, m, s, seed, 0, 0, 0, 0, 1, NA)
      }
      else if(data_res == 'file_error'){
        #cat(paste0('\t','File error', '\n'))
        data_agg[nrow(data_agg) + 1,] = c(exp, k, m, s, seed, 0, 0, 0, 0, 0, 'File error')
      }
      else{
        #cat(paste0('\t','Unknown error', '\n'))
        data_agg[nrow(data_agg) + 1,] = c(exp, k, m, s, seed, 0, 0, 0, 0, 0, 'Unknown error')
      }
    }
    else{ # No error, clean up results and store in aggregate data frame
      #cat(paste0('Seed: ', seed, '\n'))
      for(row_idx in nrow(data_res)){
        row = data_res[row_idx,]
        #if(row$pass){
        #  cat(paste0('\t', row$bin_name, ': ', row$count, ', ', row$count, '\n'))
        #}
      }
      default = nrow(data_res[data_res$pass == T & data_res$bin_name == 'Default',]) > 0
      prey_1 = nrow(data_res[data_res$pass == T & data_res$bin_name == 'Specialist: Prey 1',]) > 0
      prey_2 = nrow(data_res[data_res$pass == T & data_res$bin_name == 'Specialist: Prey 2',]) > 0
      plastic = nrow(data_res[data_res$pass == T & data_res$bin_name == 'Plastic',]) > 0
      #print(paste0('seed: ',seed))
      #print(c(exp, k, m, s, seed, default, prey_1, prey_2, plastic, F, NA))
      data_agg[nrow(data_agg) + 1,] = c(exp, k, m, s, seed, default, prey_1, prey_2, plastic, F, NA)
    }
  }
}
# Clean up and group data by counting the number of seeds with each result profile
data_agg$had_error = !is.na(data_agg$error)
data_agg_grouped = dplyr::group_by(data_agg, experiment, k, m, s, default, prey_1, prey_2, plastic, extinction, had_error)
data_agg_summary = dplyr::summarize(data_agg_grouped, count = n())

# Assign each profile to a category (to match the paper)
data_agg_summary$category = 'Error'
data_agg_summary[data_agg_summary$extinction != 0,]$category = 'Extinction'
data_agg_summary[data_agg_summary$default != 0,]$category = 'Default'
data_agg_summary[data_agg_summary$prey_1 != 0,]$category = 'Adapt to prey 1'
data_agg_summary[data_agg_summary$prey_2 != 0,]$category = 'Adapt to prey 2'
data_agg_summary[data_agg_summary$plastic != 0,]$category = 'Plastic'
data_agg_summary[data_agg_summary$plastic != 0,]$category = 'Plastic'
data_agg_summary[data_agg_summary$plastic != 0 & data_agg_summary$prey_1 != 0,]$category = 'Two branches'
data_agg_summary[data_agg_summary$plastic != 0 & data_agg_summary$prey_2 != 0,]$category = 'Two branches'
data_agg_summary[data_agg_summary$prey_1 != 0 & data_agg_summary$prey_2 != 0,]$category = 'Two branches'
data_agg_summary[data_agg_summary$plastic != 0 & data_agg_summary$prey_1 != 0 & data_agg_summary$prey_2 != 0,]$category = 'Three branches'

# Save out raw + cleaned data
write.csv(data_agg, './data_agg.csv')
write.csv(data_agg_summary, './data_agg_summary.csv')

# Plot!
ggplot(data_agg_summary, aes(x = as.factor(m), y = count, fill = category)) + 
  geom_bar(position = 'stack', stat = 'identity')
print(data_agg_summary)

print('Done!')
