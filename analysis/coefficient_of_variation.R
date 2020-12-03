rm(list = ls())
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  print('Error! We expect two command line arguments: the directory to find the data and the directory to output the plot')
  quit()
}
in_dir = args[1]
out_dir = args[2]

if(substr(in_dir, nchar(in_dir), nchar(in_dir)) != '/'){
  in_dir = paste0(in_dir, '/')
}
if(substr(out_dir, nchar(out_dir), nchar(out_dir)) != '/'){
  out_dir = paste0(out_dir, '/')
}

print(paste('Loading pop_sizes.csv from', in_dir))
data = read.csv(paste0(in_dir, 'pop_sizes.csv'))
data = data[data$name == 'predator',]

arr = 0:500
for(time_unit in 0:500){
  arr[time_unit + 1] = data[data$time == min(data[data$time > time_unit,]$time),]$count
}
mu = mean(arr)
sigma = sd(arr)
cv = sigma / mu
print(paste('Mean:', mu))
print(paste('Standard deviation:', sigma))
print(paste('Coefficient of variance:', cv))

cv_data = data.frame(cv = cv)
print(paste('Saving cv.csv in', out_dir))
write.csv(cv_data, paste0(out_dir, 'cv.csv'))
