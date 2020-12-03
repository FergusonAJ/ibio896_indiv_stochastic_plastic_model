rm(list = ls())
library(ggplot2)

setwd('~/documents/school/cur_fa20/pop_comm/ibio896_indiv_stochastic_plastic_model/cpp/')
data = read.csv('pop_sizes.csv')
data = data[data$name == 'predator',]

arr = 0:500
for(time_unit in 0:500){
  #print(paste(time_unit, data[data$time == min(data[data$time > time_unit,]$time),]$count))
  arr[time_unit + 1] = data[data$time == min(data[data$time > time_unit,]$time),]$count
}
mu = mean(arr)
sigma = sd(arr)
cv = sigma / mu
print(paste('Mean:', mu))
print(paste('Standard deviation:', sigma))
print(paste('Coefficient of variance:', cv))
