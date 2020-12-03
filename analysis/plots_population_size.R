rm(list = ls())
library(ggplot2)

setwd('~/documents/school/cur_fa20/pop_comm/ibio896_indiv_stochastic_plastic_model/cpp/')
data = read.csv('pop_sizes.csv')
data$name_str = 'Predator'
data[data$name == 'prey1',]$name_str = 'Prey 1'
data[data$name == 'prey2',]$name_str = 'Prey 2'

text_size = 14
color_pred   = '#b2df8a'
color_prey_1 = '#1f78b4'
color_prey_2 = '#a6cee3'
max_time_val = max(data$time)

ggplot(data, aes(x = time, y = count)) + 
  geom_hline(yintercept = 2333, linetype = 'dashed', alpha = 0.5) +
  annotate(geom = 'text', x = (9/10) * max_time_val, y = 2333 + 50, label = "Prey carrying capacity", alpha = 0.5, size = (5/14) * text_size) +
  geom_line(aes(x = time, y = count, color = name_str)) + 
  scale_color_manual(values = c(color_pred, color_prey_1, color_prey_2)) + 
  xlab('Time units') + 
  ylab('Population size') +
  theme(
    legend.position = 'bottom',
    axis.text = element_text(size = text_size), 
    axis.title = element_text(size = text_size),
    legend.text = element_text(size = text_size),
    legend.title = element_blank()
  ) + 
  ggsave('~/documents/school/cur_fa20/pop_comm/ibio896_indiv_stochastic_plastic_model/analysis/plots/pop_sizes.png', 
         width = 20, height = 6, units = 'in')

