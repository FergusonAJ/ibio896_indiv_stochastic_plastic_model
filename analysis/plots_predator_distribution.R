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
print(paste('Loading final_pred_genotypes.csv from', in_dir))

data = read.csv(paste0(in_dir, 'final_pred_genotypes.csv'))

x_min = min(min(data$x), -1)
x_max = max(max(data$x), 1)
y_min = min(min(data$y), 0)
y_max = max(max(data$y), 1)
text_size = 14
print(paste('Saving final_predaotor_distribution.png in', out_dir))
ggplot(data, aes(x = x, y = y)) + 
  #geom_point(alpha = 0.1) + 
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0,0)) + 
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0,0)) + 
  stat_density2d(aes(fill=..density..), geom = "tile", contour = FALSE) +
  scale_fill_gradient(low = "white", high = "black") + 
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.5) + 
  geom_vline(xintercept = 1, linetype = 'dotted', alpha = 0.5) + 
  geom_vline(xintercept = -1, linetype = 'dotted', alpha = 0.5) + 
  theme_bw() + 
  xlab('Base value') + 
  ylab('Plasticity') + 
  labs(fill = 'Density') +
  theme(
    axis.text = element_text(size = text_size), 
    axis.title = element_text(size = text_size),
    legend.text = element_text(size = text_size), 
    legend.title = element_text(size = text_size)
  ) +
  ggsave(paste0(out_dir, 'final_predator_distribution.png'), 
         width = 8, height = 6, units = 'in')
print('Done!')
