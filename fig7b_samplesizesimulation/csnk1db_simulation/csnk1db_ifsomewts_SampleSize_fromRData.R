# 26/03/2020
# from csnk1db data
# plot: % of wt in KO pool vs statistical power to detect a difference

# >> version if RData is ready; needs to be loaded to polish plot

library(ggpubr)
library(esc)
library(pwr)
library(dplyr)
library(tidyr)

# v0: one simulation
# v1: do multiple simulations so can have an error bar

# import ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load(file='csnk1db_ifsomewts_SampleSize_100fish10sim.RData') # takes a couple of minutes

# polish plot -------------------------------------------------------------
# highest bar is 5% KO; close to 1,000,000!

colfunc <- colorRampPalette(c('#697a87', '#f1876b'))
mycols <- colfunc(nfish+1)

# first few are > 100; will plot them as 'infinite'
# will artificially change mean to 999 / sd to NA

Simstats[which(Simstats$mean>100), 2] <- 999
Simstats[which(Simstats$mean>100), 3] <- NA

# error bar should not be < 0, as minimum sample size is 2
Simstats[which(Simstats$mean - Simstats$sd < 0), 3] <- Simstats[which(Simstats$mean - Simstats$sd < 0), 2] - 2
# Ymax value of this error bar is above 100 (i.e. is not visible on the plot); so OK to artificially reduce so does not go below 0 on plot

Simplot <- ggplot(Simstats, aes(x=nko, y=mean, fill=as.factor(nko))) +
  # geom_hline(yintercept=48, linetype=1, colour='black', size=0.5) +
  # geom_vline(xintercept=59, linetype=1, colour='black', size=0.5) + # minimum % KO if 96-well plate
  geom_bar(stat='identity', width=0.9) +
  geom_hline(yintercept=16, linetype=2, colour='#ebebeb', size=0.5) +
  geom_linerange(aes(ymin=mean-sd, ymax=mean+sd), size=0.25) +
  scale_fill_manual(values=mycols) +
  theme_minimal() +
  theme(
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    legend.position='none') +
  xlab('') + ylab('') +
  coord_cartesian(ylim=c(0, 100)) +
  scale_x_continuous(breaks=c(0, 25, 50, 75, 100), labels=c('0', '25', '50', '75', '100'))
Simplot

# export ------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='f0_csnk1db_samplesize.pdf', plot=Simplot, width=100, height=70, units='mm', useDingbats=FALSE)