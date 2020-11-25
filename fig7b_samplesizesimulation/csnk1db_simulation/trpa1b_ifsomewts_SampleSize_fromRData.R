# 26/03/2020
# from trpa1b data
# plot: % of wt in KO pool vs statistical power to detect a difference

# >> version if RData is ready; needs to be loaded to polish plot

library(ggpubr)
library(esc)
library(pwr)
library(dplyr)
library(tidyr)

# import ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load(file='trpa1b_ifsomewts_SampleSize_100fish10sim.RData') # takes a couple of minutes

# polish plot -------------------------------------------------------------
# highest bar is 5% KO; close to 1,000,000!

colfunc <- colorRampPalette(c('#697a87', '#f1876b'))
mycols <- colfunc(nfish+1)

# first few are > 100; will plot as 'infinite'
# will artificially change mean to 999 / sd to NA
firstsim <- 14
Simstats[1:firstsim, 2] <- 999
Simstats[1:firstsim, 3] <- NA


Simplot <- ggplot(Simstats, aes(x=nko, y=mean, fill=as.factor(nko))) +
  # geom_hline(yintercept=48, linetype=1, colour='black', size=0.5) +
  # geom_vline(xintercept=28, linetype=1, colour='black', size=0.5) + # minimum % KO for n = 48 (96-well plate)
  geom_bar(stat='identity', width=0.9) +
  geom_hline(yintercept=22, linetype=2, colour='#ebebeb', size=0.5) +
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
  xlab('') + ylab('minimum sample size') +
  coord_cartesian(ylim=c(0, 100)) +
  scale_x_continuous(breaks=c(0, 25, 50, 75, 100), labels=c('0', '25', '50', '75', '100'))
Simplot

# export ------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='f0_trpa1b_samplesize.pdf', plot=Simplot, width=100, height=70, units='mm', useDingbats=FALSE)