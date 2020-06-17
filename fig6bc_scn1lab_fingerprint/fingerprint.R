# replicate plots from Vp_Analyse_MD
# (MatLab script from Marcus)
library(tidyr)
library(dplyr)
library(ggplot2)

sem <- function(x) sd(x)/sqrt(length(x))
sem_narm <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))


# functions for half correlation heatmap ----------------------------------

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# file paths ------------------------------------------------------------------
fingerprint_file <- '~/.../fingerprint.csv'

exptags_file <- '~/.../exptags.csv'
grptags_file <- '~/.../grptags.csv'

# import & clean ------------------------------------------------------------------
fgr <- read.csv(fingerprint_file, header=FALSE)

exptags <- read.csv(exptags_file, header=FALSE)
grptags <- read.csv(grptags_file, header=FALSE)

exptags[which(exptags==1),] <- 'stable'
exptags[which(exptags==2),] <- 'f06'
exptags[which(exptags==3),] <- 'f07'

grptags[which(grptags==1),] <- 'wt'
grptags[which(grptags==2),] <- 'het'
grptags[which(grptags==3),] <- 'hom'

fgr <- cbind(exptags, grptags, fgr)
fgr <- cbind( paste(fgr[,1], fgr[,2], sep='_'), fgr)

# conversions might save some hassle later
fgr[1:3] <- apply(fgr[1:3], 2, as.factor)
fgr[4:ncol(fgr)] <- apply(fgr[4:ncol(fgr)], 2, as.numeric)

# good to keep but will not use this to plot; as need parameters to be numerics
fgrClean <- fgr
colnames(fgrClean) <- c('expgrp', 'exp', 'grp', 
                        
                        'D_active_length', 'D_active_mean', 'D_active_std', 'D_active_total', 'D_active_min', 
                        'D_active_max','D_active_num', 'D_totaltimeactive', 'D_totalactivity', 'D_inactive_length',
                        
                        'N_active_length', 'N_active_mean', 'N_active_std', 'N_active_total', 'N_active_min', 
                        'N_active_max','N_active_num', 'N_totaltimeactive', 'N_totalactivity', 'N_inactive_length')

# will use this version forward
colnames(fgr) <- c('expgrp', 'exp', 'grp',
                   c(1:20))


# summarise ---------------------------------------------------------------
fgr <- tbl_df(fgr) # converts to a tibble

fgrsumm <- fgr %>% 
  gather(key=param, value=zscore, 4:ncol(fgr)) %>%
  # key will take all columns that are not factors; value the values in these columns
  # param & zscore can be any string
  group_by(expgrp, exp, grp, param) %>%
  summarise(mean=mean(zscore), sd=sd(zscore), sem=sem(zscore), n=n())

# need to change the parameter numbers to actual numbers so plots it like a timeseries
fgrsumm$param <- as.numeric(as.character(fgrsumm$param)) 

# plot --------------------------------------------------------------------

# does not plot the wt lines (all over 0 with sd +- 1)
# does not plot Ellen's Het either

fgrsumm$expgrp <- as.factor(fgrsumm$expgrp)
fgrsumm$expgrp <- factor(fgrsumm$expgrp, levels=c('stable_wt', 'stable_het', 'stable_hom', 
                                                  'f06_wt', 'f06_hom', 
                                                  'f07_wt', 'f07_hom'))

grps2plot <- c('stable_hom', 'stable_het', 'f06_hom', 'f07_hom')

# for colours;
# original orange = #f1876b; stable KO
# 8 shades darker = #c23712; F0 exp1
# 13 shades darker = #74210b; F0 exp2

mycols <- c('#b3bdc4', '#f1876b', '#c23712', '#74210b') # in order of levels, i.e. stable_het; stable_hom; f06_hom; f07_hom

fingerprint <- ggplot(filter(fgrsumm, expgrp %in% grps2plot), aes(x=param, y=mean, col=expgrp)) +
  geom_hline(yintercept=0, linetype=1, colour='#EBEBEB') +
  geom_line(size=0.2) +
  geom_pointrange(aes(ymin=mean-sem, ymax=mean+sem), size=0.2) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    
    #panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    
    axis.title.x=element_text(size=9, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    
    axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    
    legend.position='none') +
  xlab('parameter id') + ylab('') +
  coord_cartesian(xlim=c(0,20), ylim=c(-3, 12)) +
  scale_x_continuous(breaks=1:20, labels=rep(1:10, 2)) +
  scale_y_continuous(breaks=c(-3, 0, 3))
  
   #geom_vline(xintercept=10.5, linetype='dotted') # commentout; just to place the night rectangle in Illustrator

fingerprint

# correlation -------------------------------------------------------------

fgrsumm <- ungroup(fgrsumm)

fgrsumm2 <- fgrsumm[,-c(2,3,6,7,8)] # won't pivot correctly if still unique values per group/parameter

fgrsumm3 <- fgrsumm2 %>%
  group_by(expgrp) %>%
  pivot_wider(names_from=expgrp, values_from=mean)

fgrsumm3 <- fgrsumm3[order(fgrsumm3$param),]
fgrsumm3 <- fgrsumm3[, -which(names(fgrsumm3) %in% c('param', 'f06_wt', 'f07_wt', 'stable_het', 'stable_wt'))]

fgrcor <- as.data.frame(cor(fgrsumm3))

fgrcor <- get_upper_tri(fgrcor)

# plot heatmap correlation ------------------------------------------------

fgrcor <- cbind(rownames(fgrcor), fgrcor)
colnames(fgrcor)[1] <- 'trace'

fgrcor <- tbl_df(fgrcor)

fgrcorm <- fgrcor %>% 
  pivot_longer(-trace, names_to='vtrace', values_to='corr')

fgrcorm$corr <- round(fgrcorm$corr, digits=2)

corrheat <- ggplot(fgrcorm, aes(x=trace, y=vtrace, fill=corr)) + 
  geom_tile(colour=NA) + # edges around tiles
  geom_text(aes(x=trace, y=vtrace, label=corr), colour='#ffffff', size=2) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid=element_blank(),
    axis.ticks=element_blank(),
    legend.position='none',
    axis.text.x = element_text(size=7, margin = margin(t=-1, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin = margin(t=0, r=-1, b=0, l=0))
  ) +
  # scale_fill_continuous(na.value=NA) +
  scale_fill_gradient2(low='#ffffff', high='#EE7163', 
                       midpoint=0.5, limit = c(0.5,1), na.value=NA) +
  coord_fixed() +
  scale_x_discrete(labels=c('exp1', 'exp2', '-/-')) +
  scale_y_discrete(labels=c('exp1', 'exp2', '-/-'))
corrheat

# maybe delete half of it?
# see http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization


# statistics for general behaviour ----------------------------------------
# from activity plot; knockouts hypoactive during day / hyperactive during night
# is it confirmed statistically?

# probably parameter that makes most sense is totaltimeactive or totalactivity
# these should be:
# day: parameter 8 & 9
# night: parameter 18 & 19
# looks correct on plot

# extracting the datapoints I need
fgrClean <- tbl_df(fgrClean)
pairwise.t.test(fgrClean$D_totaltimeactive, fgrClean$expgrp)
pairwise.t.test(fgrClean$N_totaltimeactive, fgrClean$expgrp)

fgrClean <- tbl_df(fgrClean)
pairwise.t.test(fgrClean$D_totalactivity, fgrClean$expgrp)
pairwise.t.test(fgrClean$N_totalactivity, fgrClean$expgrp)

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

ggsave (filename='f0_fingerprint.pdf', plot=fingerprint, width=100, height=100, units='mm', useDingbats=FALSE) # f0 paper
ggsave (filename='f0_corrheat.pdf', plot=corrheat, width=31, height=31, units='mm', useDingbats=FALSE) # f0 paper
