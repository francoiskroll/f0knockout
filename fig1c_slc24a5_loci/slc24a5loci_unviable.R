# slc24a5 loci unviable
# unviable barplot

# libraries ---------------------------------------------------------------
library (ggplot2)
library (reshape2)
library (openxlsx)

# input -------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # wherever the script is
myfile <- "221019_slc24a5loci_fish.xlsx"
fish <- read.xlsx (myfile)

# remove comments ---------------------------------------------------
comrows <- which(fish$plate=='Comments') : nrow(fish)
fish <- fish[-comrows,]

# remove same plates as in _score.R ----------------------------------------
# excluding:
# plate 1019_1049_4x; was injected when inflated (labelled as such at the time of injections)
# plate 1252_1321_2x; I am using only first injection of each clutch whenever possible
# (second injection is unlikely to have been injected while cell still flat)
plates2remove <- c('1019_1049_4x', '1252_1321_2x')
fish <- fish[-match(plates2remove, fish$plate), ]

# fix factors -------------------------------------------------------------

# easier to change the factor levels early on
# that is to put them in right order
goodorder <- c( which(fish$condition=='uninjected'), 
                which(fish$condition=='1x'), 
                which(fish$condition=='2x'), 
                which(fish$condition=='3x'), 
                which(fish$condition=='4x'))

fish$plate <- factor(fish$plate,
                     levels=fish$plate[goodorder])

# barplot per plate - not normalised ----------------------------------------
ggplot() + 
  geom_bar(aes(x = plate, y = totalunviable), data = fish,
           stat="identity")

# in percentage of 1dpf --------------------------------------------------------
perc <- (fish$totalunviable / fish$`1dpf`) * 100
fish <- cbind(fish, perc)
colnames(fish)[ncol(fish)] <- 'totalunviable_perc'

# barplot per plate - percentage of 1dpf --------------------------------------------

full <- ggplot (data = fish, aes (x = plate, y = totalunviable_perc)) +
  geom_bar (stat = 'identity') +
  theme_minimal() +
  ylab ('% unviable') + xlab ('Plate') +
  coord_cartesian(ylim=c(0, 100))
full

# pool clutches -----------------------------------------------------------
unviable_n <- aggregate(fish$totalunviable,
                     by=list(condition=fish$condition),
                     FUN=function(x) c(sum=sum(x)))

tots_n <- aggregate(fish$`1dpf`,
                    by=list(condition=fish$condition),
                    FUN=function(x) c(sum=sum(x)))

fishpool <- merge(unviable_n, tots_n, by='condition')
colnames(fishpool) <- c('condition', 'unviable', 'tot')

# percentage & normalise to uninjected ------------------------------------
# alternative would be to normalise to uninjected of same clutch, but not always have it
perc <- (fishpool$unviable / fishpool$tot) * 100
fishpool <- cbind(fishpool, perc)
colnames(fishpool)[ncol(fishpool)] <- 'unviable_perc'

conperc <- fishpool$unviable_perc [which(fishpool$condition=='uninjected')]
unviable_perc_norm <- fishpool$unviable_perc - conperc
unviable_perc_norm[unviable_perc_norm<0] <- 0

fishpool <- cbind(fishpool, unviable_perc_norm)

# fix factor levels order -------------------------------------------------
goodorder <- c( which(fishpool$condition=='uninjected'), 
                which(fishpool$condition=='1x'), 
                which(fishpool$condition=='2x'), 
                which(fishpool$condition=='3x'), 
                which(fishpool$condition=='4x'))

fishpool$condition <- factor(fishpool$condition,
                             levels=fishpool$condition[goodorder])

# plot ------------------------------------------------------

# a lot of factors seem to change the width of the bars or the plot, making it impossible to align it correctly with the scores plots

# 1/ if some values are 0; width of the plot will change
fishpool[which(fishpool$unviable_perc_norm == 0),5] <- 0.001 # adds dummy null value

# 2/ if different number of digits in the Y axis (eg. 20/40/60 vs 0/100) will change width of the bars;
  # I will just remove Y axis labels, leaving enough space to add them in Illustrator

poolbar <- ggplot (data = fishpool, aes (x = condition, y = unviable_perc_norm)) +
  geom_bar (stat = 'identity', colour='#939393', fill='white', width=0.9) +
  theme_minimal() +
  theme(
    # grid
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    # axis text/title size/positions
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 8, b = 0, l = 0)), # leaving enough space to add Y axis labels
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    legend.position='none') +
  ylab ('% unviable') + xlab ('number of guides') +
  coord_cartesian(ylim=c(0, 40)) +
  scale_y_continuous(breaks=c(0, 20, 40))
poolbar

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # wherever script is
ggsave('f0_slc24a5loci_unviable.pdf', plot=poolbar, width=65, height=20, unit='mm', useDingbats=FALSE)