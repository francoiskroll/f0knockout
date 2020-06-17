# 111019 crystal
# unviable barplot

# libraries ---------------------------------------------------------------
library (ggplot2)
library (reshape2)
library (openxlsx)

# input -------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # wherever the script is
myfile <- "111019_crystal_fish.xlsx"
fish <- read.xlsx (myfile)

# remove comments ---------------------------------------------------
comrows <- which(fish$plate=='Comments') : nrow(fish)
fish <- fish[-comrows,]

# fix factors -------------------------------------------------------------
# easier to change the factor levels early on
# that is to put them in right order
goodorder <- c( which(fish$condition=='uninjected'), 
                which(fish$condition=='crystal'))

fish$plate <- factor(fish$plate,
                     levels=fish$plate[goodorder])

# in percentage of 1dpf --------------------------------------------------------
perc <- (fish$totalunviable / fish$`1dpf`) * 100
fish <- cbind(fish, perc)
colnames(fish)[ncol(fish)] <- 'totalunviable_perc'

# barplot - normalised in % --------------------------------------------

full <- ggplot (data = fish, aes (x = plate, y = totalunviable_perc)) +
  geom_bar (stat = 'identity') +
  theme_minimal() +
  ylab ('% not OK') + xlab ('Plate') +
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
                which(fishpool$condition=='crystal'))

fishpool$condition <- factor(fishpool$condition,
                             levels=fishpool$condition[goodorder])


# plot --------------------------------------------------------------------
# not actually used as just a single value

poolbar <- ggplot (data = fishpool, aes (x = condition, y = unviable_perc_norm)) +
  geom_bar (stat = 'identity', colour='#595E60', size=0.5) +
  theme_minimal() +
  theme(
    # grid
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    # axis text/title size/positions
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=5, margin = margin(t = 0, r = -1, b = 0, l = 0))) +
  ylab ('% toxicity') + xlab ('number of guides') +
  coord_cartesian(ylim=c(0, 100)) +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100))
poolbar