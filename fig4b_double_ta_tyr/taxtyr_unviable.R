# 200919 ta x tyr
# unviable barplot

# ! last timepoint for this experiment was 2dpf, as ta KO larvae die shortly after

# libraries ---------------------------------------------------------------
library (ggplot2)
library (reshape2)
library (openxlsx)
# input -------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # wherever the script is
myfile <- "200919_taxtyr_fish.xlsx"
fish <- read.xlsx (myfile)

# remove comments ---------------------------------------------------
comrows <- which(fish$plate=='Comments') : nrow(fish)
fish <- fish[-comrows,]

# easier to change the factor levels early on
# that is to put them in right order
goodorder <- c( which(fish$condition=='uninjected'), 
                which(fish$condition=='tyr'), 
                which(fish$condition=='ta'), 
                which(fish$condition=='both'))

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
colnames(fishpool) <- c('condition', 'notok', 'tot')

# normalise percentage & to uninjected ------------------------------------
# alternative would be to normalise to uninjected of same clutch, but might not always have it
perc <- (fishpool$notok / fishpool$tot) * 100
fishpool <- cbind(fishpool, perc)
colnames(fishpool)[ncol(fishpool)] <- 'notok_perc'

conperc <- fishpool$notok_perc [which(fishpool$condition=='uninjected')]
notok_perc_norm <- fishpool$notok_perc - conperc
notok_perc_norm[notok_perc_norm<0] <- 0

fishpool <- cbind(fishpool, notok_perc_norm)


# fix factor levels order -------------------------------------------------
goodorder <- c( which(fishpool$condition=='uninjected'), 
                which(fishpool$condition=='tyr'), 
                which(fishpool$condition=='ta'), 
                which(fishpool$condition=='both'))

fishpool$condition <- factor(fishpool$condition,
                             levels=fishpool$condition[goodorder])

# f0 paper ----------------------------------------------------------------

poolbar <- ggplot (data = fishpool, aes (x = condition, y = notok_perc_norm)) +
  geom_bar (stat = 'identity', colour='#939393', fill='white', width=0.9) +
  theme_minimal() +
  theme(
    # grid
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    # axis text/title size/positions
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 12, b = 0, l = 0)), # leaving enough space to add Y axis labels
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    legend.position='none') +
  ylab ('% unviable') +
  coord_cartesian(ylim=c(0, 40)) +
  scale_y_continuous(breaks=c(0, 20, 40))
poolbar

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_taxtyr_unviable.pdf', plot=poolbar, width=65, height=20, unit='mm', useDingbats=FALSE)
