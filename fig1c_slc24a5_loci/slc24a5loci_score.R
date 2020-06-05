# slc24a5 loci scores
# stacked barplot of eye scores

# libraries ---------------------------------------------------------------
library (ggplot2)
library (reshape2)
library (openxlsx)

# input -------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # wherever the script is
myfile <- "221019_slc24a5loci_score.xlsx"
scores <- read.xlsx (myfile)

# remove comments ---------------------------------------------------
comrows <- which(scores$plate=='Comments') : nrow(scores)
scores <- scores[-comrows,]

# remove 2dpf_counts (same as TOTAL)
scores$`2dpf_count` <- NULL

# clean ---------------------------------------------------
# excluding:
  # plate 1019_1049_4x; was injected when inflated (labelled as such at the time of injections)
  # plate 1252_1321_2x; I am using only first injection of each clutch whenever possible
    # (second injection is unlikely to have been injected while cell still flat)
plates2remove <- c('1019_1049_4x', '1252_1321_2x')
scores <- scores[-match(plates2remove, scores$plate), ]

sco <- melt(scores, id.vars=c('plate', 'condition', 'clutch', 'TOTAL'))

colnames(sco)[5:6] <- c('score', 'n')

sco$score <- factor(sco$score, levels=c(1,2,3,4,5))
sco$condition <- factor(sco$condition, levels=c('uninjected', '1x', '2x', '3x', '4x'))

# stacked barplot - not normalised ----------------------------------------
ggplot() + 
  geom_bar(aes(y = n, x = plate, fill = score), data = sco,
           stat="identity")

# pool clutches -----------------------------------------------------------
scores_n <- aggregate(sco$n,
                      by=list(condition=sco$condition, score=sco$score),
                      FUN=function(x) c(sum=sum(x)))

tots_n <- aggregate(sco$TOTAL,
                    by=list(condition=sco$condition, score=sco$score),
                    FUN=function(x) c(sum=sum(x)))

scopool <- merge(scores_n, tots_n, by=c('condition', 'score'))
colnames(scopool) <- c('condition', 'score', 'n', 'tot')

# in percentage of 2dpf --------------------------------------------------------
perc <- (scopool$n / scopool$tot) * 100
scopool <- cbind(scopool, perc)

# plot ---------------------------------------------------------------------
scorestack <- ggplot (data = scopool, aes (x = condition, y = perc, fill=score)) +
  geom_bar (stat = 'identity', width=0.9) +
  scale_fill_manual(drop = FALSE, values = c('#eceff1', '#bdc5cb', '#8d9ca6', '#697a87', '#4b5860')) +
  # theme; including legend
  theme_minimal() +
  theme(
    # grid
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    # axis text/title size/positions
    axis.title.x=element_text(size=9, margin = margin(t = 4, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 8, b = 0, l = 0)), # leaving enough space to add Y axis labels
    axis.text.x = element_text(size=7, margin = margin(t = -4, r = 0, b = 0, l = 0)),
    axis.text.y = element_blank(),
    legend.position='none') +
  # axis
  ylab ('%') + xlab ('number of loci') +
  coord_cartesian(ylim=c(0,100)) +
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100)) +
  scale_x_discrete(labels=c('uninjected', '1', '2', '3', '4'))
scorestack

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_slc24a5loci_score.pdf', plot=scorestack, width=65, height=57.9, unit='mm', useDingbats=FALSE) # f0paper
