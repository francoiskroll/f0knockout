# model of biallelic knockout by frameshift
# compare different efficiencies
# 30/06/19

library(ggplot2)
library(reshape2)

# function ---------------------------------------------------------------
guide2prob <-  function (n, ploidy, efficiency) {
  pmiss <- 1 - (0.66 * efficiency)
  biallel <- (1 - pmiss^n) ^ ploidy
  return (biallel)
}

ploidy <- 2

# 80% cutting -------------------------------------------------------------
efficiency <- 0.8

nguides <- c()
probs <- c()

for (i in 1:5) {
  nguides[i] <- i
  probs[i] <- guide2prob(i, ploidy, efficiency)
}

probs80 <- as.data.frame (cbind (nguides, probs))

# 100% cutting ------------------------------------------------------------
efficiency <- 1.0

nguides <- c()
probs <- c()

for (i in 1:5) {
  nguides[i] <- i
  probs[i] <- guide2prob(i, ploidy, efficiency)
}

probs100 <- as.data.frame (cbind (nguides, probs))

# together ----------------------------------------------------------------
probs2 <- cbind (probs80, probs100)
probs2 <- probs2[,-3]
colnames(probs2) <- c('nguides', 'prob80', 'prob100')


# melt for plotting -------------------------------------------------------
probs2m <- melt(probs2, id='nguides')
colnames(probs2m) <- c('nguides', 'grp', 'value')


# plot --------------------------------------------------------------------
mycols <- c('#B1C0C8', '#EE7163')
axislabs <- c('number of loci', 'probability of knockout')

biallelic_80vs100 <- ggplot (probs2m, aes(x=nguides, y=value, colour=grp)) +
  geom_hline(yintercept=0.9, linetype=2, colour='#ebebeb', size=0.5) +
  geom_line(size=0.5) +
  geom_point(size=0.6) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  theme( # remove the vertical grid lines
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y=element_blank(),
    legend.position='none') +
  coord_cartesian(xlim=c(1,5), ylim=c(0,1)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  ylab(axislabs[2]) + xlab(axislabs[1])
biallelic_80vs100

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # where the script is
ggsave (filename='f0_f0theory.pdf', plot=biallelic_80vs100, width=75, height=59, units='mm', useDingbats=FALSE) # f0 paper