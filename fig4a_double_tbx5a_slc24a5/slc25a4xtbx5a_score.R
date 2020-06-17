# 261119 slc24a5 x tbx5a
# barplot

# libraries ---------------------------------------------------------------
library (ggplot2)
library (reshape2)
library (openxlsx)
# input -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # wherever the script is
myfile <- "261119_slc24a5xtbx5a_score.xlsx"
fish <- read.xlsx (myfile)


# fix factor levels -------------------------------------------------------
# easier to change the factor levels early on
# that is to put them in right order
goodorder <- c( which(mosaicism$condition=='uninjected'), 
                which(mosaicism$condition=='slc24a5'), 
                which(mosaicism$condition=='tbx5a'), 
                which(mosaicism$condition=='both'))

mosaicism$plate <- factor(mosaicism$plate,
                          levels=mosaicism$plate[goodorder])

# only biallelic ----------------------------------------------------------
# it is not possible currently to do stacked + dodged
# would probably be too cluttered anyway
# I will do dodged barplot

eye_ko <- mosaicism$eye1
fin_ko <- mosaicism$fin0
mosaicism <- cbind(mosaicism, eye_ko, fin_ko)

# pool clutches -----------------------------------------------------------
mosapool <- aggregate(cbind(mosaicism$eye_ko, mosaicism$fin_ko, mosaicism$eyeTOTAL, mosaicism$finTOTAL),
                      by=list(condition=mosaicism$condition),
                      FUN=function(x) c(sum=sum(x)))

colnames(mosapool) <- c('condition', 'eye_ko', 'fin_ko', 'eye_tot', 'fin_tot')

eye_kop <- (mosapool$eye_ko/mosapool$eye_tot) * 100
fin_kop <- (mosapool$fin_ko/mosapool$fin_tot) * 100
# p is for percent

mosapool <- cbind(mosapool, eye_kop, fin_kop)

mosa <- melt(mosapool, id.vars=c('condition', 'eye_tot', 'fin_tot', 'eye_ko', 'fin_ko'))

colnames(mosa)[6:7] <- c('pheno', 'perc')

mosa$pheno <- factor(mosa$pheno, levels=c('eye_kop', 'fin_kop'))
mosa$condition <- factor(mosa$condition, levels=c('uninjected', 'slc24a5', 'tbx5a', 'both'))

# f0 paper ---------------------------------------------------------------------

stacka4 <- ggplot (data = mosa, aes (x = condition, y = perc, fill = pheno)) +
  geom_bar (stat = 'identity', width=0.8, position=position_dodge(width=0.9)) + 
  # spacing between conditions: width
  # spacing between dodged bars: position_dodge(width)
  scale_fill_manual(drop = FALSE, values = c('#B3BDC4', '#98a285')) +
  scale_x_discrete(labels=c('uninjected', expression(italic('slc24a5')), expression(italic('tbx5a')), 'both')) +
  # theme; including legend
  theme_minimal() +
  theme(
    # grid
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    # axis text/title size/positions
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 12, b = 0, l = 0)),
    axis.text.x = element_text(size=7, margin = margin(t = -4, r = 0, b = 0, l = 0)),
    axis.text.y = element_blank(),
    legend.position='none') +
  # axis
  coord_cartesian(ylim=c(0,100))+
ylab ('%')
stacka4

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_slc24a5xtbx5a.pdf', plot=stacka4, width=65, height=57.9, unit='mm') # f0 paper