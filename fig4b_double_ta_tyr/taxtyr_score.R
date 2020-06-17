# 200919 ta x tyr
# barplot

# libraries ---------------------------------------------------------------
library (ggplot2)
library (reshape2)
library (openxlsx)
# input -------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # wherever the script is
myfile <- "200919_taxtyr_score.xlsx"
scores <- read.xlsx (myfile)

# remove comments ---------------------------------------------------
comrows <- which(scores$plate=='Comments') : nrow(scores)
scores <- scores[-comrows,]

# remove 2dpf_counts (same as TOTAL)
scores$`2dpf_count` <- NULL

# easier to change the factor levels early on
# that is to put them in right order
goodorder <- c( which(scores$condition=='uninjected'), 
                which(scores$condition=='tyr'), 
                which(scores$condition=='ta'), 
                which(scores$condition=='both'))

scores$plate <- factor(scores$plate,
                          levels=scores$plate[goodorder])

# only biallelic ----------------------------------------------------------
# it is not possible currently to do stacked + dodged on ggplot
# would probably be too cluttered anyway
# I will do dodged barplot

eye_ko <- scores$eye1
fin_ko <- scores$tail0
scores <- cbind(scores, eye_ko, fin_ko)

# pool clutches -----------------------------------------------------------
scopool <- aggregate(cbind(mosaicism$eye_ko, mosaicism$fin_ko, mosaicism$eyeTOTAL, mosaicism$finTOTAL),
                      by=list(condition=mosaicism$condition),
                      FUN=function(x) c(sum=sum(x)))

colnames(scopool) <- c('condition', 'eye_ko', 'tail_ko', 'eye_tot', 'tail_tot')

eye_kop <- (scopool$eye_ko/scopool$eye_tot) * 100
tail_kop <- (scopool$tail_ko/scopool$tail_tot) * 100
# p is for percent

scopool <- cbind(scopool, eye_kop, tail_kop)

sco <- melt(scopool, id.vars=c('condition', 'eye_tot', 'tail_tot', 'eye_ko', 'tail_ko'))

colnames(sco)[6:7] <- c('pheno', 'perc')

sco$pheno <- factor(sco$pheno, levels=c('eye_kop', 'tail_kop'))
sco$condition <- factor(sco$condition, levels=c('uninjected', 'tyr', 'ta', 'both'))

# f0 paper ----------------------------------------------------------------

scorestack <- ggplot (data = sco, aes (x = condition, y = perc, fill = pheno)) +
  geom_bar (stat = 'identity', width=0.8, position=position_dodge(width=0.9)) + 
  # spacing between conditions: width
  # spacing between dodged bars: position_dodge(width)
  scale_fill_manual(drop = FALSE, values = c('#B3BDC4', '#98a285')) +
  scale_x_discrete(labels=c('uninjected', expression(italic('tyr')), expression(italic('ta')), 'both')) +
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
scorestack

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_tyrxta.pdf', plot=stacka4, width=65, height=57.9, unit='mm', useDingbats=FALSE) # f0 paper
