# replicate plots from Vp_Analyse_MD
# (MatLab script from Marcus)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm) # for better positions of jitter dots
library(esc)

sem <- function(x) sd(x)/sqrt(length(x))
sem_narm <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

# file paths ------------------------------------------------------------------
scatter_file <- '~/.../eucdists.csv'

exptags_file <- '~/.../exptags.csv'
grptags_file <- '~/.../grptags.csv'

# import ------------------------------------------------------------------
sct <- read.csv(scatter_file, header=FALSE)

# tidy up -----------------------------------------------------------------

# MatLab uses NaN for NA; change them

NaN2NA <- function(rc) { # rc = row or col; i.e. run on a dataframe with apply
  rc[is.nan(rc)] <- NA
  return(rc)
}

sct <- apply(sct, 2, NaN2NA)

sct <- sct[,colSums(is.na(sct))<nrow(sct)] # delete empty columns

colnames(sct) <- c('stable_wt', 'stable_het', 'stable_hom', 'f06_wt', 'f06_hom', 'f07_wt', 'f07_hom')

sct <- tbl_df(sct)

# prepare to plot --------------------------------------------------------------------

sctm <- sct %>%
  gather(key=grp, value=eucdist)

sctm <- sctm[ - which(is.na(sctm$eucdist)),] # delete NA rows; they were just paddings


# add grouping variables --------------------------------------------------

# exp: stable / f06 / f07
sctm$exp <- as.vector(sapply(as.character(sctm$grp), function(x) strsplit(x, '_')[[1]][1]))
sctm$exp <- as.factor(sctm$exp)
sctm$exp <- factor(sctm$exp, levels=c('stable', 'f07', 'f06'))

# geno: wt / het (only for stable) / hom
sctm$geno <- as.vector(sapply(as.character(sctm$grp), function(x) strsplit(x, '_')[[1]][2]))
sctm$geno <- as.factor(sctm$geno)
sctm$geno <- factor(sctm$geno, levels=c('wt', 'het', 'hom'))

# order group level (if using them)
sctm$grp <- as.factor(sctm$grp)
sctm$grp <- factor(sctm$grp, levels=c('stable_wt', 'stable_het', 'stable_hom', 
                                      'f07_wt', 'f07_hom',
                                      'f06_wt', 'f06_hom'))

# plot --------------------------------------------------------------------

mycols <- c('#697a87', '#b3bdc4', '#f1876b', '#697a87', '#74210b', '#697a87', '#c23712')

eucdistscatter <- ggplot(data=sctm, aes(x=exp, y=eucdist, colour=grp)) +
  # geom_jitter(size=2.0, position=position_jitter(0.2)) +
  geom_quasirandom(dodge.width=0.7, varwidth=TRUE, size=0.2) + # places dots better, more like they were inside a violin plot
  stat_summary(fun.y=mean, aes(group=grp), position=position_dodge(width=0.7), geom='point', colour='#020304', shape=3, size=1, stroke=1) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x=element_text(size=7, margin = margin(t = 12, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.title.x=element_blank(),
    legend.position='none') +
  ylab('Euclidean distance from wild-type') + # could remove in theme() but it is to keep plot the same width with traces
  coord_cartesian(ylim=c(0, 13)) +
  scale_x_discrete(labels=c('stable', 'F0 exp1', 'F0 exp2')) +
  scale_y_continuous(breaks=c(0, 5, 10))

eucdistscatter


# statistics --------------------------------------------------------------

# density plot
mycols <- c('#b3bdc4', '#83939f', '#5a6974', '#b3bdc4', '#F39D80', '#b3bdc4', '#EE7163')
eucdistdensity <- ggplot(sctm, aes(x=eucdist, fill=grp, colour=grp)) +
  geom_density(alpha=0.2) +
  # geom_vline(data=stats, aes(xintercept=mean, color=grp),
  #            linetype='dashed') +
  scale_fill_manual(values=mycols) +
  scale_colour_manual(values=mycols) +
  theme_minimal()
eucdistdensity

# differences between groups
summary(aov(eucdist ~ grp, data=sctm)) # p ~ 0
pairwise.t.test(sctm$eucdist, sctm$grp)
# Het vs Wt NS
# Hom vs Wt Si
# rest looks all as expected

# is F0 phenotype stronger?
# look at effect size vs wild-type
stats <- sctm %>%
  group_by(grp) %>%
  summarise_at(vars(eucdist),
               list(
                 mean=mean,
                 sd=sd,
                 sem=sem,
                 median=median,
                 min=min,
                 max=max,
                 spln=length
               ))

# need effect sizes:
  # stable_hom vs stable_wt
  # f06_hom vs f06_wt
  # f07_hom vs f07_wt

esCalc <- function(grp1, grp2, stats) { # grp1 = test / grp2 = control / stats = where to find the summary stats
  ef <- esc_mean_sd(grp1m = stats$mean[which(stats$grp==grp1)],
              grp1sd = stats$sd[which(stats$grp==grp1)],
              grp1n = stats$spln[which(stats$grp==grp1)],
              grp2m = stats$mean[which(stats$grp==grp2)],
              grp2sd = stats$sd[which(stats$grp==grp2)],
              grp2n = stats$spln[which(stats$grp==grp2)],
              es.type = "d")
  return(ef)
}

ef_stable <- esCalc('stable_hom', 'stable_wt', stats)
ef_stableHet <- esCalc('stable_het', 'stable_wt', stats)

ef_06 <- esCalc('f06_hom', 'f06_wt', stats)

ef_07 <- esCalc('f07_hom', 'f07_wt', stats)


# efStatTest --------------------------------------------------------------

# approach from here: https://stats.stackexchange.com/questions/77269/statistical-comparison-of-2-independent-cohens-ds

efStatTest <- function(grp11, grp12, grp21, grp22, stats) {
  
  # v for ef1
  ef1 <- esCalc(grp11, grp12, stats)
  spln11 <- stats$spln[which(stats$grp==grp11)]
  spln12 <- stats$spln[which(stats$grp==grp12)]
  v1 <- (1/spln11) + (1/spln12) + ef1$es^2 / sqrt((2*(spln11 + spln12)))
  
  # v of ef2
  ef2 <- esCalc(grp21, grp22, stats)
  spln21 <- stats$spln[which(stats$grp==grp21)]
  spln22 <- stats$spln[which(stats$grp==grp22)]
  v2 <- (1/spln21) + (1/spln22) + ef2$es^2 / sqrt((2*(spln21 + spln22)))
  
  # z
  z <- (ef1$es - ef2$es) / sqrt(v1 + v2)
  
  # pval
  p = 2 * pnorm(-abs(z))
  
  # reject?
  if (abs(z) > 1.96) {
    cat( '\t\t\t', grp11, ' vs ', grp12, '  //VS//  ', grp21, ' vs ', grp22, '\n',
        '\t\t\tz = ', abs(z), ' // pval = ', p,' >> effect sizes are statistically different\n', sep='')
  } else {
    cat( '\t\t\t', grp11, ' vs ', grp12, '  //VS//  ', grp21, ' vs ', grp22, '\n',
        '\t\t\tz = ', abs(z), ' // pval = ', p, ' >> effect sizes are NOT statistically different\n', sep='')
  }
}


# test statically effect sizes --------------------------------------------

# question: Is the phenotype stronger in F0?
# i.e. is effect size F0 vs SCR > effect size HOM vs WT
efStatTest('stable_hom', 'stable_wt', 'f06_hom', 'f06_wt', stats)
efStatTest('stable_hom', 'stable_wt', 'f07_hom', 'f07_wt', stats)
# >> it is not; effect size not statistically different

efStatTest('stable_het', 'stable_wt', 'f06_hom', 'f06_wt', stats)
efStatTest('stable_het', 'stable_wt', 'f07_hom', 'f07_wt', stats)
efStatTest('stable_het', 'stable_wt', 'stable_hom', 'stable_wt', stats)
# >> but significant Het vs Wt VS f0 vs scr


# test statistically standard deviations ----------------------------------
# question: Is the phenotype more variable in F0?

# sample sizes are similar which is good; should not affect too much std

# within stable
var.test(sct$stable_wt, sct$stable_het)
var.test(sct$stable_wt, sct$stable_hom) # std is larger in Hom vs Wt

# within f06 & f07
var.test(sct$f06_wt, sct$f06_hom) # f06 Hom more variant than f06 WT
var.test(sct$f07_wt, sct$f07_hom) # f07 Hom more variant than f07 WT

# as negative control: stable WT vs SCR
var.test(sct$stable_wt, sct$f06_wt)
var.test(sct$stable_wt, sct$f07_wt)

# F06/07 Hom vs Stable Hom
var.test(sct$stable_hom, sct$f06_hom) # f06 NOT more variant than stable KO
var.test(sct$stable_hom, sct$f07_hom) # f07 MORE VARIANT than stable KO

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='f0_eucdist.pdf', plot=eucdistscatter, width=75, height=85, units='mm', useDingbats=FALSE) # f0 paper