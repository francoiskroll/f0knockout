# headloop score
# correlation with MiSeq results

library(openxlsx)
library(dplyr)
library(ggplot2)

# import ------------------------------------------------------------------
path <- '~/.../HLscores_all.xlsx'
hls <- read.xlsx(path) # hls = headloop scores

path <- '~/.../MiSeq_amplicanresults.xlsx'
mis <- read.xlsx(path) # mis = miseq

# is check column ok ------------------------------------------------------

if (!identical(hls$gDNA_well, hls$gDNA_wellCheck)) stop('Error')
# all good

if (identical(hls$gDNA_well, hls$gDNA_wellCheck)) { # if all good, can delete the Check column
  hls <- subset(hls, select=-gDNA_wellCheck)
}

# mean of replicates ------------------------------------------------------
hls <- tbl_df(hls)

hlsm <- hls %>% # hlsm = headloop scores mean
  group_by(sampleID) %>%
  mutate(hlmean=mean(headloop_score)) %>%
  distinct(sampleID, .keep_all=TRUE)

# get MiSeq data ----------------------------------------------------------
misd <- mis[match(hlsm$sampleID, mis$ID),c('ID', 'Reads_Filtered', 'Reads_Del', 'Reads_In', 'Reads_Edited', 'Reads_Frameshifted')] # misd = miseq data
# should be the data I need to add

# check all correct
identical(as.character(unlist(subset(hlsm, dilution=='no', sampleID))) , misd[!is.na(misd$ID),'ID']) # OK
colnames(misd) <- c('sampleIDCheck', 'reads_filtered', 'reads_del', 'reads_in', 'reads_edited', 'reads_frameshifted')

# put them together
hlsm <- tbl_df(cbind(as.data.frame(hlsm), misd))

# check sampleID ----------------------------------------------------------

if (! sum(subset(hlsm, dilution=='no', sampleID) ==
          subset(hlsm, dilution=='no', sampleIDCheck)) ==
    nrow(subset(hlsm, dilution=='no'))) 
  stop('Error')
# all good

# add locus column --------------------------------------------------------
hlsm$locus <- paste ( sapply(hlsm$sampleID, function(x) strsplit(x, '_')[[1]][1]) , 
                      sapply(hlsm$sampleID, function(x) strsplit(x, '_')[[1]][2]) ,
                      sep='_') # using a regex and gsub would be cleaner but does the job

# calculate mutated -------------------------------------------------------
hlsm$mutated <- hlsm$reads_edited / hlsm$reads_filtered

# fix dilution samples ----------------------------------------------------
hlsm[which(hlsm$dilution!='no'),'mutated'] <- 
  hlsm[match(unlist(subset(hlsm, dilution!='no', ifdilution_whichsample)) , hlsm$sampleID) , 'mutated'] *
  as.numeric(unlist(subset(hlsm, dilution!='no', dilution)))
# subset(): gives the original samples that were diluted
# match(): found their rows
# *  : multiplies by dilution factor

# remove degraded samples -------------------------------------------------
hlsm <- subset(hlsm, sampleID != 'tbx16_AB_3' & sampleID != 'tbx16_AB_4')
# clear on gel that these samples were problematic
# can also be seen from the plot that biological replicates are very close to each other in terms of headloop scores
# but was not the case for these 2 samples


# linear regression -------------------------------------------------------
# ! lm (y ~ x)
# x = predictor // y = to be predicted
# here: x = HL score // y = mutated
hlfit <- lm(mutated ~ hlmean, data=hlsm)
summary(hlfit)
# formula is y = intercept + slope * x
# here: mutated = 0.3266 + 0.6915 * HLscore

# HLscore is a significant predictor of mutated; p-value = 9.229e-5 = 0.00009229

# correlation -------------------------------------------------------------
cor(hlsm$mutated, hlsm$hlmean, method='pearson') # r=0.662

# plot --------------------------------------------------------------------

hlsm$locus <- as.factor(hlsm$locus)
hlsm$dilution <- as.factor(hlsm$dilution)

mycols <- c('#b2bbc4', '#732210', '#f0937a',
            '#60717e', '#859170', '#c23a19', '#050607')

transparent <- rgb(255,255,255,max=255,alpha=0) # gives complete transparent code

hlcorr <- ggplot(hlsm, aes(x=hlmean, y=mutated, colour=locus, shape=dilution)) +
  geom_point(size=2.0) +
  geom_segment(aes(x=-hlfit[[1]][1]/hlfit[[1]][2], xend=1.00,
                   y=0, yend=hlfit[[1]][1] + hlfit[[1]][2]), # checked; it gives the same result as geom_smooth()
               size=0.5, colour='#5a6974', linetype=2) +
  # geom_vline(xintercept=0.6) + # temporarily to add line in Illustrator
  # geom_hline(yintercept=0.733) + # mutated intercept = 73.3%
  scale_shape_manual(values=c(17, 15, 19)) + # order = 0.25 = triangle / 0.5 = square / no = standard point
  scale_colour_manual(values=mycols) +
  theme_minimal() + 
  theme(
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=7),
    axis.title.x=element_text(size=9),
    axis.title.y=element_text(size=9),
    legend.position='none'
  ) +
  ylab('% mutated reads') + xlab('headloop score') +
  coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
  scale_y_continuous(labels=c('0', '25', '50', '75', '100'))
hlcorr

# export plot -------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_hlcorr.pdf', plot=hlcorr, width=92.3, height=92.3, unit='mm', useDingbats=FALSE)

# what you would calculate as mutated reads with hl score alone -----------
# i.e. you take hl score, and guess the % mutated reads from the line of best fit
estimateMutReads <- function(hl) {
  mut <- hlfit[[1]][1] + hlfit[[1]][2] * hl
  if(mut > 1.00) {mut <- 1.00} # ! > 100% mutated reads cannot exist, so if any estimate above 1.00, bring it down to 1.00
  if(mut < 0) {mut <- 0} # ! negative mutated reads cannot exist, so if any estimate below 0.00, bring it up to 0.00
  return(as.numeric(mut))
}

hlsm$mutated_estimate <- sapply(hlsm$hlmean, estimateMutReads)
hlsm$error <- abs(hlsm$mutated - hlsm$mutated_estimate)
mean(hlsm$error)
sd(hlsm$error)
nrow(subset(hlsm, dilution=='no' & error<0.2)) / nrow(subset(hlsm, dilution=='no'))
# 18/23 = 0.78 = above 3 times out of 4

# ns ----------------------------------------------------------------------
nrow(hlsm) # total n of samples
nrow(subset(hlsm, dilution=='no')) # n original samples
nrow(subset(hlsm, dilution!='no')) # n diluted samples
