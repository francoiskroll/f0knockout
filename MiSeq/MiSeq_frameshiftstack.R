# MiSeq ampliCan: frameshift stacked barplots

# Figure 2A & Figure 3B (right)
# also some summary statistics reported in text

library (openxlsx)
library (ggplot2)
library (reshape2)
library(tidyr)
library(dplyr)

# import ------------------------------------------------------------------

myfile <- '~/.../MiSeq_amplicanresults.xlsx'
summ0 <- read.xlsx (myfile)

# remove low coverage samples ---------------------------------------------
# 30x is probably a good threshold
# it sounds like the correct way to report coverage for paired-end reads is to NOT double count
# i.e. what ampliCan does, but not IGV

# 6 samples have lower than 30x coverage

rows2remove <- c()

lowcov_samples <- summ0$ID[which(summ0$Reads_Filtered < 30)]

rows2remove <- c(rows2remove, which(summ0$Reads_Filtered < 30))

# remove other samples ----------------------------------------------------

# remove off-target samples
off_rows <- which(summ0$ID=='off_AA1_0') : which(summ0$ID=='off_AD3_4') 
# 44 samples: correct (5 samples x 3 off-targets x 3 guides; minus 1 sample low coverage)
rows2remove <- c(rows2remove, off_rows)

# remove slc24a5 AC; was sequenced is for headloop PCR comparisons
slc24a5ac_rows <- which(summ0$ID=='slc24a5_AC_0') : which(summ0$ID=='slc24a5_AC_4')
rows2remove <- c(rows2remove, slc24a5ac_rows)

summ <- summ0[-rows2remove,]

# fix factor levels -------------------------------------------------------

# this is 'in order of appearance'
goodorder <- c(which(summ$gene=='slc24a5'), which(summ$gene=='tyr'), which(summ$gene=='tbx16'), which(summ$gene=='tbx5a'), which(summ$gene=='ta'),
               which(summ$gene=='slc45a2'), which(summ$gene=='mitfa'), which(summ$gene=='mpv17'), which(summ$gene=='csnk1db'), which(summ$gene=='scn1lab'))

if(length(goodorder) != nrow(summ)) stop('problem')

summ$ID <- as.factor(summ$ID) # by default levels are alphabetical

summ$ID <- factor(summ$ID,
                  levels=summ$ID[goodorder])

summ <- summ[goodorder,] # get the whole dataframe in that order

# compute % ---------------------------------------------------------------
# total reads is Reads_Filtered
# reads cannot be in two categories at once; barplot adds it on top
# need to compute reads not frameshifted, and add that category on top of reads frameshifted
edit <- (summ$Reads_Edited/summ$Reads_Filtered) * 100 # % of filtered reads that have edits
frameshift <- (summ$Reads_Frameshifted/summ$Reads_Filtered) * 100 # ! % of filtered reads that have frameshift
# not to be confused with Pframeshift = Pframeshift once edit happens
frameshift[is.nan(frameshift)] <- 0 # if Reads_Edited = 0, generates NaN; replace NaN by 0
nonframeshift <- ((summ$Reads_Edited - summ$Reads_Frameshifted) / summ$Reads_Filtered) * 100 # ! edited but not frameshifted
nonframeshift[is.nan(nonframeshift)] <- 0

# Reads_Wt <- summ$Reads_Filtered - summ$Reads_Edited

summ <- cbind(summ, edit, frameshift, nonframeshift)


# add fish num column -----------------------------------------------------

summ$fish <- as.numeric(sapply(as.character(summ$ID), function(x) strsplit(x, '_')[[1]][3])) # using a regex and gsub would be cleaner but does the job

# barplot %edit per fish per locus ----------------------------------------

stack1 <- ggplot (data=summ, aes (x=ID, y=edit)) +
  geom_bar(stat='identity')
stack1

# melt for stacked barplots -----------------------------------------------
summ2 <- cbind.data.frame(summ$ID, summ$gene, summ$crRNA_name, summ$fish, summ$Direction, summ$pam_seq, summ$Group, summ$Control, summ$Reads_Filtered, summ$edit, summ$frameshift, summ$nonframeshift)
colnames(summ2) <- c('ID', 'gene', 'locus', 'fish', 'direction', 'pam_seq', 'group', 'control', 'cov', 'edit', 'frameshift', 'nonframeshift')

summ_m <- melt(summ2, id.vars=c('ID', 'gene', 'locus', 'fish', 'direction', 'pam_seq', 'group', 'control', 'cov', 'edit'))
colnames(summ_m) <- c('ID', 'gene', 'locus', 'fish', 'direction', 'pam_seq', 'group', 'control', 'cov', 'edit', 'type', 'prop')

summ_m$type <- factor(summ_m$type, levels=c('nonframeshift', 'frameshift'))

# f0 paper ----------------------------------------------------------------

controlsamples <- as.character(summ_m$ID[summ_m$control])

# get fish number in right order for x axis
levels(summ_m$ID)

stack2 <- ggplot (data=summ_m, aes (x=ID, y=prop, fill=type)) +
  geom_bar(stat='identity') +
  scale_fill_manual(drop = FALSE, values = c('#5a6974', '#f1876b'), labels=c('edited', 'frameshifted')) +
  # theme; including legend
  theme_minimal() +
  theme(
    # grid
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    # axis text/title size/positions
    axis.title.x=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(size=5, margin = margin(t = -8, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    legend.position='none') +
  # axis
  coord_cartesian(ylim=c(0,100)) +
  ylab ('% reads') + xlab ('') +
  scale_x_discrete(labels=summ_m$fish)
stack2

# what is distribution % frameshift ---------------------------------------

# question is: what is the % of frameshift once edit has happened
# need to compute % frameshift in edit

# summ; summ2 have same samples
# summ0 is + low_cov; off; slc24a5 AC

# remove control samples
rows2remove <- c()
rows2remove <- c(rows2remove, which(summ2$control)) # 32 control samples; OK

rows2remove <- unique(sort(rows2remove))
summ_inj <- summ2[-rows2remove,] # 123 samples now; OK

# compute % edited in injected
mean(summ_inj$edit) # 86.66 %
sd(summ_inj$edit) # 17.85 %

length(which(summ_inj$edit == 100)) # 8 samples reached 100%
length(which(summ_inj$edit > 95)) # 53 samples over 95% (of 123)
length(which(summ_inj$edit > 90)) # 74 samples over 90% (of 123)


# is there any sample where edit = 0?
# (as question = what is the probability of frameshift when edit happens)
rows2remove <- c()
rows2remove <- c(rows2remove, which(summ_inj$edit==0)) # 1 sample
# there is 1 non-control samples where edit/frameshift/nonframeshift = 0
# ta_AA_3
# >> keeping it would not change anything as frameshift / edit would return Inf
summ_inj2 <- summ_inj[-rows2remove,]

# compute % frameshift within edit reads
frameshift_perc <- (summ_inj2$frameshift / summ_inj2$edit) * 100
frameshift_perc[is.nan(frameshift_perc)] <- 0
summ_inj2 <- cbind (summ_inj2, frameshift_perc)

mean(summ_inj2$frameshift_perc) # mean is 65.4
# so assumption that frameshift happens 2/3 seems accurate

# summary stats on dataset ---------------------------------------

# TOTAL NUMBER OF SAMPLES/READS

# with controls + injected (not low cov/off-targets/AC)
length(summ$ID) # 155
sum(summ$Reads_Filtered) # 150,351

# only injected samples
length(summ_inj$ID) # 123
sum(summ_inj$cov) # 122,360
# >> to be consistent with below (i.e. not taking off-targets), I will report 100,000+

# NUMBER OF LOCI
length(unique(summ0$locus)) # 42
length(unique(summ$locus)) # 32
length(unique(summ_inj$locus)) # 32
# >> I will say 30+ (depends if include off-targets or no)

length(unique(summ0$gene)) # 19
length(unique(summ$gene)) # 10
length(unique(summ_inj$gene)) # 19
# will report 10

# COVERAGE
mean(summ_inj$cov) # 994.8
sd(summ_inj$cov) # 630.7

# MINIMUM MUTATION %
# ! removed some low edit samples
# use summ, then exclude controls

summarystats <- summ[!summ$Control,] %>% # i.e. all except controls
  tbl_df(.) %>%
  group_by(locus) %>%
  summarise_at(vars(edit),
               list(
                 mean=mean,
                 sd=sd,
                 median=median,
                 min=min,
                 max=max
               ))

summarystats$locus[which(summarystats$mean==min(summarystats$mean))]
min(summarystats$mean)
# tbx16 AB is 40.6%

# MUTATION % IN SLC24A5 AC
tmp <- subset(summ0, locus=='slc24a5_AC' & !Control)
(tmp$Reads_Edited / tmp$Reads_Filtered) * 100

# barplots for headloop ---------------------------------------------------

# get the samples from summ0

hlsamples <- c('slc24a5_AA_0', 'slc24a5_AA_3',
               'slc24a5_AB_0', 'slc24a5_AB_4',
               'slc24a5_AC_0', 'slc24a5_AC_3',
               'slc24a5_AD_0', 'slc24a5_AD_3',
               'slc24a5_AG_0', 'slc24a5_AG_4')

summhl <- summ0[match(hlsamples, summ0$ID),]
summhl$ID <- as.factor(summhl$ID) # order seems right

# same maths as above
edit <- (summhl$Reads_Edited/summhl$Reads_Filtered) * 100 # % of filtered reads that have edits
frameshift <- (summhl$Reads_Frameshifted/summhl$Reads_Filtered) * 100 # ! % of filtered reads that have frameshift
frameshift[is.nan(frameshift)] <- 0
nonframeshift <- ((summhl$Reads_Edited - summhl$Reads_Frameshifted) / summhl$Reads_Filtered) * 100
nonframeshift[is.nan(nonframeshift)] <- 0
summhl <- cbind(summhl, edit, frameshift, nonframeshift)

# add fish column
summhl$fish <- as.numeric(sapply(as.character(summhl$ID), function(x) strsplit(x, '_')[[1]][3])) # using a regex and gsub would be cleaner but does the job

# melt for stack barplot
summhl2 <- cbind.data.frame(summhl$ID, summhl$gene, summhl$crRNA_name, summhl$fish, summhl$Direction, summhl$pam_seq, summhl$Group, summhl$Control, summhl$Reads_Filtered, summhl$edit, summhl$frameshift, summhl$nonframeshift)
colnames(summhl2) <- c('ID', 'gene', 'locus', 'fish', 'direction', 'pam_seq', 'group', 'control', 'cov', 'edit', 'frameshift', 'nonframeshift')

summhl_m <- melt(summhl2, id.vars=c('ID', 'gene', 'locus', 'fish', 'direction', 'pam_seq', 'group', 'control', 'cov', 'edit'))
colnames(summhl_m) <- c('ID', 'gene', 'locus', 'fish', 'direction', 'pam_seq', 'group', 'control', 'cov', 'edit', 'type', 'prop')

summhl_m$type <- factor(summhl_m$type, levels=c('nonframeshift', 'frameshift'))

# PLOT

controlsamples <- as.character(summ_m$ID[summ_m$control])

# get fish number in right order for x axis
levels(summ_m$ID)

stackhl <- ggplot (data=summhl_m, aes (x=ID, y=prop, fill=type)) +
  geom_bar(stat='identity', width=0.9) +
  scale_fill_manual(drop = FALSE, values = c('#5a6974', '#f1876b'), labels=c('edited', 'frameshifted')) +
  # theme; including legend
  theme_minimal() +
  theme(
    # grid
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor = element_blank(),
    # axis text/title size/positions
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    # axis.text.x = element_text(size=7, margin = margin(t = -4, r = 0, b = 0, l = 0), angle=45, hjust=1),
    axis.text.x=element_blank(),
    axis.text.y = element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    legend.position='none') +
  # axis
  coord_cartesian(ylim=c(0,100)) +
  # scale_x_discrete(labels=rep(c('uninjected', 'injected'), 5)) +
  ylab ('% reads')
stackhl

# export ------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set directory to wherever the script is
ggsave('frameshift_barplot.pdf', plot=stack2, width=165, height=65, unit='mm', useDingbats=FALSE) # f0 paper
ggsave('frameshift_headloop.pdf', plot=stackhl, width=45, height=54, unit='mm', useDingbats=FALSE)
