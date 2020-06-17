# MiSeq: offtarget
# scatterplot; each datapoint = 1 fish/1 off-target

library(openxlsx)
library(ggplot2)
library (reshape2)

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

summ <- summ0[-rows2remove,]

# only keep the slc24a5 on/off samples -----------------------------------------
# crRNA_name AA / AB / AD
# AC and AG's off-targets were not sequenced
rows2keep <- c( which(summ$crRNA_name == 'slc24a5_AA') , which(summ$crRNA_name == 'slc24a5_AB') , 
                which(summ$crRNA_name == 'slc24a5_AD') )

summoff <- summ[rows2keep,]

# fix factor levels -------------------------------------------------------

summoff$ID <- as.factor(summoff$ID) # by default levels are alphabetical

summoff$locus <- as.factor(summoff$locus)

goodorder <- levels(summoff$locus)[c(10, 11, 12, 1:9)]

summoff$locus <- factor(summoff$locus,
                        levels=goodorder)

# compute % ---------------------------------------------------------------
# total reads is Reads_Filtered
# reads cannot be in two categories at once; barplot adds it on top
# need to compute reads not frameshifted, and add that category on top of reads frameshifted
edit <- (summoff$Reads_Edited/summoff$Reads_Filtered) # % of filtered reads that have edits
frameshift <- (summoff$Reads_Frameshifted/summoff$Reads_Filtered) # ! % of filtered reads that have frameshift
frameshift[is.nan(frameshift)] <- 0 # if Reads_Edited = 0, generates NaN; replace NaN by 0
nonframeshift <- ((summoff$Reads_Edited - summoff$Reads_Frameshifted) / summoff$Reads_Filtered)
nonframeshift[is.nan(nonframeshift)] <- 0

# off-targets do not have a ON score, so can classify based on that
isOFF <- is.na(summoff$on_score)

summoff <- cbind(summoff, edit, frameshift, nonframeshift, isOFF)

# change so adds info ON/OFF/Control; will put legend based on that
summoff$isOFF[which(summoff$isOFF==TRUE)] <- 'off'
summoff$isOFF[which(summoff$isOFF==FALSE)] <- 'on'

summoff$isOFF[which(summoff$Control==1)] <- 'control'

summoff$isOFF <- as.factor(summoff$isOFF)
summoff$isOFF <- factor(summoff$isOFF,
                        levels=c('on', 'off', 'control'))

# scatterplot — colour is on/off ---------------------------------------------
# for f0 paper

mycols <- c('#EE7163', '#B1C0C8', 'white')

offdotplot <- ggplot(summoff, aes(x=locus, y=edit, fill=isOFF)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.2, dotsize=0.8, stroke=0.3) +
  theme_minimal() +
  theme(
        # panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.x = element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=7),
        legend.position='none'
        ) +
  xlab('') + ylab('% mutated reads') +
  scale_fill_manual(values=mycols, labels=c('on', 'off', 'control'))
offdotplot

# AD1 is the only off-target that is positive
# anything special about it?

# actually lowest off score within off-targets (58 vs 59-87)
# and only one with 2 mismatches in sgRNA sequence (vs 3-4)

# numbers to report -------------------------------------------------------

summoff$edit[which(summoff$isOFF=='off')] * 100 # apart from off_AD1, maximum mutation % is 0.15%
summoff$frameshift[which(summoff$isOFF=='off')] * 100

# results for off-target AD1
summoff[44:48, c(1, 58:61)]

# % frameshifted?
summoff[44:48, 59]
mean(summoff[44:48, 59]*100)

# number of mismatches & score
which(colnames(summoff)=='off_numMismatch') # col38
summoff[44:55, c(1, 37, 38)]

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_offtarget.pdf', plot=offdotplot, width=92.3, height=60, unit='mm', useDingbats=FALSE)
