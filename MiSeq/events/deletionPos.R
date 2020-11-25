# septMiSeq amplican: all events/histogram

# first ~ 85 rows are pasted from allevents_histo.R

library (ggplot2)
library (reshape2)
library (schoolmath)
library(dplyr)

all.neg <- function(x) -1*abs(x)
# from https://stackoverflow.com/questions/51178646/assign-a-negative-sign-to-every-value-in-a-vector-in-r;
# add negative sign to every positive value in a vector

# import ------------------------------------------------------------------

myfile <- '~/.../events_all.csv'
allevents0 <- read.csv (myfile)


# add a locus column ------------------------------------------------------
# will help with filtering below

locus <- sub('_([^_]*)$', '', allevents0$seqnames)
gene <- sub('_([^_]*)$', '', locus) # for off-targets, gene will be 'off', helps with filtering below
fish <- as.numeric(sapply(as.character(allevents0$seqnames), function(x) strsplit(x, '_')[[1]][3])) # using a regex and gsub would be cleaner but does the job
allevents0 <- cbind(gene, locus, fish, allevents0)

# clean up ----------------------------------------------------------------

rows2remove <- c()

# remove any controls (should be few as they have less mutations)
rows2remove <- c(rows2remove, which(allevents0$fish==0)) # 24684 rows
# surprisingly high
# but 94% mismatches (see below)
length( which(allevents0[which(allevents0$fish==0), 11] == 'mismatch') ) / length(which(allevents0$fish==0)) # column 11 is type

# remove all FALSE consensus
rows2remove <- c(rows2remove, which(allevents0$consensus==FALSE))

# remove all mismatch
rows2remove <- c(rows2remove, which(allevents0$type=='mismatch'))

# remove off-target samples
rows2remove <- c(rows2remove, which(allevents0$gene=='off'))

# remove slc24a5 AC samples
rows2remove <- c(rows2remove, which(allevents0$locus=='slc24a5_AC'))

# remove low coverage samples
lowcovrows <- c( which(allevents0$seqnames=='slc24a5_AB_2') , which(allevents0$seqnames=='off_AB1_4') , 
                 which(allevents0$seqnames=='tyr_AA_3') , which(allevents0$seqnames=='scn1lab_AC_2'), 
                 which(allevents0$seqnames=='scn1lab_AC_3'), which(allevents0$seqnames=='scn1lab_AC_4') )
rows2remove <- c(rows2remove, lowcovrows)

rows2remove <- unique(sort(rows2remove))
length(rows2remove)

allevents <- allevents0[-rows2remove,]

# remove strand -
# after filtering above, there remain 3% of mutations on negative strand
# they may represent mutations that were found solely on the negative strand (not certain) -- safer to delete them
length(which(allevents$strand=='-')) / nrow(allevents) * 100 # ~ 3% on negative strand
allevents <- allevents[-which(allevents$strand=='-'),]

# keep only unique events ----------------------------------------------
# closest I can to unique repair events
# i.e. keep only unique event in each sample

# duplicates defined as two events that have:
# same gene (col1)
# same locus (col2)
# same fish (col3)
# same sample (col4)
# same start position (col5)
# same end position (col6)
# same width (col7)
# same originally (col9) (NA ok)
# same replacement (col10) (NA ok)
# same type (col11)

dups <- duplicated(allevents[, c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11)])
# same as nrow(allevents) because for each line: duplicated TRUE or FALSE
duplicates <- allevents[which(dups),]
nrow(duplicates)
uni <- allevents[-which(dups),]
nrow(uni)

# keep only deletions -----------------------------------------------------
unidel <- subset(uni, type=='deletion')

# prepare position deleted vector -----------------------------------------

# idea is: pre-allocate a vector of length 401 (from -200 to 200) with all 0
# then loop through each deletion (each row of unidel)
# and add +1 to each deleted position

# pre-allocate
posdel <- rep(0, length(-200:200))

# ! element 201 corresponds to position 0
# i.e. need to do + 201 to any position to match it to the vector

# loop through each deletion (each row in unidel)
for (de in 1:nrow(unidel)){
  delb <- unidel[de,c('start')] : unidel[de,c('end')] # deleted bases
  posdel[delb+201] <- posdel[delb+201]+1
}

# put it nicely in a dataframe
posDel <- as.data.frame(matrix(nrow=length(posdel), ncol=2))
colnames(posDel) <- c('pos', 'ndel') # position / number of times deleted
posDel$pos <- -200:200
posDel$ndel <- posdel

# better Y axis: normalised in frequency deleted
posDel$freq <- posDel$ndel / nrow(unidel)

# plot --------------------------------------------------------------------

posDel[which.max(posDel$ndel),] # peak is at -4

mycols <- '#83939f'

posdelPlot <- ggplot(posDel, aes(x=pos, y=freq)) +
  geom_bar(stat='identity', width=0.8, fill=mycols) +
  coord_cartesian(xlim=c(-50, 50)) +
  theme_minimal() + 
  theme(
    # grid
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    # axis text
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=7),
    legend.position='none'
  ) +
  ylab('frequency of deletion')
posdelPlot

# note --------------------------------------------------------------------

# positions are given based on the PAM

# example slc24a5_AD, binding site is on negative strand
tmp <- subset(unidel, seqnames=='slc24a5_AD_2')
tmp[order(-tmp$counts),]
# and look in parallel at id_report.html Variants plot for slc24a5_AD_2
# 1- deleted nucleotides are included in positions
# 2- if gRNA binding site is on -: everything is inverted (which in practice means only the sign of the positions are changed)
# inverted; i.e. read Variants plot from right to left to read from 5' to 3'

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_posdel.pdf', plot=posdelPlot, width=110, height=60, unit='mm', useDingbats=FALSE) # f0 paper