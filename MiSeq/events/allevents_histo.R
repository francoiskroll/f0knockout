# septMiSeq amplican: all events/histogram

library (ggplot2)
library (reshape2)
library (schoolmath)

all.neg <- function(x) -1*abs(x)
# from https://stackoverflow.com/questions/51178646/assign-a-negative-sign-to-every-value-in-a-vector-in-r;
# add negative sign to every positive value in a vector

# import ------------------------------------------------------------------

myfile <- '~/Dropbox/phd/septMiSeq/septMiSeq_amplican/events_histo/events_all.csv'
allevents0 <- read.csv (myfile)


# add a locus column ------------------------------------------------------
# will help with filtering below

locus <- sub('_([^_]*)$', '', allevents0$seqnames)
gene <- sub('_([^_]*)$', '', locus) # for off-targets, gene will be 'off', helps with filtering below
fish <- as.numeric(sapply(as.character(allevents0$seqnames), function(x) strsplit(x, '_')[[1]][3])) # using a regex and gsub would be cleaner but does the job
allevents0 <- cbind(gene, locus, fish, allevents0)

# clean up ----------------------------------------------------------------

rows2remove <- c()

# remove any controls (should be very few as they have less mutations)
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

# by the way: same as when starting from unfolded (quantifyDiversity.R), which is a good sign

# in vs del ratio --------------------------------------------------------------

length(which(uni$type=='insertion'))
length(which(uni$type=='deletion'))
length(which(uni$type=='insertion')) / nrow(uni)
length(which(uni$type=='deletion')) / nrow(uni)

# if del > width should be negative ---------------------------------------

uni$width[which(uni$type=='deletion')] <- all.neg(uni$width[which(uni$type=='deletion')])

# collapse into counts per width ----------------------------------------------

countsperwidth <-aggregate(uni$width, by=list(uni$width), FUN=length)
colnames(countsperwidth) <- c('width', 'count')

# looks correct, tests below
length(which(uni$width==1))
length(which(uni$width==2))

# add back type, based on sign
uni[is.negative(countsperwidth$width),]


# convert counts to probability ---------------------------------------------
# probability as = proportion of the dataset
countsperwidth$prob <- countsperwidth$count / sum(countsperwidth$count)

# add deletion or insertion column ----------------------------------------
# will use it for colours
countsperwidth$deletion <- is.negative(countsperwidth$width)

# width vs freq histo -----------------------------------------------------

mycols <- c('#83939f', '#5a6974')

freqvwidth <- ggplot (data=countsperwidth, aes (x=width, y=prob, fill=deletion)) +
  geom_bar(stat='identity', width=0.7) +
  scale_fill_manual(values=mycols) +
  theme_minimal() +
  theme(
    # grid
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    # axis text
    axis.title.x=element_text(size=9, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x = element_text(size=7),
    axis.text.y = element_text(size=7),
    legend.position='none'
  ) +
  coord_cartesian(xlim=c(-50, 50)) +
  xlab('indel length (bp)') + ylab('frequency')
freqvwidth

# export ------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_freqvwidth.pdf', plot=freqvwidth, width=90, height=60, unit='mm', useDingbats=FALSE) # f0 paper
