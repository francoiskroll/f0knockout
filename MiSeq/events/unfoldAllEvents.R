# MiSeq amplican: all events

library (ggplot2)
library (reshape2)
library (schoolmath)

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
fish <- as.numeric(sapply(as.character(allevents0$seqnames), function(x) strsplit(x, '_')[[1]][3])) 
# using a regex and gsub would be cleaner but does the job
allevents0 <- cbind(gene, locus, fish, allevents0)

# clean up ----------------------------------------------------------------

rows2remove <- c()

# remove all FALSE consensus
rows2remove <- c(rows2remove, which(allevents0$consensus==FALSE))

# remove all mismatch
rows2remove <- c(rows2remove, which(allevents0$type=='mismatch'))

# remove off-target samples
rows2remove <- c(rows2remove, which(allevents0$gene=='off'))

# remove low coverage samples
lowcovrows <- c( which(allevents0$seqnames=='slc24a5_AB_2') , which(allevents0$seqnames=='off_AB1_4') , 
                 which(allevents0$seqnames=='tyr_AA_3') , which(allevents0$seqnames=='scn1lab_AC_2'), 
                 which(allevents0$seqnames=='scn1lab_AC_3'), which(allevents0$seqnames=='scn1lab_AC_4') )
rows2remove <- c(rows2remove, lowcovrows)

length(unique(rows2remove))
nrow(allevents0)

# 125,811 rows - 102,558 rows = 23,253 rows
allevents <- allevents0[-rows2remove,]

# >>> I think I do not need to remove duplicates here
# as will do histogram per sample; so differences in coverage is not an issue

# if del >> width should be negative ---------------------------------------

allevents$width[which(allevents$type=='deletion')] <- all.neg(allevents$width[which(allevents$type=='deletion')])
# looks good
# could add a line to check

# unfold allevents --------------------------------------------------------
# i.e. when counts > 1, copy/paste multiple times the row
# essentially: so counts can only be 1

# avoid running again; takes at least ~ 30 minutes

fold <- as.data.frame(matrix(ncol=ncol(allevents)))
colnames(fold) <- colnames(allevents)

for (rown in 1:nrow(allevents)) { # for each row
  
  cat('Row', rown, 'out of', nrow(allevents), '\n')

  rowx <- allevents[rown,] # store the line
  
  if (rowx$counts > 1) { # if counts is > 1
    times <- rowx$counts - 1 # ! will already have one in allevents; so just add counts - 1
    rowx$counts <- 1 # can replace counts by 1 as 'unfolded'
    
    for (t in 1:times) {
      fold <- rbind(fold, rowx) # each 'times' times the line to fold
    }
  }
}

# first row of fold is all NA
fold <- fold[-1,]

# check all counts in fold are 1
# then will rbind 'fold' to all events
# then can replace all the counts by 1

# number of rows should be same as total of counts before


# different checks --------------------------------------------------------

# all counts in fold should be 1
length(which(fold$counts != 1)) # 0; ok

# number of rows in fold should be sum(counts) - nrow(counts)
nrow(fold) == sum(allevents$counts) - nrow(allevents) # ok


# merge fold with allevents ----------------------------------------------

# can first replace all counts in allevents by 1
allevents$counts <- 1

alluf <- rbind(allevents, fold)

alluf <- alluf[order(alluf$seqnames),]

# export alluf ------------------------------------------------------------------
# better to export and import again in other scripts as takes a long time

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
write.csv(alluf, 'events_all_unfolded.csv', row.names=FALSE)