# MiSeq amplican: all events/how similar different samples are

library(reshape2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(schoolmath)
library(ggbeeswarm)

# get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# import ------------------------------------------------------------------

myfile <- '~/.../MiSeq/events/events_all_unfolded.csv'
alluf0 <- read.csv (myfile)

# remove some rows --------------------------------------------------------
rows2remove <- c()

# remove control samples
rows2remove <- c(rows2remove, which(alluf0$fish == 0))

# remove slc24a5 AC samples
rows2remove <- c(rows2remove, which(alluf0$locus == 'slc24a5_AC'))

# remove negative strand
(nrow(subset(alluf0, strand=='-')) / nrow(subset(alluf0, strand=='+'))) * 100 # very few mutations are on negative strand
# and seem to be duplicated on the positive strand
# not sure what they represent;
# I think safer to delete mutations called on negative strand
rows2remove <- c(rows2remove, which(alluf0$strand == '-'))

rows2remove <- unique(rows2remove)
length(rows2remove)

# remove these rows
alluf <- alluf0[- rows2remove ,]

# off-targets, low-coverage, consensus FALSE, mismatch all removed already

# create composite ID column ---------------------------------------------
# i.e. one ID per mutation
alluf$mutID <- paste(alluf$start, alluf$end, alluf$width, alluf$originally, alluf$replacement, alluf$type, sep='_')

# frequency/probability of each indel for each sample -------------------------------------------

# will store each dataframe in a list
# so length of the list will be number of samples
length(unique(alluf$seqnames)) # 138 samples
ppw_all <- vector(mode='list', length=length(unique(alluf$seqnames)))
names(ppw_all) <- unique(alluf$seqnames)

for (spl in unique(alluf$seqnames)) {
  
  alluf_spl <- alluf[which(alluf$seqnames==spl),] # take the subset of alluf regarding that sample
  
  ppw <- aggregate(alluf_spl$mutID, by=list(alluf_spl$mutID), FUN=length)
  colnames(ppw) <- c('mutID', 'count')
  ppw$prob <- ppw$count / sum(ppw$count)
  if (round(sum(ppw$prob), 6) != 1) stop(cat('total probability for sample', spl, 'is not 1')) 
  # need to round, sometimes weird thing because of how it stores numbers
  
  ppw_all[[spl]] <- ppw
  
}

# pick top10 mutation per sample -----------------------------------

# loop thru the list
# pick the top 10 mutation ID based on P
# store them in a df sample x top1-10

# first for one element of the list

topn <- 10

topall <- as.data.frame(matrix(nrow=length(ppw_all), ncol=4+topn)) # 4 metadata about sample + n top indels

for (spln in 1:length(ppw_all)) { # for sample n
  ppw <-  ppw_all[[spln]] # take the dataframe at that position
  ppw <- ppw[order(ppw$prob, decreasing=TRUE),] # sort the mutation IDs based on probabilities (or counts, would do the same ranking)
  topnw <- ppw[1:topn,1] # top n mutation ID
  meta <- alluf[which(alluf$seqnames == names(ppw_all)[spln])[1], 1:4] # take the metada of that sample
  meta <- c( as.character(unlist(meta[1:2])) , as.integer(unlist(meta[3])),  as.character(unlist(meta[4])) ) # dirty fix to conversion problem
  rowtopnw <- c(meta, topnw) # stick the row together: metadata + top n mutation ID
  
  topall[spln,] <- rowtopnw # adding that row to the dataframe
}

colnames(topall) <- c('gene', 'locus', 'fish', 'seqnames', sprintf('top%s', 1:topn))

# improve formats, helps after
topall[5:(5+topn-1)] <- apply(topall[5:(5+topn-1)], 2, as.character)

# pairwise intersects -------------------------------------------------------------

# NA (when less than n indels in total) would affect a bit the comparisons...

# how many samples have less than 10 variants total
length( which( (rowSums(!is.na(topall)) - 4) < 10 )) # 8 samples have less than 10 variants
length( which( (rowSums(!is.na(topall)) - 4) < 5 )) # 2 samples have less than 5 variants

# remove samples with less than topn indels ------------------------------
topall <- topall [ - which(apply(topall[5:(5+topn-1)], 1, function(x) sum(is.na(x))) > 0) , ]
# inside gives the number of NA per row; then which one have more than 0 NA

# check
nrow(topall) # 114 samples = 122 - 8; ok


# pairwise comparisons between samples

# eg. for just two rows
length(intersect(as.character(topall[1,5:(5+topn-1)]) , as.character(topall[2,5:(5+topn-1)])))

# pairwise comparison > store how many wifths are same
topallpw <- as.data.frame(matrix(nrow=nrow(topall), ncol=(nrow(topall)))) # preallocate topall pairwise

for (i in 1:nrow(topall)) { # will loop thru rows of topall
  rowi <- as.character(topall[i, 5:(5+topn-1)]) # take the top10 indels
  
  for (j in 1:nrow(topall)) { # will loop again thru rows of topall to find intersects
    rowj <- as.character(topall[j, 5:(5+topn-1)])
    topallpw[i, j] <- length(intersect(rowi, rowj))
  }
}
# >> topallpw is symmetrical and main diagonal is all 10 so looks accurate

# take only upper half; as symmetrical would not want to repeat dots
topallpw <- get_upper_tri(topallpw) # looks correct

# columns & rows (first column; do not use rownames) are samples
topallpw <- cbind(topall$seqnames, topallpw)
colnames(topallpw) <- c('sample', topall$seqnames)

# pivot longer -------------------------------------------------

topallpw <- tbl_df(topallpw)

topallpwm <- topallpw %>% 
  pivot_longer(-sample, names_to='v_sample', values_to='int')
# result is three columns: sample; versus sample; number of common indels within top10

# check
nrow(topallpw) * (ncol(topallpw) - 1) == nrow(topallpwm) # expected number of rows

# is heatmap interesting? -------------------------------------------------

intsheat <- ggplot(topallpwm, aes(x=sample, y=v_sample, fill=int)) + 
  geom_tile(colour=NA) # edges around tiles
intsheat
# might be interesting, but too many comparisons

# pool in same locus vs different loci ------------------------------------

# add metadata
topallpwm$locus <- sub('_([^_]*)$', '', topallpwm$sample)
topallpwm$v_locus <- sub('_([^_]*)$', '', topallpwm$v_sample)

topallpwm <- topallpwm[, c(1, 4, 2, 5, 3)]

# add a column for the type of comparison
topallpwm$compar <- NA

# fill in that column
# there are three types of comparisons

# 1/ within sample (not interesting)
topallpwm$compar[topallpwm$sample==topallpwm$v_sample] <- 'same_sample'
# check
check1 <- nrow(topallpwm[topallpwm$sample==topallpwm$v_sample,])

# 2/ locus is the same; but sample is different
topallpwm$compar[ topallpwm$locus==topallpwm$v_locus & topallpwm$sample!=topallpwm$v_sample] <- 'same_locus'
# check
check2 <- nrow(topallpwm[topallpwm$locus==topallpwm$v_locus & topallpwm$sample!=topallpwm$v_sample,])

# 3/ locus is different
topallpwm$compar[ topallpwm$locus!=topallpwm$v_locus ] <- 'diff_locus'
# check
check3 <- nrow(topallpwm[topallpwm$locus!=topallpwm$v_locus,])

# final checks
if(check1 + check2 + check3 != nrow(topallpwm)) stop('problem') # are all rows falling in one of these categories
if(sum(is.na(topallpwm$compar)) != 0) stop('problem') # is there any rows that was not filled in

# convert to factor for plot
topallpwm$compar <- as.factor(topallpwm$compar)
# I think order is good

# plot ------------------------------------------------------------

# dot plot is more accurate as no jitter; but difficult to avoid overplotting

# ggplot(topallpwm, aes(x=compar, y=int, fill=compar)) +
#   geom_dotplot(binaxis='y', stackdir='center',
#                stackratio=0.3, dotsize=0.3, stroke=0.3, alpha=0.5)

# here with jitter
mycols <- c('#708490', '#4b5860')

# adding a bit of jitter to Y axis
topallpwm_jit <- topallpwm
topallpwm_jit$int <- jitter(topallpwm$int)

compardot <- ggplot(topallpwm_jit[topallpwm_jit$sample!=topallpwm_jit$v_sample ,], aes(x=compar, y=int)) + # data = only rows where sample is different
  geom_quasirandom(aes(colour=compar), stat='identity', groupOnX=TRUE, width=0.4, size=0.0, shape=20) +
  stat_summary(fun.y=mean, geom='point', colour='#020304', shape=3, size=2, stroke=1) +
  theme_minimal() +
  theme(
    panel.grid.minor=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size=7),
    legend.position='none'
  ) +
  scale_colour_manual(values=mycols) +
  scale_y_continuous(breaks=0:10) +
  scale_x_discrete(labels=c('different locus', 'same locus')) +
  coord_cartesian(ylim=c(0, 10)) +
  ylab('number of indels in common')
compardot

# summary stat ------------------------------------------------------------

topallpwm %>%
  group_by(compar) %>% # group by comparison type
  summarise_at(vars(int), # then summarise the int values, computing mean, sd and median
               list(
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 median= ~ median(., na.rm=TRUE),
                 ns=~ sum(!is.na(.)))) # ! do not count NA (i.e. reciprocal intersection; X vs Y then Y vs X; one was set to NA above)

t.test(int ~ compar, data = topallpwm[topallpwm$sample!=topallpwm$v_sample ,])
# p-val = 0

# export ------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='f0_top10indel.pdf', plot=compardot, width=90, height=60, units='mm', useDingbats=FALSE)