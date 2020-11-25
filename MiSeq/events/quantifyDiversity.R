# MiSeq amplican: quantify diversity in null alleles

library(dplyr)
library(ggplot2)
library(openxlsx)

# import ------------------------------------------------------------------

myfile <- '~/Dropbox/phd/septMiSeq/septMiSeq_amplican/events_histo/events_all_unfolded.csv'
alluf0 <- read.csv (myfile)

myfile <- '~/Dropbox/phd/septMiSeq/septMiSeq_amplican/septMiSeq_amplicanresults.xlsx'
amp <- read.xlsx (myfile)

# remove some rows --------------------------------------------------------
rows2remove <- c()

# remove control samples
rows2remove <- c(rows2remove, which(alluf0$fish == 0)) # 5895 rows

# remove slc24a5 AC samples
rows2remove <- c(rows2remove, which(alluf0$locus == 'slc24a5_AC'))

# remove negative strand
(nrow(subset(alluf, strand=='-')) / nrow(subset(alluf, strand=='+'))) * 100 # very few mutations on negative strand
# and seem to be duplicated on the positive strand
# not entirely sure what they represent;
# I think safer to delete all mutations called on negative strand, especially if counting coverage as proxy for allele frequency
rows2remove <- c(rows2remove, which(alluf0$strand == '-'))

rows2remove <- unique(rows2remove)
length(rows2remove)

# remove these rows
alluf <- alluf0[- rows2remove ,]

# off-target samples, low-coverage samples, consensus FALSE events, mismatch events all removed already

# only unique events ------------------------------------------------------
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

dups <- duplicated(alluf[, c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11)]) 
# 113,065; same as nrow(allevents) because for each line: duplicated TRUE or FALSE
# duplicated says FALSE first instance; then TRUE next instances, so will keep one copy of all elements

duprows <- alluf[which(dups),] # 106,050 rows duplicates
uni <- alluf[-which(dups),] # 7015 left; looks ok

# count number of unique mutations per locus ------------------------------
uni <- tbl_df(uni)

nmuts <- uni %>% # n mutations
  group_by(gene, locus, fish, seqnames) %>% 
  tally()

mean(nmuts$n) # 57.5
median(nmuts$n) # 42
sd(nmuts$n) # 45.8
min(nmuts$n) # 3
max(nmuts$n) # 192

# this estimate is likely an overevaluation
# if one read has 2 indels; it will count 2, as they are detected as unique events
# however, it is really 1 allele, as it is the same copy of the genome

# one way to check is to count number of unique read_id;
# each read_id should roughly correspond to one copy of the genome (assuming no sequencing error etc)
# does it give a similar estimate?
nids <- uni %>% # n read IDs
  group_by(gene, locus, fish, seqnames) %>% 
  summarise(nid=n_distinct(read_id))

mean(nids$nid) # 55
median(nids$nid) # 40.5
sd(nids$nid) # 43.9
min(nids$nid) # 3
max(nids$nid) # 185
# very close, which is great news
# >>> will use this estimate in text
# will use median ± median absolute deviation
mad(nmuts$n)
# !!! R's mad function is not computing what you think it does
# manually calculating MAD:
# 1) residuals absolute value of values - median
abs(nids$nid - median(nids$nid))
# 2) median of the residuals
median(abs(nids$nid - median(nids$nid)))
# 27.5
# so: 40.5 ± 27.5

ggplot(data=nids, aes(x=nid)) +
  geom_density()

# extra way to check further the estimate: delete the duplicated read_id and count again
iddups <- duplicated(uni[, c(1, 2, 3, 4, 12)]) # read ID duplicated
length(which(iddups)) # 261 / 7015
# discard the duplicate read_id TRUE
# If a read has this indel alone; its read_id should be different; so should not create any issues
uniuni <- uni[-which(iddups),]
uniuni <- tbl_df(uniuni)
nmutid <- uniuni %>% # n mutations
  group_by(gene, locus, fish, seqnames) %>% 
  tally()
sum(nmutid$n==nids$nid)==nrow(nmutid) # as expected, both methods are identical

# safe solution to confirm would be to see if I get similar estimates by eye in IGV
# >> checked 3 samples; it looks good

# number of possible gene-wide genotypes, if random
# i.e. combinatorial

# it should be a simple multiplication
# for instance:
# locus A has 2 mutations / locus B has 3 mutations >> number of possible genotypes = 6
# locus A has 2 mutations / locus B has 3 mutations / locus C has 2 mutations >> number of possible genotypes = 12

nids <- as.data.frame(nids)

# remove 4th locus slc24a5 / tyr
nids <- nids[- which(nids$locus=='tyr_AD' | nids$locus=='slc24a5_AG')]

nidsdiv <- nids %>% # number of mutations/diversity
  group_by(gene, fish) %>%
  summarise_at(vars(nid),
               prod)

mean(nidsdiv$nid)
sd(nidsdiv$nid)
min(nidsdiv$nid)
max(nidsdiv$nid)
# sd > mean because distribution is very right-skewed
ggplot(data=nmutsdiv, aes(x=n)) +
  geom_density()

# check correct with eg. csnk1db fish 1
prod(subset(nids, gene=='csnk1db' & fish==1, nid))
# OKAY

# frequency of alleles ----------------------------------------------------
# in events_all_unfolded: each mutation/read is one row
# so can simply take the number of rows for each unique mutation within each sample to estimate its frequency (i.e. abundance in the knockout fish)
# similar logic as allevents_histo.R
# interpretation slightly different:
# histogram = type of mutation Cas9 generates (this is why taking unique, goal is to get as close as possible to unique repair events)
# here = frequency of alleles in each animal

# as above, will do with unique mutation then read_id
# I think read_id is likely to be a better estimate, as explained above

# 1- unique mutations
alfreq <- alluf %>%
  group_by(seqnames,
           start,
           end,
           width,
           strand,
           originally,
           replacement,
           type) %>%
  summarise(nread=n())
# each row is one mutation, nread is number of reads carrying that mutation

alfreq$coverage <- NA

# get reads_filtered (i.e. coverage) from amplican results
for (r in 1:nrow(alfreq)){
  spl <- alfreq[r,'seqnames']
  spl <- as.character(unlist(spl))
  alfreq[r, 'coverage'] <- subset(amp, ID==spl, Reads_Filtered)
}

# calculate freq
alfreq$freq <- alfreq$nread/alfreq$coverage

alfreq[which.max(alfreq$freq),]
alfreq[which.min(alfreq$freq),]

# how many alleles below 10% frequency
nrow(subset(alfreq, freq<=0.1)) / nrow(alfreq) # 96% of mutations are at < 10% abundance in the sample

# 2- read_id
# now same with read_id to compare (I think probably better approximation of alleles)
# ! number of times each read_id is repeated does not work as an approximation of abundance
# because it can be repeated many times if single read has multiple indels, which does not represent more coverage of that read
# solution: you want a 1-to-1 match 1 read_id : 1 mutation;
# then delete all rows which mention mutations that are not in that 'simplified' version
alluf11 <- alluf %>% # alluf 1 read_id : 1 mutation
  group_by(gene, locus, fish, seqnames,
           read_id) %>% # group by read_id
  sample_n(1) # then take one row at random for each read_id
# check
nrow(alluf11)
# should be same as
nidsinsample <- alluf %>%
  group_by(gene, locus, fish, seqnames) %>%
  summarise(nid=n_distinct(read_id))
sum(nidsinsample$nid)
# correct

# now remove from alluf everything that is not mentioned in alluf11
# I think semi_join is what I need: return all rows of alluf where they are matching values in alluf11; keeping only columns from alluf
# and never duplicate rows of alluf
allufsing <- semi_join(alluf, alluf11, by=c('seqnames', 'start', 'end', 'width', 'strand', 'originally', 'replacement', 'type', 'read_id'))
# all events single = read_id can only have a single indel

anti_join(alluf, alluf11, by=c('seqnames', 'start', 'end', 'width', 'strand', 'originally', 'replacement', 'type', 'read_id')) # these are the rows removed

# check
nrow(alluf) - nrow(allufsing) # 5446 rows were removed
# example1
# this mutation is in alluf11
subset(alluf11, seqnames=='csnk1db_AA_1' & start==-5 & end==-4 & width==-2 & read_id==1)
# so it is also in alluf, but multiple times
subset(allufsing, seqnames=='csnk1db_AA_1' & start==-5 & end==-4 & width==-2 & read_id==1)

# example2
# this mutation was in original alluf
subset(alluf, seqnames=='csnk1db_AA_1' & start==-25 & end==-7 & type=='deletion' & read_id==78)
# but was not selectd in alluf11
subset(alluf11, seqnames=='csnk1db_AA_1' & start==-25 & end==-7 & type=='deletion' & read_id==78) # no match
# consequence: it is not kept in allufsing
subset(allufsing, seqnames=='csnk1db_AA_1' & start==-25 & end==-7 & type=='deletion' & read_id==78) # no match

# now working on allufsing: count number of times each read_id is repeated
# as proxy for abundance of that allele
alfreqid <- allufsing %>% # allele frequency based on read ID
  group_by(seqnames, read_id) %>%
  summarise(nread=n())

# get reads_filtered (i.e. coverage) from amplican results
for (r in 1:nrow(alfreqid)){
  spl <- alfreqid[r,'seqnames']
  spl <- as.character(unlist(spl))
  alfreqid[r, 'coverage'] <- subset(amp, ID==spl, Reads_Filtered)
}
# calculate freq
alfreqid$freq <- alfreqid$nread/alfreqid$coverage

alfreqid[which.max(alfreqid$freq),]
alfreqid[which.min(alfreqid$freq),]

# compare both
ggplot(data=alfreqid, aes(x=freq)) +
  geom_density() +
  coord_cartesian(xlim=c(0, 0.05))

ggplot(data=alfreq, aes(x=freq)) +
  geom_density() +
  coord_cartesian(xlim=c(0, 0.05))
# similar idea -- allele frequency from unique mutations has higher frequencies
# these may be because same mutation is on reads that are slightly different
# eg. 2 bp del which is the same on 2 reads, but each read is slightly different (eg. another SNP from a sequencing error on one of the read)
# >> read_id stays unique (i.e. more unique read_ids = lower frequency of each) but mutation is merged as above (i.e. less unique mutations = higher frequency of each)
# it is hard to tell which one represent best alleles

# are there some summary insights that are common for both?

# how many alleles below 10% frequency
nrow(subset(alfreq, freq<=0.1)) / nrow(alfreq)
nrow(subset(alfreqid, freq<=0.1)) / nrow(alfreqid)
# > 95% of alleles are present at frequencies lower than 10%

# plot it well for author response
# I think will use frequency from mutation
# matches better what can be concisely explained in text
# if one read has 2 indels: their frequencies should be the same (as taking coverage as total) so hopefully a fair estimate
mycols <- c('#5a6974')
alfreqbar <- ggplot(data=alfreq[order(alfreq$freq),], aes (x=nrow(alfreq):1, y=freq)) +
  geom_bar(stat='identity', colour=mycols, width=0) + # minimum width
  theme_minimal() +
  theme(
    # grid
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    # axis text
    axis.title.x=element_text(size=9, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=7),
    legend.position='none'
  ) +
  xlab('mutation') + ylab('frequency in individual animal')
alfreqbar

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave('f0_authorresponse_mutfreq.pdf', plot=alfreqbar, width=110, height=60, unit='mm', useDingbats=FALSE) # f0 paper