# MiSeq amplican: quantify diversity in null alleles

library(dplyr)
library(ggplot2)

# import ------------------------------------------------------------------

myfile <- '~/Dropbox/phd/septMiSeq/septMiSeq_amplican/events_histo/events_all_unfolded.csv'
alluf0 <- read.csv (myfile)
# I do not actually need the unfolded version as I will look at unique variant; but less filtering to do

# remove some rows --------------------------------------------------------
rows2remove <- c()

# remove control samples
rows2remove <- c(rows2remove, which(alluf0$fish == 0))

# remove slc24a5 AC samples
rows2remove <- c(rows2remove, which(alluf0$locus == 'slc24a5_AC'))

rows2remove <- unique(rows2remove)
length(rows2remove) # 17233 rows to be removed

# remove these rows
alluf <- alluf0[- rows2remove ,]

# off-targets, low-coverage, consensus FALSE, mismatch all removed already

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
# 102,963; same as nrow(allevents) because for each line: duplicated TRUE or FALSE
# duplicated says FALSE first instance; then TRUE next instances, so will keep one copy of all elements

duprows <- alluf[which(dups),] # 102,963 rows duplicates
uni <- alluf[-which(dups),] # 6971 left; looks ok

# count number of unique mutations per locus ------------------------------

uni <- tbl_df(uni)

nmuts <- uni %>% # n mutations
  group_by(gene, locus, fish, seqnames) %>% 
  tally()

mean(nmuts$n) # 58
sd(nmuts$n) # 45.7
min(nmuts$n) # 3
max(nmuts$n) # 194


# w/ read_id --------------------------------------------------------------
# above estimate might be overevaluating
# because: if one read has 2 indels; it counts 2 (as they are both unique events)
# but it is really 1 allele, as it is the same copy of the genome

# how often does this happen
# I think simply looking at read_id in uni (col 12)
# ! I am guessing read_id starts back to 1 for each sample; so look for duplicated read_id and seqnames
readdups <- duplicated(uni[, c(1, 2, 3, 4, 12)])

# eg. row 4 vs row 49
# csnk1db AA 1
# row4: deletion -25 to -7; read id 78
# row49: insertion -3 to -1; read id 78
# >> good example; should be counted ONCE; not twice

length(which(readdups)) # 262 / 6971
# >> not a massive issue
# but easily solvable

# I can simply discard the duplicate TRUE; as I am simply counting them
# If a read has this indel alone; its read_id should be different; so should not create any issues

uniuni <- uni[-which(readdups),]

# count number of unique mutations per locus ------------------------------

uniuni <- tbl_df(uniuni)

nmuts <- uniuni %>% # n mutations
  group_by(gene, locus, fish, seqnames) %>% 
  tally()

ggplot(data=nmuts, aes(x=n)) +
  geom_density()

mean(nmuts$n)
sd(nmuts$n)

median(nmuts$n)
mad(nmuts$n)
# !!! R's mad function is not computing what you think it does
# confusing; check this to understand where it comes from

# manually calculating MAD:
# 1) residuals absolute value of values - median
abs(nmuts$n - median(nmuts$n))
# 2) median of the residuals
median(abs(nmuts$n - median(nmuts$n)))

min(nmuts$n)
max(nmuts$n)

length(which(nmuts$n < 50))
nrow(nmuts)

# maybe a safe solution to confirm would be to see if I get similar estimates by eye in IGV
# >> checked 3 samples; it looks good

# number of possible gene-wide genotypes if random
# i.e. combinatorial

# it should be a simple multiplication
# for instance:
# locus A has 2 mutations / locus B has 3 mutations >> number of possible genotypes = 6
# locus A has 2 mutations / locus B has 3 mutations / locus C has 2 mutations >> number of possible genotypes = 12

nmuts <- as.data.frame(nmuts)

# remove 4th locus slc24a5 / tyr
nmuts <- nmuts[- which(nmuts$locus=='tyr_AD' | nmuts$locus=='slc24a5_AG')]

nmutsdiv <- nmuts %>% # number of mutations/diversity
  group_by(gene, fish) %>%
  summarise_at(vars(n),
               prod)

mean(nmutsdiv$n)
sd(nmutsdiv$n)

# sd > mean because distribution is very right-skewed
ggplot(data=nmutsdiv, aes(x=n)) +
  geom_density()

min(nmutsdiv$n)
max(nmutsdiv$n)