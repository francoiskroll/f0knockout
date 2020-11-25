# scn1lab comparisons/time spent active
# for eLife reviews

# nomenclature for experiments

  # 1. F0 exp1 = 190813_07
  # 2. F0 exp2 = 190813_06
  
  # 3. stable exp1 = scn1lab del44 = 140218_34
  # 4. stable exp2 = scn1lab didy = 180528_0A

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

msgspace <- c('\n \t \t \t \t \t \t \t \t \t \t \t \t \t \t') # to use when printing messages

# import ------------------------------------------------------------------

where <- '~/.../'

# paths to DATA files
f01_path <- paste(where, '190813_07_DATA.txt', sep='')
f02_path <- paste(where, '190813_06_DATA.txt', sep='')
st1_path <- paste(where, '140218_34_DATA.txt', sep='')
st2_path <- paste(where, '180528_0A_DATA.txt', sep='')

# paths to genotype files
f01_genopath <- paste(where, '190813_07genotype.txt', sep='')
f02_genopath <- paste(where, '190813_06genotype.txt', sep='')
st1_genopath <- paste(where, '140218_34genotype.txt', sep='')
st2_genopath <- paste(where, '180528_0Agenotype.txt', sep='')

# import DATA files
f01 <- read.table(f01_path, header=TRUE, skip=1, sep='\t')
f02 <- read.table(f02_path, header=TRUE, skip=1, sep='\t')
st1 <- read.table(st1_path, header=TRUE, skip=1, sep='\t')
st2 <- read.table(st2_path, header=TRUE, skip=1, sep='\t')

# import genotype files
f01geno <- read.table(f01_genopath, header=TRUE, skip=1, sep='\t')
f02geno <- read.table(f02_genopath, header=TRUE, skip=1, sep='\t')
st1geno <- read.table(st1_genopath, header=TRUE, skip=1, sep='\t')
st2geno <- read.table(st2_genopath, header=TRUE, skip=1, sep='\t')

# vector of exp names (useful later)
exps <- c('190813_07',
          '190813_06',
          '140218_34',
          '180528_0A')

# clean up  ---------------------------------------------------------------

# put DATAs in a list
Mid <- list(f01, f02, st1, st2) # Mid for middurs (the data contained in each DATA file = time spent active)
Gen <- list(f01geno, f02geno, st1geno, st2geno) # Gen for genotypes

# pad Gen with NA rows so they all have 96 rows (will help later when adding the empty wells as a genotype)
padWithNArowsUntil <- function(x, until) { # x = dataframe, until = pad until how many rows
  padding <- as.data.frame(matrix(nrow=until-nrow(x), ncol=ncol(x)))
  colnames(padding) <- colnames(x) # if colnames are different, rbind freaks out
  x <- rbind(x, padding)
  return(x)
}

Gen <- lapply(Gen, padWithNArowsUntil, until=96) # all dataframes in Gen will now be 96 rows (filled in with NA)

# extract each CLOCK column (it is the last one) from Mid dataframes (df) and store elsewhere
Clock <- lapply(Mid, function(x) {x[,ncol(x)]}) # list of Clock values (each is one big vector)

# each df (dataframe) in Mid has 2 extra columns: start, end (do not need them -- delete)
# and delete last 3 columns (2 NA columns + 1 CLOCK column)
Mid <- lapply(Mid, function(x) {x[,3:(ncol(x)-3)]})

# check all 96 columns
if (! identical(as.numeric(sapply(Mid, ncol)) , rep(96, length(Mid)))) stop('>>> ERROR: Not all DATA are 96 columns')

# improve column names
colfish <- sprintf('f%i', 1:96) # from above, we are sure all df have 96 columns
Mid <- lapply(Mid, setNames, colfish)

# will keep list numbered as above (1 to 6)
# (instead of changing name of element)


# sum data in each window -------------------------------------------------

Midsum <- lapply(1:length(Mid), function(x) as.data.frame(matrix(nrow=4, ncol=96))) # preallocate list of empty dataframes

# structure should be:
  # \start\
    # DAY0 -- do not use
  # sunset1
    # NIGHT0 -- do not use
  # sunrise1
    # DAY1
  # sunset2
    # NIGHT1
  # sunrise2
    # DAY2
  # sunset3
    # NIGHT2
  # sunrise3
    # DAY3 -- do not use
  # \stop\

# ! for files where the CLOCK value was manually corrected; Closest to 14 will only return 1 value due to tiny rounding differences
  # (while expect multiple for the different days)
  # solution: look within a very small interval just after 14 or 0 (for 24H)
# (see below)

# will also check how many rows of data there are in each timewindow to make sure it is the same between exp (or close)
windowCheck <- as.data.frame(matrix(nrow=4, ncol=length(Mid)))

for (exp in 1:length(Mid)) { # for each df in Mid
  
  # get sunset/sunrise row numbers
  se <- which(Clock[[exp]] >= 14 & Clock[[exp]] < 14+0.015) # find the sunsets = 14H after 9am = 11pm
  cat(msgspace, '>> Found', length(se), 'sunsets\n')
  si <- which(Clock[[exp]] >= 0 & Clock[[exp]] < 0+0.015) # find the sunrises = 24H after 9am = 9am next day
  cat(msgspace, '>> Found', length(si), 'sunrises\n')
  
  # fill in Midsum with sum of data within these windows
  Midsum[[exp]][1,] <- colSums(Mid[[exp]][c((si[1]+1): se[2]),]) # DAY1 = from sunrise1 to sunset2
  windowCheck[1, exp] <- length((si[1]+1): se[2])
  
  Midsum[[exp]][2,] <- colSums(Mid[[exp]][c((si[2]+1): se[3]),]) # DAY2 = from sunrise2 to sunset3
  windowCheck[2, exp] <- length((si[2]+1): se[3])
  
  Midsum[[exp]][3,] <- colSums(Mid[[exp]][c((se[2]+1): si[2]),]) # NIGHT1 = from sunset2 to sunrise2
  windowCheck[3, exp] <- length((se[2]+1): si[2])
  
  Midsum[[exp]][4,] <- colSums(Mid[[exp]][c((se[3]+1): si[3]),]) # NIGHT2 = from sunset3 to sunrise3
  windowCheck[4, exp] <- length((se[3]+1): si[3])

}

# windowCheck OK
  # each day = 840 min = 14H
  # each night = 600 min = 10H

# I think topup is not too much of a worry as gets normalised within experiment to Wt by Z-score
# (I also don't see a solution that is significantly better than leaving it there)

# convert Midsum to Z-score vs Wt -----------------------------------------

MidsumZ <- lapply(1:length(Mid), function(x) as.data.frame(matrix(nrow=4, ncol=96))) # preallocate list of empty dataframes

for (exp in 1:length(Mid)) {
  
  Wtmean <- apply(Midsum[[exp]][, Gen[[exp]]$wt[!is.na(as.numeric(Gen[[exp]]$wt))] ], # takes the columns of the Wt fish
                  1,
                  mean) # mean of their total middur
  
  Wtsd <- apply(Midsum[[exp]][, Gen[[exp]]$wt[!is.na(as.numeric(Gen[[exp]]$wt))] ], # takes the columns of the Wt fish
                1,
                sd) # mean of their total middur
  
  for (win in 1:length(Wtmean)) {
    MidsumZ[[exp]][win,] <- (Midsum[[exp]][win,] - Wtmean[win]) / Wtsd[win] # Zscore is (x - mean) / sd
  }
  
}


# add identifier ----------------------------------------------------------

# will add unique ID as column name
  # exp_genotype_well

# before: add 'empty well' column to genotypes
for (exp in 1:length(Gen)) {
  empties <- c(1:96)[! 1:96 %in% sort(as.vector(unlist(Gen[[exp]])[!is.na(unlist(Gen[[exp]]))]))] # gives the empty wells
  # after %in%: vector of all wells from genotype files (w/o NA)
  # in 1, 2, 3, ..., 96 which ones are NOT (!) in genotype file wells >> these are the empty wells (i.e. the wells that are not mentioned)
  # >> vector 96 x TRUE (well is empty) or FALSE (well is not empty as was mentioned in genotype file)
  Gen[[exp]]$empty <- NA # first preallocate the column as empty
  Gen[[exp]]$empty[1:length(empties)] <- empties# then fill it in with the empty wells
}

# create/add identifier
for (exp in 1:length(Mid)) { # loops thru experiments
  
  idsexp <- vector()
  # columns in Mid are always fish 1:96, 
  # so can just make a big bag of IDs from the different genotype and will then order them as 1:96 to add them as column names
  
  for (gen in 1:ncol(Gen[[exp]])) { # loops through genotype

    # below builds id for the fish in that genotype
    idsexp <- c(idsexp, paste(exps[[exp]],
                              colnames(Gen[[exp]])[gen],
                              Gen[[exp]][!is.na(Gen[[exp]][,gen]),gen],
                              sep='_')) # append the identifiers for these fish
    
  }
  
  if(length(idsexp) != 96) stop('>>> Error: not all identifiers are 96')
  # now need to put the identifiers in the right order
  idsexp <- idsexp[order(as.numeric(sapply(idsexp, function(x) strsplit(x, '_')[[1]][4])))] # within [...]: gives the right order (ranking) to follow
  # takes the fish numbers of the identifier and order it
  
  # now can swap the column names in Mid
  colnames(MidsumZ[[exp]]) <- idsexp

}


# merge in one dataframe --------------------------------------------------

# as column have now unique identifier; can collapse in one big dataframe

midz <- do.call('cbind', MidsumZ)
midz$window <- c('day1', 'day2', 'night1', 'night2')
# midz is :
  # one column per fish, i.e. 96 fish * 6 experiments = 576 columns
  # one row per day/night = 4 rows = day1/day2/night1/night2 (in that order)

midz <- tbl_df(midz)

midzl <- midz %>% # midz long
  pivot_longer(-window, names_to='id', values_to='zmid')

# split the id column
midzl$exp <- paste(sapply(midzl$id, function(x) strsplit(x, '_')[[1]][1]),
                  sapply(midzl$id, function(x) strsplit(x, '_')[[1]][2]),
                  sep='_')

midzl$gen <- as.vector(sapply(midzl$id, function(x) strsplit(x, '_')[[1]][3]))

midzl$well <- as.numeric(sapply(midzl$id, function(x) strsplit(x, '_')[[1]][4]))

# also need a column exp_gen for colour on plot
midzl$expgen <- paste(midzl$exp, midzl$gen, sep='_')

# put columns in nice order: id, expgen, exp, gen, well, window, zmid
midzl <- midzl[,c(2, 7, 4, 5, 6, 1, 3)]

# fix factor levels for plot
midzl$exp <- factor(midzl$exp, levels=exps) # put exps in order they were uploaded in
midzl$gen <- factor(midzl$gen, levels=c('empty', 'wt', 'het', 'hom'))
midzl$expgen <- factor(midzl$expgen,
                       levels=c('190813_07_empty', '190813_07_wt', '190813_07_hom',
                                '190813_06_empty', '190813_06_wt', '190813_06_hom',
                                '140218_34_empty', '140218_34_wt', '140218_34_het', '140218_34_hom',
                                '180528_0A_empty', '180528_0A_wt', '180528_0A_het', '180528_0A_hom'))

# plot DAY --------------------------------------------------------------------

# mycolsWt <- c('#697a87', '#74210b', # F0exp1: wt / hom
#               '#697a87', '#c23712', # F0exp2: wt / hom
#               '#697a87', '#b3bdc4', '#f1876b', # stableexp1: wt/het/hom
#               '#697a87', '#b3bdc4', '#f1876b', # stableexp2
#               '#697a87', '#b3bdc4', '#f1876b', # stableexp3
#               '#697a87', '#b3bdc4', '#f1876b') # stableexp4

mycolsWt <- c('#697a87', '#74210b', # F0exp1: wt / hom
              '#697a87', '#c23712', # F0exp2: wt / hom
              '#697a87', '#b3bdc4', '#f1876b', # stableexp1: wt/het/hom
              '#697a87', '#b3bdc4', '#f1876b') # stableexp2

# mycols <- c('#74210b', # F0exp1: hom
#             '#c23712', # F0exp2: hom
#             '#b3bdc4', '#f1876b', # stableexp1: het/hom
#             '#b3bdc4', '#f1876b', # stableexp2
#             '#b3bdc4', '#f1876b', # stableexp3
#             '#b3bdc4', '#f1876b') # stableexp4

mycols <- c('#74210b', # F0exp1: hom
            '#c23712', # F0exp2: hom
            '#b3bdc4', '#f1876b', # stableexp1: het/hom
            '#b3bdc4', '#f1876b') # stableexp2

# WITH WT
midzlDayWt <- subset(midzl, (window=='day1' | window=='day2') & # plot only days
                   (gen=='wt' | gen=='het' | gen=='hom')) # plot only wt, het, hom (not empty)

midDayWt <- ggplot(data=midzlDayWt, aes(x=exp, y=zmid, colour=expgen)) +
  geom_quasirandom(dodge.width=0.6, varwidth=TRUE, size=0.2) + # dots like they were inside a violin plot + dodge_width splits by genotype
  stat_summary(fun.y=mean, aes(group=expgen), position=position_dodge(width=0.6), geom='point', colour='#020304', shape=3, size=1, stroke=1) +
  # ! position_dodge(width=) should be same as dodge.width above
  scale_colour_manual(values=mycolsWt) +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.title.x=element_blank(),
    legend.position='none') +
  ylab('') + # will add in Illustrator so around 0
  coord_cartesian(ylim=c(-3.2, 4.0)) +
  scale_x_discrete(labels=c('F0 exp1', 'F0 exp2', 'stable exp1', 'stable exp2')) +
  scale_y_continuous(breaks=c(-3, 0, 3))
midDayWt
# plotting with Wt allows to check that the Wt mean is always exactly 0
  # (by z-score normalisation)

# WITHOUT WT
midzlDay <- subset(midzl, (window=='day1' | window=='day2') & # plot only days
                  (gen=='het' | gen=='hom')) # plot only het, hom (not wt or empty)

midDay <- ggplot(data=midzlDay, aes(x=exp, y=zmid, colour=expgen)) +
  geom_quasirandom(dodge.width=0.6, varwidth=TRUE, size=0.2) + # dots like they were inside a violin plot + dodge_width splits by genotype
  stat_summary(fun.y=mean, aes(group=expgen), position=position_dodge(width=0.6), geom='point', colour='#020304', shape=3, size=1, stroke=1) +
  # ! position_dodge(width=) should be same as dodge.width above
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.title.x=element_blank(),
    legend.position='none') +
  ylab('') + # will add in Illustrator so around 0
  coord_cartesian(ylim=c(-3.2, 4.0)) +
  scale_x_discrete(labels=c('F0 exp1', 'F0 exp2', 'stable exp1', 'stable exp2', 'stable exp3', 'stable exp4')) +
  scale_y_continuous(breaks=c(-3, 0, 3))
midDay

# plot NIGHT --------------------------------------------------------------

# with Wt
midzlNightWt <- subset(midzl, (window=='night1' | window=='night2') & # plot only days
                         (gen=='wt' | gen=='het' | gen=='hom')) # plot only wt, het, hom (not empty)

midNightWt <- ggplot(data=midzlNightWt, aes(x=exp, y=zmid, colour=expgen)) +
  geom_quasirandom(dodge.width=0.6, varwidth=TRUE, size=0.2) + # dots like they were inside a violin plot + dodge_width splits by genotype
  stat_summary(fun.y=mean, aes(group=expgen), position=position_dodge(width=0.6), geom='point', colour='#020304', shape=3, size=1, stroke=1) +
  # ! position_dodge(width=) should be same as dodge.width above
  scale_colour_manual(values=mycolsWt) +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.title.x=element_blank(),
    legend.position='none') +
  ylab('') + # will add in Illustrator so around 0
  coord_cartesian(ylim=c(-2.0, 10)) +
  scale_x_discrete(labels=c('F0 exp1', 'F0 exp2', 'stable exp1', 'stable exp2', 'stable exp3', 'stable exp4')) +
  scale_y_continuous(breaks=c(-2, 0, 2))
midNightWt

# without Wt
midzlNight <- subset(midzl, (window=='night1' | window=='night2') & # plot only days
                       (gen=='het' | gen=='hom')) # plot only wt, het, hom (not empty)

midNight <- ggplot(data=midzlNight, aes(x=exp, y=zmid, colour=expgen)) +
  geom_quasirandom(dodge.width=0.6, varwidth=TRUE, size=0.2) + # dots like they were inside a violin plot + dodge_width splits by genotype
  stat_summary(fun.y=mean, aes(group=expgen), position=position_dodge(width=0.6), geom='point', colour='#020304', shape=3, size=1, stroke=1) +
  # ! position_dodge(width=) should be same as dodge.width above
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.x=element_blank(), # easier to add the genotype names in Illustrator
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.title.x=element_blank(),
    legend.position='none') +
  ylab('') + # will add in Illustrator so around 0
  coord_cartesian(ylim=c(-2.0, 10)) +
  scale_x_discrete(labels=c('F0 exp1', 'F0 exp2', 'stable exp1', 'stable exp2', 'stable exp3', 'stable exp4')) +
  scale_y_continuous(breaks=c(-2, 0, 2))
midNight

# >> could only plot second Half of night but seems arbitrary; I think better to explain


# statistics --------------------------------------------------------------

# DAY
summary(aov(zmid ~ expgen, midzlDayWt)) # p ~ 0

pairwise.t.test(midzlDayWt$zmid, midzlDayWt$expgen)
# the p-val adjustment is slightly too stringent as I am not looking at all these comparisons
# but fine to leave it (as more stringent than necessary)

# p-values I am interested in

# F0 exp1 (box7): wt vs hom = < 2e-16 ***
# F0 exp2 (box6): wt vs hom = 2.0e-07 ***

# del44 (140218): wt vs het = 1.0 ns
# del44 (140218): wt vs hom = 1.5e-09 ***

# didy (180528): wt vs het = 0.02126 *
# didy (180528): wt vs hom = 0.00167 **


# NIGHT
summary(aov(zmid ~ expgen, midzlNightWt)) # p ~ 0

pairwise.t.test(midzlNightWt$zmid, midzlNightWt$expgen)
# the p-val adjustment is slightly too stringent as I am not looking at all these comparisons
# but fine to leave it (as more stringent than necessary)

# p-values I am interested in

# F0 exp1 (box7): wt vs hom = 1.0 ns
# F0 exp2 (box6): wt vs hom = 0.00388 **

# del44 (140218): wt vs het = 0.18642 ns
# del44 (140218): wt vs hom = 1.0 ns

# didy (180528): wt vs het = 1.0 ns
# didy (180528): wt vs hom = 4.1e-12 ***


# Ns ----------------------------------------------------------------------

length(unique(subset(midzlNightWt, expgen=='190813_07_wt')$well))
length(unique(subset(midzlNightWt, expgen=='190813_07_hom')$well))

length(unique(subset(midzlNightWt, expgen=='190813_06_wt')$well))
length(unique(subset(midzlNightWt, expgen=='190813_06_hom')$well))

length(unique(subset(midzlNightWt, expgen=='140218_34_wt')$well))
length(unique(subset(midzlNightWt, expgen=='140218_34_het')$well))
length(unique(subset(midzlNightWt, expgen=='140218_34_hom')$well))

length(unique(subset(midzlNightWt, expgen=='180528_0A_wt')$well))
length(unique(subset(midzlNightWt, expgen=='180528_0A_het')$well))
length(unique(subset(midzlNightWt, expgen=='180528_0A_hom')$well))

# check day is the same
length(unique(subset(midzlDayWt, expgen=='190813_07_wt')$well))
length(unique(subset(midzlDayWt, expgen=='190813_07_hom')$well))

length(unique(subset(midzlDayWt, expgen=='190813_06_wt')$well))
length(unique(subset(midzlDayWt, expgen=='190813_06_hom')$well))

length(unique(subset(midzlDayWt, expgen=='140218_34_wt')$well))
length(unique(subset(midzlDayWt, expgen=='140218_34_het')$well))
length(unique(subset(midzlDayWt, expgen=='140218_34_hom')$well))

length(unique(subset(midzlDayWt, expgen=='180528_0A_wt')$well))
length(unique(subset(midzlDayWt, expgen=='180528_0A_het')$well))
length(unique(subset(midzlDayWt, expgen=='180528_0A_hom')$well))

# OK

# export DAY ------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='DayWt.pdf', plot=midDayWt, width=85, height=90, units='mm', useDingbats=FALSE) # f0 paper review

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='Day.pdf', plot=midDay, width=85, height=90, units='mm', useDingbats=FALSE) # f0 paper review

# export NIGHT ------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='NightWt.pdf', plot=midNightWt, width=85, height=90, units='mm', useDingbats=FALSE) # f0 paper review

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='Night.pdf', plot=midNight, width=85, height=90, units='mm', useDingbats=FALSE) # f0 paper review