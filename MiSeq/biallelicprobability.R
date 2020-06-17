# septMiSeq amplican: biallelic probability plot

# similar plots as for theoretical model
# but can now do it per gene per fish

library (openxlsx)
library (ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# import ------------------------------------------------------------------

myfile <- '~/.../MiSeq_amplicanresults.xlsx'
summ0 <- read.xlsx (myfile)

# remove some samples ---------------------------------------------

rows2remove <- c()

# 1/ remove samples whose coverage is below 30X
  # 30x is probably a good threshold
  # it sounds like the correct way to report coverage for paired-end reads is to NOT double count (so 1 pair of reads = 1X)
  # i.e. what ampliCan does, but not IGV

rows2remove <- c(rows2remove, which(summ0$Reads_Filtered < 30))

# 2/ remove off-target samples
  # 44 samples: sounds correct (5 samples x 3 off-targets x 3 guides; minus 1 sample low coverage)
rows2remove <- c(rows2remove, which(summ0$ID=='off_AA1_0') : which(summ0$ID=='off_AD3_4'))

# 3/ remove slc24a5 AC samples
rows2remove <- c(rows2remove, which(summ0$crRNA_name=='slc24a5_AC'))

# 4/ remove control samples
rows2remove <- c(rows2remove, which(summ0$Control==TRUE)) # Control columns makes it easy to delete control fish for now

# >> remove the rows
rows2remove <- sort(unique(rows2remove)) # total of 102
summ <- summ0[-rows2remove,]

rownames(summ) <- NULL # resets row numbers

# compute Ps ---------------------------------------------------------------

# total number of reads is Reads_Filtered
# reads cannot be in two categories at once; barplot adds it on top
# need to compute reads not frameshifted, and add that category on top of reads frameshifted

# below can be shorter , but choosing clarity

edit <- (summ$Reads_Edited/summ$Reads_Filtered) # proportion of filtered reads that have edits
frameshift <- (summ$Reads_Frameshifted/summ$Reads_Filtered) # proportion of filtered reads that have frameshift
editbutnotfshift <- ((summ$Reads_Edited - summ$Reads_Frameshifted) / summ$Reads_Filtered) # essentially proportion of reads that have indels of lengths multiple of 3

# computing proportion of reads that do not lead to knockout, i.e. wildtype + edited but not frameshifted
readsWt <- summ$Reads_Filtered - summ$Reads_Edited # number of WT reads
readsEditedNotFrameshifted <- summ$Reads_Edited - summ$Reads_Frameshifted # number of edited but not frameshifted reads (indel multiples of 3)
readsNotKO <- readsWt + readsEditedNotFrameshifted # number of reads that do not lead to a knockout, i.e. are either Wt or have an indel multiple of 3
summ$notko <- (readsNotKO / summ$Reads_Filtered) # proportion of reads that do not lead to a knockout

summ <- cbind(summ, edit, frameshift, editbutnotfshift)

# check that SNPs are not counted
  # i.e. that number of reads edited = number of reads with insertions + number of reads with deletions

notaddition <- summ[summ$Reads_In + summ$Reads_Del != summ$Reads_Edited , ] # there are some samples where reads edited is NOT reads with in + reads with del

sum( notaddition$Reads_Edited - (notaddition$Reads_In + notaddition$Reads_Del) > 0 ) 
  # >> ALL of them are cases where there are LESS edited reads than reads with del + reads with in

# if SNPs were counted, there would be MORE edited reads than reads with del + reads with in
# here are cases where reads have both in and del
# eg. read X has a 5bp in & a 5bp del >> it is counted once in reads del & once in reads in, but only once in reads edited
# summary: correct; should leave it as it is

# ! that also means I should think carefully about it if tempting anywhere to sum Reads_Del and Reads_In, as that would double count reads

# rank the crRNAs per gene ------------------------------------------------------------

# when plotting number of guides as X axis, will assume user picks them going down the list from IDT
# I do not know exactly how IDT rank them
# but here can just follow alphabetical order to reproduce this
  # see biallelicprobability_rankedON.R for same section but classifies crRNAs based on ON scores

crRanked <- data.frame(matrix(data=NA, ncol=length(unique(summ$gene)), nrow=4)) # all genes have 3 guides, except slc24a5 which has 4
colnames(crRanked) <- unique(summ$gene)

for (g in colnames(crRanked)) {
  rows <- which(sub("\\_.*", "", summ$ID) == g) # which rows are about that gene g
  
  cr_names <- sort(unique(summ$crRNA_name [rows])) # picks the unique crRNA names of these
  # these are most likely classified correctly already (i.e. sort() probably unnecessary), but better be safe
  
  crRanked[1:length(cr_names), g] <- cr_names
}

# fish lookup table ------------------------------------------------------

# add a fish column
  # for most: within each gene, fish are the same
  # eg. tbx5a_AA fish1 is same as tbx5a_AB fish1 is same as tbx5a_AD fish1
    # ! slc24a5 is not exactly like this:
      # AA & AB & AD reads are from fish injected with AA + AB + AD
      # AC reads are from fish injected with AA + AB + AC
      # AG reads are from fish injected with AA + AB + AD + AG
summ$fish <- as.numeric(sapply(summ$ID, function(x) strsplit(x, '_')[[1]][3])) # using a regex and gsub would be cleaner but does the job

fislookup <- as.data.frame(matrix(nrow=4, ncol=length(unique(summ$gene))))
colnames(fislookup) <- unique(summ$gene)

for (g in colnames(fislookup)) {
  
  fisl <- list() # will store in there 3 or 4 vectors for each gene: fish numbers for crRNA AA; fish numbers for crRNA AB; etc.
  commonfis <- c()
  
  # annoyingly, it takes NA in the loop so need to clean up first
  crs2loop <- crRanked[g][!is.na(crRanked[g])]
  
  for (cr in 1:length(crs2loop)) {
    fisl[[cr]] <- summ$fish[which(summ$crRNA_name==crRanked[cr,g])]
    commonfis <- Reduce(intersect, fisl)
  }
  
  fislookup[1:length(commonfis), g] <- commonfis
}

# calculate Pbiallelic for all ----------------------------------------------------

pkosAll <- list() # will store dataframes, one per gene
# each dataframe will be: cols = fish; rows = Pko one guide/two guides/...

for (g in unique(summ$gene)) {
  
  pkos <- as.data.frame(matrix(nrow=nrow(crRanked[g]), ncol=sum(!is.na(fislookup[g])))) # preallocate the Pko dataframe of that gene
  
  # takes NA in the loop, so need to clean up first
  fis2loop <- fislookup[g][!is.na(fislookup[g])]
  
  for (f in 1:length(fis2loop)) { # ! loops thru fish index (so 1st fish, then 2nd, ...); not the actual fish ID (eg. 1, 3, 4)
    
    previousPnotko <- 1 # if user uses ONE guide: should not change the Pko
    
    cr2loop <- crRanked[g][!is.na(crRanked[g])] # takes NA in the loop, so need to clean up first
    for (rk in 1:length(cr2loop)) {
      # P of not making a KO using only that crRNA
        # i.e. proportion of reads not frameshifted by that crRNA
      notfm <- summ[which(summ$gene==g & 
                            summ$crRNA_name==crRanked[rk, g] & 
                            summ$fish==fislookup[f, g]) , 'notko']
      Pko <- 1 - (notfm * previousPnotko) # probability of having made a KO so far is 
      pkos[rk, f] <- Pko
      previousPnotko <- 1 - Pko # probability of not having made a KO so far
    }
  }
  
  # just before switching to next gene:
  colnames(pkos) <- sprintf('f%s', fis2loop) # so keep track of which fish is which
  pkosAll[[g]] <- pkos # add the dataframe to the list
}

# plot grid --------------------------------------------------------------------

# loop through list
# plot each
# add to list of ggplots
# ggarrange

pkosAll <- lapply(pkosAll, tbl_df)

mycols <- c('#b3bdc4', '#8d9ca6', '#697a87', '#4b5860')
mycols2 <- c('#4b5860', '#4b5860', '#4b5860', '#4b5860')

pkoPlots <- list() # will store the ggplots in there

for (pkog in 1:length(pkosAll)) {
    
  pkoPlots[[pkog]] <- 
      pkosAll[[pkog]] %>%
      drop_na() %>%
      mutate(guides=c(1:nrow(.))) %>%
      pivot_longer(-guides, names_to='fish', values_to='pko') %>%
    
      ggplot(., aes(x=guides, y=pko, col=fish)) +
        geom_line(size=0.5) +
        geom_point(size=0.6) +
        # ggtitle(names(pkosAll)[pkog]) +
        scale_colour_manual(values=mycols2) + # can revert back to mycols for shades of blue
        theme_minimal() +
        theme(
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
      
          # axis.title.x=element_text(size=9, margin = margin(t = 2, r = 0, b = 0, l = 0)),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
      
          axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          # axis.text.y=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text.y=element_blank(),
          
          plot.margin=margin(t = 4, r = 4, b = 4, l = 4), # controls spacing between plots in grid
      
          legend.position='none') +
      xlab('number of guides') + ylab('probability of knockout') +
      coord_cartesian(xlim=c(1,3), ylim=c(0,1)) + # if will not delete slc24a5/tyr; use xlim=c(1,4)
      scale_x_continuous(breaks=c(1, 2, 3, 4)) +
      scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1))
}

pkoPlot <- ggarrange(plotlist=pkoPlots, nrow=5, ncol=2) # better for A4
pkoPlot <- ggarrange(plotlist=pkoPlots, nrow=2, ncol=5) # better for ppt

# version for f0 paper
pkoPlots <- pkoPlots[c(1, 3, 2, 4, 5, 10, 8, 9, 6, 7)]
pkoPlots <- pkoPlots[-c(1, 2)] # 1 is slc24a5 / 2 is tyr

pkoPlot <- ggarrange(plotlist=pkoPlots, nrow=4, ncol=2) # f0 paper
pkoPlot_h <- ggarrange(plotlist=pkoPlots, nrow=1, ncol=8) # f0 paper; horizontal version

# plot slc24a5/tyr only -------------------------------------------------------

plotSingle <- function(gene) {
  
  pkosAll[[gene]] %>%
    drop_na() %>%
    mutate(guides=c(1:nrow(.))) %>%
    pivot_longer(-guides, names_to='fish', values_to='pko') %>%
    
    ggplot(., aes(x=guides, y=pko, col=fish)) +
    geom_hline(yintercept=0.9, linetype=2, colour='#ebebeb', size=0.5) +
    geom_line(size=0.5) +
    geom_point(size=0.6) +
    # ggtitle(gene) +
    scale_colour_manual(values=mycols2) +
    theme_minimal() +
    theme(
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      
      axis.title.x=element_text(size=9, margin = margin(t = 3, r = 0, b = 0, l = 0)),
      axis.title.y=element_text(size=9, margin = margin(t = 0, r = 3, b = 0, l = 0)),
      
      axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text.y=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      
      legend.position='none') +
    xlab('number of loci') + ylab('proportion of frameshift alleles') +
    coord_cartesian(xlim=c(1,4), ylim=c(0,1)) +
    scale_x_continuous(breaks=c(1, 2, 3, 4)) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1))
}

pkoplot_slc24a5 <- plotSingle('slc24a5')
pkoplot_tyr <- plotSingle('tyr')

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

ggsave (filename='f0_pKO_h.pdf', plot=pkoPlot_h, width=175, height=46, units='mm', useDingbats=FALSE) # f0 paper; horizontal version

ggsave (filename='f0_pKOslc24a5.pdf', plot=pkoplot_slc24a5, width=68, height=57, units='mm', useDingbats=FALSE)
ggsave (filename='f0_pKOtyr.pdf', plot=pkoplot_tyr, width=68, height=57, units='mm', useDingbats=FALSE)

# Pko after 3rd guide ---------------------------------------------
thirdguide <- as.data.frame(matrix(nrow=4, ncol=length(pkosAll)))
colnames(thirdguide) <- names(pkosAll)

for (gene in 1:length(pkosAll)) {
  thirdguide[1:ncol(pkosAll[[gene]]), gene] <- t(pkosAll[[gene]][3,])
}
# >> dataframe thirdguide is fish x gene.
  # data = P after 3rd guide

thirdguide <- tbl_df(thirdguide)
mean(unlist(thirdguide), na.rm=TRUE)
sd(unlist(thirdguide), na.rm=TRUE)
# can do more summary statistics if needed

apply(thirdguide, 2, function(x) mean(x, na.rm=TRUE))
length(which(unlist(thirdguide) > 0.8))
sum(!is.na(unlist(thirdguide)))


# probability increase when 4th guide added -------------------------------

# slc24a5
mean(as.numeric(pkosAll$slc24a5[3,]))
sd(as.numeric(pkosAll$slc24a5[3,]))

mean(as.numeric(pkosAll$slc24a5[4,]))
sd(as.numeric(pkosAll$slc24a5[4,]))

( mean(as.numeric(pkosAll$slc24a5[4,])) - mean(as.numeric(pkosAll$slc24a5[3,])) ) * 100

# tyr
mean(as.numeric(pkosAll$tyr[3,]))
sd(as.numeric(pkosAll$tyr[3,]))

mean(as.numeric(pkosAll$tyr[4,]))
sd(as.numeric(pkosAll$tyr[4,]))

( mean(as.numeric(pkosAll$tyr[4,])) - mean(as.numeric(pkosAll$tyr[3,])) ) * 100


# proportion of frameshift alleles for slc24a5/tyr/tbx16 -----------------------------

mean(unlist(thirdguide[,1:3]), na.rm=TRUE)
