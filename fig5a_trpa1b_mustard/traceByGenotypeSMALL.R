# traceByGenotypeSMALL

# based on TraceByGenotype
# simplified here to be used in short experiments (no light boundaries for instance)

# ! Smoothing is very slightly different than MatLab smoothing, differs at the end points

# SMALL_v1: from delta_px_sq_sec & can control smoothing (easier than going back to MatLab to trial/error smoothing)

# packages ----------------------------------------------------------------
library(reshape)
library(ggplot2)
library (data.table)
library(tictoc)
library(DescTools)
library(gtools)
library(zoo)
library(schoolmath)
library(pracma)

sem <- function(x) sd(x)/sqrt(length(x))
sem_narm <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

### Smoothing parameters
bin <- 60 # in seconds

### Time window parameters

# PRE mustard
pre_start <- 60 * 10 # 10 minutes (into secs) into the video
pre_end <- 60 * 20 # 20 minutes (into secs) into the video

pre_start <- 60 * 17 # 10 minutes (into secs) into the video
pre_end <- 60 * 20 # 20 minutes (into secs) into the video

# POST mustard
post_start <- 60 * 21 + 20 # 10 minutes (into secs) into the video
post_end <- 60 * 31 + 20 # 20 minutes (into secs) into the video

post_start <- 60 * 21 + 22 # 10 minutes (into secs) into the video
post_end <- 60 * 24 + 22 # 20 minutes (into secs) into the video

# import ------------------------------------------------------------------
file <- '~/.../070220_07_deltapxsqsec.csv'

# guess genotype path/file from file name
datebox <- substr(file, which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+1, 
                  which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+9)
path <- substr (file, 1, which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))])

boxnum <- as.numeric(substr(datebox, 9, nchar(datebox)))

genofile <- paste(datebox, 'genotype.txt', sep='')
genopath <- paste(path, genofile, sep='')

# import
dpx <- read.csv(file, header=FALSE)
geno <- read.delim(genopath, header=TRUE, skip = 1, na.strings=c("","NA"))

# sort the genotype file ----------------------------------------------------------------
genofis <- vector(mode = "list", length = ncol(geno))

for (g in 1:ncol(geno)) {
  genofis[[g]] <- geno[,g][!is.na(c(geno[,g]))]
}

fis <- sort(unlist(genofis)) # pool & sort all fish ids from geno

colnames(dpx) <- c(sprintf('f%s', fis))

if (length(fis) != ncol(dpx)) stop('something about fish IDs') # check

# pick only time window we need --------------------------------------------------
dpxpre <- dpx[pre_start:pre_end,] # comment out to plot everything
dpxpost <- dpx[post_start:post_end,]

# add second column --------------------------------------------------------
dpxpre <- cbind(rownames(dpxpre), dpxpre)
colnames(dpxpre)[1] <- 'sec'

dpxpost <- cbind(rownames(dpxpost), dpxpost)
colnames(dpxpost)[1] <- 'sec'

# smoothing ---------------------------------------------------------------

# pre
tic('smoothing')
dpxsmoothpre <- apply(dpxpre[,2:ncol(dpxpre)], 2,
              FUN=function(x) movavg(x, bin, type='s'))
toc()

dpxsmoothpre <- as.data.frame(cbind(as.numeric(as.character(dpxpre[,1])), dpxsmoothpre)) # puts back sec column
colnames(dpxsmoothpre)[1] <- 'sec'

# post
tic('smoothing')
dpxsmoothpost <- apply(dpxpost[,2:ncol(dpxpost)], 2,
                      FUN=function(x) movavg(x, bin, type='s'))
toc()

dpxsmoothpost <- as.data.frame(cbind(as.numeric(as.character(dpxpost[,1])), dpxsmoothpost)) # puts back sec column
colnames(dpxsmoothpost)[1] <- 'sec'


# arrange by genotype -----------------------------------------------------

# pre
dpxsubspre <- vector(mode='list', length=ncol(geno))
colcheck <- vector()
for (g in 1:ncol(geno)) {
  cols <- match(sprintf('f%s', genofis[[g]]), colnames(dpxsmoothpre))
  dpxsubpre <- data.table(cbind(dpxsmoothpre$sec, dpxsmoothpre[,cols]))
  colnames (dpxsubpre)[1] <- 'sec'
  dpxsubmpre <- melt(dpxsubpre, id='sec') # data.table::melt on data.table object is 40x faster
  dpxsubmpre <- cbind(dpxsubmpre, rep(colnames(geno)[g], nrow(dpxsubmpre)))
  dpxsubmpre <- dpxsubmpre[,c(1,2,4,3)]
  colnames(dpxsubmpre) <- c('sec', 'fish', 'grp', 'px')
  dpxsubspre[[g]] <- dpxsubmpre
  colcheck <- c(colcheck, cols)
}

if (!(identical(sort(colcheck),2:ncol(dpxsmoothpre)))) stop('seems like not all columns were taken') # check

# post
dpxsubspost <- vector(mode='list', length=ncol(geno))
colcheck <- vector()
for (g in 1:ncol(geno)) {
  cols <- match(sprintf('f%s', genofis[[g]]), colnames(dpxsmoothpost))
  dpxsubpost <- data.table(cbind(dpxsmoothpost$sec, dpxsmoothpost[,cols]))
  colnames (dpxsubpost)[1] <- 'sec'
  dpxsubmpost <- melt(dpxsubpost, id='sec') # data.table::melt on data.table object is 40x faster
  dpxsubmpost <- cbind(dpxsubmpost, rep(colnames(geno)[g], nrow(dpxsubmpost)))
  dpxsubmpost <- dpxsubmpost[,c(1,2,4,3)]
  colnames(dpxsubmpost) <- c('sec', 'fish', 'grp', 'px')
  dpxsubspost[[g]] <- dpxsubmpost
  colcheck <- c(colcheck, cols)
}

if (!(identical(sort(colcheck),2:ncol(dpxsmoothpost)))) stop('seems like not all columns were taken') # check


# melt and summarise for plotting -----------------------------------------

# pre
dpxmpre <- data.table(rbindlist(dpxsubspre)) # 3x faster than rbind

dpxsumpre <- dpxmpre[,
               j=list(mean = mean(px), sd = sd(px), sem = sem(px), rep = length(px)),
               by=list(grp, sec)] # 2x faster than aggregate (still main bottleneck)

colnames(dpxsumpre) <- c('grp', 'sec', 'mean', 'sd', 'sem', 'n')

dpxsumpre <- as.data.frame(dpxsumpre)

# post
dpxmpost <- data.table(rbindlist(dpxsubspost)) # 3x faster than rbind

dpxsumpost <- dpxmpost[,
                     j=list(mean = mean(px), sd = sd(px), sem = sem(px), rep = length(px)),
                     by=list(grp, sec)] # 2x faster than aggregate (still main bottleneck)

colnames(dpxsumpost) <- c('grp', 'sec', 'mean', 'sd', 'sem', 'n')

dpxsumpost <- as.data.frame(dpxsumpost)

# plot for f0 paper --------------------------------------------------------------------
dpxsum$grp <-  as.factor(dpxsum$grp)
dpxsum$grp <- factor(dpxsum$grp, levels=colnames(geno))
mycols <- c('#697a87', '#f1876b')
ribboncols <- c('#D8DDE1', '#F8C2B6') # manually setting ribbon colours to covers the grid

xbreakspre <- c(dpxsumpre$sec[1]+60, dpxsumpre$sec[1]+120, dpxsumpre$sec[1]+180)

dpxtracepre <- ggplot (dpxsumpre, aes(x=sec, y=mean, group=grp, colour=grp)) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem, fill=grp), colour=NA) +
  geom_line(size=0.8) +
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=ribboncols) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank() ,
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 12, b = 0, l = 0)),
    axis.text.x=element_text(size=7, margin=margin(t = -8, r = 0, b = 0, l = 0)),
    axis.text.y = element_blank(),
    legend.position='none') +
  coord_cartesian(ylim=c(0, 60)) + # can remove this
  xlab('') + ylab('') +
  scale_x_continuous(breaks=xbreakspre, labels=c('1', '2', '3'))
dpxtracepre

# last second of pre is 1200
# first of post is 1282
# minute 4 will be during plate switch; then can label:
# minute 5 = 1320
# minute 6 = 1380
# minute 7 = 1440

xbreakspost <- c(dpxsumpre$sec[1]+(5*60), dpxsumpre$sec[1]+(6*60), dpxsumpre$sec[1]+(7*60))

dpxtracepost <- ggplot (dpxsumpost, aes(x=sec, y=mean, group=grp, colour=grp)) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem, fill=grp), colour=NA) +
  geom_line(size=0.8) +
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=ribboncols) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank() ,
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 12, b = 0, l = 0)),
    axis.text.x=element_text(size=7, margin=margin(t = -8, r = 0, b = 0, l = 0)),
    axis.text.y = element_blank(),
    legend.position='none') +
  coord_cartesian(ylim=c(0, 60)) + # can remove this
  xlab('') + ylab('') +
  scale_x_continuous(breaks=xbreakspost, labels=c('5', '6', '7'))
dpxtracepost

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='f0_mustard_pre.pdf', plot=dpxtracepre, width=55, height=90, units='mm')
ggsave (filename='f0_mustard_post.pdf', plot=dpxtracepost, width=55, height=90, units='mm')