# TRACE BY GENOTYPE

# packages ----------------------------------------------------------------
library(reshape)
library(ggplot2)
library (data.table)
library(tictoc)
library(DescTools)
library(gtools)
library(zoo)

sem <- function(x) sd(x)/sqrt(length(x))
sem_narm <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

# import ------------------------------------------------------------------
file <- '~/.../140218_34_deltapxsecsmooth.csv'

# guess genotype path/file from file name
datebox <- substr(file, which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+1, 
                  which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+9)
path <- substr (file, 1, which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))])

genofile <- paste(datebox, 'genotype.txt', sep='')
genopath <- paste(path, genofile, sep='')

lbfile <- paste(datebox, '_lbsec.csv', sep='')
lbpath <- paste(path, lbfile, sep='')

dpx <- fread(file, header=FALSE)
geno <- read.delim(genopath, header=TRUE, skip = 1, na.strings=c("","NA"))
lb <- read.csv(lbpath, header=FALSE)

dpx <- as.data.frame(cbind(as.numeric(rownames(dpx)), dpx))

# x axis should be all secs between 4th in lb & 7th in lb
# not sure, does not seem to match the number of values in matlab

# ! relationship between xtime and lb a bit more complex than expected
# because dpx is filtered and some frames are removed
dpx <- dpx[lb[3,]:lb[7,],]
# dpx <- dpx[which(xtime[,1]==lb[3,])[1]:which(xtime[,1]==lb[7,])[1],]
# xtime <- xtime[which(xtime==lb[3,])[1]:which(xtime==lb[7,])[1],]

# If you need to try first, compute on a single dataframe
# dpx <- dpx[lb[5,]:lb[6,],]


# do stuff ----------------------------------------------------------------
genofis <- vector(mode = "list", length = ncol(geno))

for (g in 1:ncol(geno)) {
  genofis[[g]] <- geno[,g][!is.na(c(geno[,g]))]
}

fis <- sort(unlist(genofis)) # pool & sort all fish ids from geno

if (length(fis)+1 != ncol(dpx)) stop('something about fish IDs') # check

colnames(dpx) <- c('sec', sprintf('f%s', fis))

# I think deltapxsecsmooth should be already smoothed
# tic('smoothing')
# dpxsmooth_raw <- as.data.frame(rollapplyr(as.zoo(dpx[,2:ncol(dpx)]), width=bin_roll, 
#                                           FUN=function(x) mean(x, na.rm=TRUE), by=1, by.column=TRUE, partial=TRUE, fill=NA, align="right"))
# toc()

dpxsubs <- vector(mode='list', length=ncol(geno))
colcheck <- vector()
for (g in 1:ncol(geno)) {
  cols <- match(sprintf('f%s', genofis[[g]]), colnames(dpx))
  dpxsub <- data.table(cbind(dpx$sec, dpx[,cols]))
  colnames (dpxsub)[1] <- 'sec'
  dpxsubm <- melt(dpxsub, id='sec') # data.table::melt on data.table object is 40x faster
  dpxsubm <- cbind(dpxsubm, rep(colnames(geno)[g], nrow(dpxsubm)))
  dpxsubm <- dpxsubm[,c(1,2,4,3)]
  colnames(dpxsubm) <- c('sec', 'fish', 'grp', 'px')
  dpxsubs[[g]] <- dpxsubm
  colcheck <- c(colcheck, cols)
}

if (!(identical(sort(colcheck),2:ncol(dpx)))) stop('seems like not all columns were taken') # check


# dpx <- dpx [1:1000,]

dpxm <- data.table(rbindlist(dpxsubs)) # 3x faster than rbind

dpxsum <- dpxm[,
               j=list(mean = mean(px), sd = sd(px), sem = sem(px), rep = length(px)),
               by=list(grp, sec)] # 2x faster than aggregate (still main bottleneck)

colnames(dpxsum) <- c('grp', 'sec', 'mean', 'sd', 'sem', 'n')

# plot --------------------------------------------------------------------
# add light boundaries vertical lines
dpxsum$grp <-  as.factor(dpxsum$grp)
dpxsum$grp <- factor(dpxsum$grp, levels=c('WT_SCN1LAB', 'HET_SCN1LAB', 'HOM_SCN1LAB'))
lbs <- as.numeric(lb[,1])
mycols <- c('#697a87', '#b3bdc4', '#f1876b')

dpxtracelines <- ggplot (dpxsum, aes(x=sec, y=mean, group=grp, colour=grp)) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem, fill=grp), colour=NA, alpha=0.5) +
  geom_line(size=0.4) +
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=mycols) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank() ,
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=7, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.y = element_text(size=5, margin = margin(t = 0, r = -1, b = 0, l = 0)),
    legend.position='none') +
  xlab('Time (days/nights)') + ylab(expression(paste('Activity (total', Delta, 'px/sec)'))) +
  geom_vline(xintercept=lbs[4:(length(lbs)-1)], linetype = 2, size = 0.1, colour = 'black') +
  scale_x_continuous(breaks=lbs[4:(length(lbs)-2)], labels=NULL)


# version w/0 day/night lines ---------------------------------------------

dpxtrace <- ggplot (dpxsum, aes(x=sec, y=mean, group=grp, colour=grp)) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem, fill=grp), colour=NA, alpha=0.5) +
  geom_line(size=0.4) +
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=mycols) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank() ,
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=7, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=7, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.y = element_text(size=5, margin = margin(t = 0, r = -1, b = 0, l = 0)),
    legend.position='none') +
  xlab('Time (days/nights)') + ylab(expression(paste('Activity (total', Delta, 'px/sec)'))) +
  scale_x_continuous(breaks=lbs[4:(length(lbs)-2)], labels=NULL)


# f0 paper ----------------------------------------------------------------

ribboncols <- c('#ADB4BA', '#D8DDE1', '#F8C2B6') # manually setting ribbon colours so cover the grid

dpxtrace <- ggplot (dpxsum, aes(x=sec, y=mean, group=grp, colour=grp)) +
  geom_ribbon(aes(ymin=mean-sem, ymax=mean+sem, fill=grp), colour=NA) +
  geom_line(size=0.4) +
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=ribboncols) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank() ,
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 2, b = 0, l = 0)),
    axis.text.y = element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    legend.position='none') +
  xlab('') + ylab('') +
  # geom_vline(xintercept=lbs[4:(length(lbs)-1)], linetype = 2, size = 0.1, colour = 'black') +
  scale_x_continuous(breaks=lbs[4:(length(lbs)-2)], labels=NULL) +
  scale_y_continuous(labels=NULL)

# export ------------------------------------------------------------------

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

tic('png export')
ggsave (filename='stablescn1labtrace_small.png', plot=dpxtrace, width=92, height=70, units='mm', bg='transparent', dpi=1000)
toc()