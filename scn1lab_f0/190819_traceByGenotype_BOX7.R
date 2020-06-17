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
file <- '~/.../190813_07_deltapxsecsmooth.csv'

# guess genotype path/file from file name
datebox <- substr(file, which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+1, 
                  which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+9)
path <- substr (file, 1, which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))])

# genofile <- paste(datebox, 'genotype.txt', sep='')
genofile <- '080819_07genotype.txt'
genopath <- paste(path, genofile, sep='')

lbfile <- paste(datebox, '_lbsec.csv', sep='')
lbpath <- paste(path, lbfile, sep='')

dpx <- fread(file, header=FALSE)
geno <- read.delim(genopath, header=TRUE, skip = 1, na.strings=c("","NA"))
lb <- read.csv(lbpath, header=FALSE)

dpx <- as.data.frame(cbind(as.numeric(rownames(dpx)), dpx))

# ! relationship between xtime and lb a bit more complex than expected
# because dpx is filtered and some frames are removed

dpx <- dpx[lb[3,]:lb[7,],]

# do stuff ----------------------------------------------------------------
genofis <- vector(mode = "list", length = ncol(geno))

for (g in 1:ncol(geno)) {
  genofis[[g]] <- geno[,g][!is.na(c(geno[,g]))]
}

fis <- sort(unlist(genofis)) # pool & sort all fish ids from geno

if (length(fis)+1 != ncol(dpx)) stop('something about fish IDs') # check

colnames(dpx) <- c('sec', sprintf('f%s', fis))

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

# f0 paper ----------------------------------------------------------------

# f0 paper; colours so matches figs B/C colour coding
mycols <- c('#697a87', '#74210b')
ribboncols <- c('#ADB4BA', '#B29186') # manually setting ribbon colours so cover the grid
# ribbons colours: Eye dropper when ribbon colours are my cols + alpha = 0.5

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
ggsave (filename='f0_scn1labtrace_07_small.png', plot=dpxtrace, width=92, height=70, units='mm', bg='transparent', dpi=1000)
toc()