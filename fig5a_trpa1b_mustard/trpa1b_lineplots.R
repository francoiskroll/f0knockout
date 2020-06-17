# trpa1b pre vs post line plots

# time windows ---------------------------------------------------------

start_pre <- 60 * 17 # 17min (into secs) into the video
end_pre <- 60 * 20 # 20min (into secs) into the video

start_post <- 60 * 21 + 22 # 21min22sec (into secs) into the video
end_post <- 60 * 24 + 22 # 24min22sec (into secs) into the video

# each window is 3 min

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
library(dplyr)
library(tidyr)
library(esc)
library(pwr)

sem <- function(x) sd(x)/sqrt(length(x))
sem_narm <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

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


# pick only time windows we need ------------------------------------------
# data type can cause issue, make sure everything is numeric
dpx <- as.data.frame(apply(dpx, 2,
                           FUN = function(x) as.numeric(as.character(x))))

dpx_pre <- dpx[start_pre:end_pre,]
dpx_post <- dpx[start_post:end_post,]

# 1/ sum total deltapx for each fish within time window -------------------
dpxtot_pre <- as.data.frame(apply(dpx_pre, 2, sum))

dpxtot_post <- as.data.frame(apply(dpx_post, 2, sum))

dpxtot <- cbind(dpxtot_pre, dpxtot_post)
colnames(dpxtot) <- c('pre', 'post')

# arrange by genotype -----------------------------------------------------
if (length(fis) != nrow(dpxtot)) stop('not the same number of rows as number of fish in genotype file') # check

genocolumn <- vector(mode='character', length=length(fis))

for (g in 1:ncol(geno)) { # for each genotype
  genorows <- match(sprintf('f%s', genofis[[g]]), rownames(dpxtot)) # find which rows are the fish of that genotype
  genocolumn[genorows] <-  colnames(geno)[g]
}

if (any(genocolumn=='')) stop('some fish were not assigned a genotype') # check

dpxtot <- as.data.frame(cbind(genocolumn, dpxtot))
dpxtot <- cbind(rownames(dpxtot), dpxtot)
colnames(dpxtot) <- c('fish', 'grp', 'pre', 'post')

dpxtot$grp <- as.character(dpxtot$grp)

dpxtot <- tbl_df(dpxtot)

# plot --------------------------------------------------------------------

dpxtot2 <- dpxtot
colnames(dpxtot2) <- c('fish', 'grp', 1, 2)

dpxtot2 <- tbl_df(dpxtot2)

mycols <- c('#697a87', '#f1876b')

slopePlot <- 
  
  dpxtot2 %>% 
  group_by(fish) %>%
  pivot_longer(-c(fish, grp), names_to='tp', values_to='px') %>%
  mutate(tp=as.numeric(tp)) %>%
  
  ggplot(., aes(x=tp, y=px, col=grp)) +
  geom_line(aes(group=fish), size=0.6) +
  geom_point(size=2.3) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    legend.position='none'
  ) +
  ylab(expression(paste('activity (total ', Delta, ' px)'))) +
  scale_x_continuous(breaks=c(1, 2), labels=c('pre', 'post'), expand=c(0.05,0.05)) # expand(x,x) increases padding on sides of the plot

slopePlot

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='trpa1b_slope.pdf', plot=slopePlot, width=65, height=88, units='mm', useDingbats=FALSE)

# STATISTICS --------------------------------------------------------------

# arrange data in delta px pre/post ---------------------------------------
dpxtot$delta <- dpxtot$post - dpxtot$pre

# t test by group ---------------------------------------------------------
t.test(delta ~ grp, data=dpxtot) # p-value = 4.894e-12
# mean scr = 6298.565 / mean trpa1b = 1078.182

stats <- dpxtot %>%
  group_by(grp) %>%
  summarise_at(vars(delta),
               list(
                 mean=mean,
                 sd=sd,
                 median=median,
                 n=length
               ))

# export results ------------------------------------------------------------------
write.csv(dpxtot, 'dpxresults.csv', row.names=FALSE)

# is data roughly normally distributed? -----------------------------------

ggplot(dpxtot, aes(x=delta, color=grp)) +
  geom_density() +
  geom_vline(data=stats, aes(xintercept=mean, color=grp),
             linetype='dashed')

# >> quite easy to see how trpa1b data is (probably) contaminated by a couple of wild-types
which( dpxtot$delta[which(dpxtot$grp=='trpa1b')] > 4000 ) # 2 of them