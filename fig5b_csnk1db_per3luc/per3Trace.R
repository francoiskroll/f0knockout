# 220719
# per3 csnk1db

# packages ----------------------------------------------------------------
library(zoo)
library(dplyr)
library(tidyr)
library(ggplot2)

# sem formula -------------------------------------------------------------
sem <- function(x) sd(x)/sqrt(length(x))
sem_narm <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

# import cps ------------------------------------------------------------------
per3_path <- '~/.../ampbaselinedtr_normmean.csv'
per3_raw <- read.csv(per3_path, skip=9) # first 7 rows are comments
per3_raw[ncol(per3_raw)] <- NULL # imports an empty column at the end

# fish to exclude ---------------------------------------------------------
# see 220719_README.txt
# 8 fish not included in period analysis / 3 fish with very important cycling

cols2remove <- which( colnames(per3_raw) %in%
                        c('X6..G2..ko_dmso', 'X7..H2..ko_dmso', 'X11..N2..scr_dmso',
                          'X16..S2..ko_dmso', 'X27..AF2..ko_dmso', 'X60..BS2..ko_pf67'))
per3 <- per3_raw[, -cols2remove]


# fix column names --------------------------------------------------------

# extract genotype tag from the column name
genotags <- as.character(sapply(colnames(per3)[2:ncol(per3)], function(x) strsplit(x, '..', fixed=TRUE)[[1]][3]))

fisid <- as.character(sapply(colnames(per3)[2:ncol(per3)], function(x) strsplit(x, '..', fixed=TRUE)[[1]][1]))
fisid <- as.numeric(substr(fisid, 2, 3))
fisid <- sprintf('f%i', fisid)

tags <- paste(fisid, genotags, sep='_')

colnames(per3) <- c('zth', tags)

# smoothing ---------------------------------------------------------------
npoints <- 20 # control smoothing here; 1 = no smoothing
per3sm <- as.data.frame(apply(per3, 2, function(x){rollmean(x,npoints)}))

# melt --------------------------------------------------------------------

per3sm <- tbl_df(per3sm)

per3smm <- per3sm %>%
  pivot_longer(-zth,
               names_to='fish',
               values_to='cps')

per3smm$id <- as.character(sapply(per3smm$fish, function(x) strsplit(x, '_', fixed=TRUE)[[1]][1]))
  
per3smm$grp <- paste( as.character(sapply(per3smm$fish, function(x) strsplit(x, '_', fixed=TRUE)[[1]][2])) , 
                      as.character(sapply(per3smm$fish, function(x) strsplit(x, '_', fixed=TRUE)[[1]][3])), sep='_')

per3smm <- per3smm[,c(2, 4, 5, 1, 3)]


# compute mean +- sd of each timepoint for plotting -------------------------------------------------

per3smmSum <- per3smm %>%
  group_by(grp, zth) %>%
  summarise_at(vars(cps),
               list(
                 mean= ~ mean(.),
                 sd= ~ sd(.),
                 sem= ~ sem(.)
               ))

# create spread -----------------------------------------------------------
# could improve this section and make it into a nice function, but it is a start

per3smmSum$grp <- factor(per3smmSum$grp,
                         levels=c('scr_dmso', 'ko_dmso', 'scr_pf67', 'ko_pf67'))

space_btwngrps <- 0.4
space_ingrps <- 0.2

per3smmSum$meanspr <- rep(NA, nrow(per3smmSum))

for (g in 1:length(levels(per3smmSum$grp))) { # assumes grps go by pair
  group <- levels(per3smmSum$grp)[g]
  rowsg <- which(per3smmSum$grp == group)
  
  if (g == 1) {
    per3smmSum$meanspr[rowsg] <- per3smmSum$mean[rowsg] + 0
  }
  
  if (g == 2) {
    per3smmSum$meanspr[rowsg] <- per3smmSum$mean[rowsg] + space_ingrps
  }
  
  if (g == 3) {
    per3smmSum$meanspr[rowsg] <- per3smmSum$mean[rowsg] + space_btwngrps + space_ingrps
  }
  
  if (g == 4) {
    per3smmSum$meanspr[rowsg] <- per3smmSum$mean[rowsg] + space_btwngrps + space_ingrps + space_ingrps
  }
  
}

# check that all rows were filled
if (sum(is.na(per3smmSum$meanspr)) != 0) stop('Error during spreading')

minmax <- per3smmSum %>%
  group_by(grp) %>%
  summarise_at(vars(meanspr),
             list(
               min=min,
               max=max
             ))

# plot --------------------------------------------------------------------

mycols <- c('#B3BDC4', '#f1876b', '#5a6974', '#a32e0f')
ribboncols <- c('#D8DDE1', '#F8C2B6', '#ADB4BA', '#D2998E')

per3trace <- ggplot (per3smmSum, aes(x=zth, y=meanspr, group=grp, colour=grp)) +
  geom_ribbon(aes(ymin=meanspr-sem, ymax=meanspr+sem, fill=grp), colour=NA) + # colour here is outline of the ribbon
  geom_line(size=0.5) +
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=ribboncols) +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.x=element_text(size=9, margin = margin(t = 4, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=9, margin = margin(t = 0, r = 4, b = 0, l = 0)),
    legend.position='none') +
  xlab('circadian time (hours)') + ylab('bioluminescence') +
  scale_x_continuous(breaks=c(0, 24, 48, 72, 96, 120))
per3trace

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='f0_per3trace.pdf', plot=per3trace, width=165, height=80, units='mm', useDingbats=FALSE) # f0 paper
