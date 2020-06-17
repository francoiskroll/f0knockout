# 220719 per3 csnk1db
# period plot

library(Hmisc) # needed for stat_summary() in plot
library(dunn.test)
library(ggbeeswarm) # for better positions of jitter dots
library(tidyr)
library(dplyr)


# import ------------------------------------------------------------------
biod0_path <- '~/.../220719_periodBoth.csv'
biod0 <- read.csv(biod0_path, skip=20, na.strings=c("","NA")) # data we want only starts at row20


# tidy up -----------------------------------------------------------------

# remove any fish with Status IGNORED (i.e. manually excluded from period analysis)

if (sum(is.na(biod0$Status)) != 0) { # if Status column is not all NA
  biod0 <- biod0[- which(biod0$Status == 'IGNORED') ,]
}

biod <- biod0[, colSums(is.na(biod0)) < nrow(biod0)] # remove empty columns
 
colnames(biod) <- c('fish', 'grp', 'period',
                       'period_err', 'circPhaseToZero', 'circPhaseToWindow',
                       'circPhaseErr', 'absPhaseToZero', 'absPhaseToWindow',
                       'absPhaseErr', 'amplitude', 'amplitudeErr', 'gof', 'err')

rownames(biod) <- NULL

# fix grp factor levels ----------------------------------------------------------------

biod$grp <- as.factor(biod$grp)
biod$grp <- factor(biod$grp, levels=c('scr_dmso', 'ko_dmso', 'scr_pf67', 'ko_pf67'))

biod$period <- as.numeric(as.character(biod$period))

# plot --------------------------------------------------------------------
mycols <- c('#B3BDC4', '#f1876b', '#5a6974', '#a32e0f')

periodplot <- ggplot(biod, aes(x=grp, y=period, colour=grp)) +
  geom_jitter(size=2.0, position=position_jitter(0.2)) +
  stat_summary(fun.y=mean, geom='point', colour='#82898b', shape=3, size=2, stroke=1) +
  # stat_summary(fun.data=mean_sdl, fun.args=list(mult=1), geom="linerange", color="#82898b") + # this adds line +- sd
  scale_colour_manual(values=mycols) +
  coord_flip() +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y=element_text(size=7, margin = margin(t = 0, r = -6, b = 0, l = 0)),
    axis.title.x=element_text(size=9, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    # axis.title.y=element_blank(),
    legend.position='none') +
  ylab('period (hours)') + xlab('') + #could remove in theme() but it is to keep plot the same width with traces
  scale_y_continuous(breaks=c(24, 27, 30, 33, 36, 39)) + scale_x_discrete(labels=c(' ', ' ', ' ', ' '))
periodplot

# f0 paper ----------------------------------------------------------------

periodf0 <- ggplot(biod, aes(x=grp, y=period, colour=grp)) +
  # geom_jitter(size=2.0, position=position_jitter(0.2)) +
  geom_quasirandom(dodge.width = 0.9, groupOnX=TRUE, size=0.8, width=0.2) + # places dots better, more like they were inside a violin plot
  stat_summary(fun.y=mean, geom='point', colour='#020304', shape=3, size=2, stroke=1) +
  # stat_summary(fun.data=mean_sdl, fun.args=list(mult=1), geom="linerange", color="#82898b") + # this adds line +- sd
  scale_colour_manual(values=mycols) +
  coord_flip() +
  theme_minimal() +
  theme(
    # panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    # panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y=element_blank(),
    axis.title.x=element_text(size=9, margin = margin(t = 2, r = 0, b = 0, l = 0)),
    # axis.title.y=element_blank(),
    legend.position='none') +
  ylab('') + xlab('') + #could remove in theme() but it is to keep plot the same width with traces
  scale_y_continuous(breaks=c(24, 27, 30, 33, 36, 39)) + scale_x_discrete(labels=c(' ', ' ', ' ', ' '))
periodf0

# density plot ------------------------------------------------------------

biod <- tbl_df(biod)

summstats <- biod %>%
  group_by(grp) %>%
  summarise_at(vars(period),
               list(
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 median= ~ median(., na.rm=TRUE),
                 nfis= ~ length
               ))

perioddensity <- ggplot(biod, aes(x=period, fill=grp, colour=grp)) +
  geom_density(alpha=0.2) +
  # geom_vline(data=stats, aes(xintercept=mean, color=grp),
  #            linetype='dashed') +
  scale_fill_manual(values=mycols) +
  scale_colour_manual(values=mycols) +
  theme_minimal()
perioddensity

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ggsave (filename='f0_per3period.pdf', plot=periodf0, width=156, height=35, units='mm') # f0 paper

# summary statistics --------------------------------------------------------------

biod <- tbl_df(biod)

summstats <- biod %>% 
  group_by(grp) %>%
  summarise_at(vars(period),
               list(
                 mean=mean, 
                 sd=sd,
                 median=median
               ))

summstats$mean[2] -  summstats$mean[1] # difference between scr_dmso and ko_dmso = 1.40h = 84 minutes
summstats$mean[3] -  summstats$mean[1] # difference between scr_dmso and ko_dmso = 8.37h = 502 minutes

kruskal.test(period ~ grp, data=biod)
dunn.test(biod$period, biod$grp)

summary(aov(period ~ grp, data=biod))
pairwise.t.test(biod$period, biod$grp)
