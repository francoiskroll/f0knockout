# 26/03/2020
# from csnk1db data/period
# plot: % of wt in KO pool vs sample size necessary to detect effect size

library(ggpubr)
library(esc)
library(pwr)
library(dplyr)
library(tidyr)

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

# keeping only DMSO-treated -------------------------------------------------------

biodo <- biod[-which(biod$grp=='scr_pf67' | biod$grp=='ko_pf67') ,] # o for dmso
rownames(biodo) <- NULL

# fix grp factor levels ----------------------------------------------------------------

biod$grp <- as.factor(biod$grp)
biod$grp <- factor(biod$grp, levels=c('scr_dmso', 'ko_dmso'))

biod$period <- as.numeric(as.character(biod$period))

# density plot ------------------------------------------------------------

mycols <- c('#f1876b', '#B3BDC4')

perioddensity <- ggplot(biodo, aes(x=period, fill=grp, colour=grp)) +
  geom_density(alpha=0.2) +
  # geom_vline(data=stats, aes(xintercept=mean, color=grp),
  #            linetype='dashed') +
  scale_fill_manual(values=mycols) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  coord_cartesian(xlim=c(23.5, 29))
perioddensity

# summary stats -----------------------------------------------------------

stats <- biodo %>%
  group_by(grp) %>%
  summarise_at(vars(period),
               list(
                 mean=mean,
                 sd=sd,
                 median=median,
                 min=min,
                 max=max,
                 spln=length
               ))

# >> I will not exclude any fish
# different than trpa1b which is not really a quantitative/population-wide trait
# as variation within the population here; cannot confidently say that a fish is probably in the SCR distribution

# >> I am assuming the current situation is the true/maximum effect size

# simulate fake data ----------------------------------------------------------------

simulateScr <- function(n) {
  sim <- rnorm(n, stats$mean[1], stats$sd[1])
  return(sim)
}

simulateKo <- function(n) {
  sim <- rnorm(n, stats$mean[2], stats$sd[2])
  return(sim)
}

# versions to draw only between min & max of each group -------------------

# # function to generate fake scr data
# simulateScr <- function(n) {
#   sim <- rnorm(n, stats$mean[1], stats$sd[1]) # pick a first set
#   i <- 1
#   while(length(which(sim < stats$min[1] | sim > stats$max[1] )) != 0) { # then replace negative values, until there is none left
#     sim[which(sim < stats$min[1] | sim > stats$max[1])] <- 
#       rnorm(length(which(sim < stats$min[1] | sim > stats$max[1])), stats$mean[1], stats$sd[1])
#     # cat('\n ran', i, 'times \n')
#     i <- i+1
#   }
#   return(sim)
# }
# 
# # function to generate fake ko data
# simulateKo <- function(n) {
#   sim <- rnorm(n, stats$mean[2], stats$sd[2]) # pick a first set
#   i <- 1
#   while(length(which(sim < stats$min[2] | sim > stats$max[2] )) != 0) { # then replace negative values, until there is none left
#     sim[which(sim < stats$min[2] | sim > stats$max[2])] <- 
#       rnorm(length(which(sim < stats$min[2] | sim > stats$max[2])), stats$mean[2], stats$sd[2])
#     # cat('\n ran', i, 'times \n')
#     i <- i+1
#   }
#   return(sim)
# }

# check if working as expected --------------------------------------------

nfish <- 1000 # a lot to see if distributions look similar

persim <- as.data.frame(cbind(
  sprintf('f%i', 1:(nfish*2)),
  c(rep('scr', nfish), rep('ko', nfish)),
  c(simulateScr(nfish), simulateKo(nfish))
))

colnames(persim) <- c('fish', 'grp', 'period')

# conversions
persim$grp <- as.factor(persim$grp)
persim$period <- as.numeric(as.character(persim$period))

persim <- tbl_df(persim)

# simulated distributions -------------------------------------------------

mycols <- c('#f1876b', '#B3BDC4')

perioddensitySim <- ggplot(persim, aes(x=period, fill=grp, colour=grp)) +
  geom_density(alpha=0.2) +
  # geom_vline(data=stats, aes(xintercept=mean, color=grp),
  #            linetype='dashed') +
  scale_fill_manual(values=mycols) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  coord_cartesian(xlim=c(23.5, 29))
perioddensitySim

# >> roughly OK; but would there be a better fit than a normal?

# GRADUALLY INCREASING NUMBER OF KNOCKOUTS ------------------------------
# = gradually decreasing number of wild-types
# i.e. simulating injections getting better

# >> effect size will get bigger and bigger as more knockouts are added to the pool
# >> more & more likely to pick it up (i.e. power increasing)

# below Function to do a full simulation
  # i.e. gradually increasing number of knockouts in group

SimulateEF <- function(nf) { # nf = number of fish in each group
  # dataframe to store each iteration of simulated data
  
  SimReturn <- list() # will return a list so can return multiple things
  
  Sim <- as.data.frame(matrix(nrow=nf*2, ncol=3))
  colnames(Sim) <- c('fish', 'grp', 'period')
  
  # list to store density plots
  SimPlots <- list() # will then be put into a grid
  
  # dataframe to store final results
  Simfinal <- as.data.frame(matrix(nrow=nf+1, ncol=7)) # n + 1 is number of iterations; as starts at 0 wt fish
  colnames(Simfinal) <- c('nko', 'nwt', 'mean_scr', 'sd_scr', 'mean_ko', 'sd_ko', 'minn')
  # number of KO in KO group /// number of WT in KO group ('escapers')
  
  i <- 1 # counter
  for (ko in 0:nf) { # ko for knockout; i.e. gradually increasing number of knockouts
    
    scr <- nf - ko # scr for scrambled; i.e. number of scr to put in the ko group
    
    cat('\n\n >> iteration ', i, ': ', ko, ' KO / ', scr, ' WT in F0 group \n', sep='')
    
    # simulate data
    simScr <- simulateScr(nf) # scr group does not change
    simKo <- c(simulateKo(ko), simulateScr(scr))
    
    # check simScr and simKo are both n fish
    if (length(simScr) != nf) stop('>>> Not the right number of simulated SCR fish')
    if (length(simScr) != nf) stop('>>> Not the right number of simulated KO fish')
    
    # place data in a dataframe
    Sim$fish <- sprintf('f%i', 1:(nf*2))
    Sim$grp <- c(rep('scr', nf), rep('ko', nf))
    Sim$period <- c(simScr, simKo)
    
    # calculate summary stats
    Simstats <- Sim %>%
      group_by(grp) %>%
      summarise_at(vars(period),
                   list(
                     mean=mean,
                     sd=sd,
                     median=median,
                     spln=length # sample n
                   ))
    
    # plot density
    mycols <- c('#B3BDC4', '#f1876b')
    SimPlots[[i]] <- ggplot(Sim, aes(x=period, fill=grp, colour=grp)) + # ! last iteration will be element n + 1 (eg. 21)
      geom_density(alpha=0.2) +
      geom_vline(data=Simstats, aes(xintercept=mean, color=grp),
                 linetype='dashed') +
      scale_fill_manual(values=mycols) +
      scale_colour_manual(values=mycols) +
      theme_minimal() +
      coord_cartesian(xlim=c(23.5, 29))
    
    # calculate effect size
    Simef <- esc_mean_sd(grp1m = Simstats$mean[1], grp1sd = Simstats$sd[1], grp1n = Simstats$spln[1], # mean, sd, n of group1
                         grp2m = Simstats$mean[2], grp2sd = Simstats$sd[2], grp2n = Simstats$spln[2], # mean, sd, n of group2
                         es.type = "d") # type is Cohen's d or Hedges' g
    
    # ef is > 5.65348; sample size calculation throws error
    # there might a theoretical maximum? Sample size is 2.0 at ef = 5.6; is that minimum theoretical sample size?
    # I will bring down any ef > 5.6 to 5.6; and report 2 as sample size
    # for confirmation; waiting for: https://stats.stackexchange.com/questions/457468/is-there-a-theoretical-maximum-to-cohens-d-effect-size
    
    if (Simef$es > 5.65348) {
      Simef$es <- 5.65348
      cat('\t \t >>> ES was > 5.6 \n')
      Simn <- list(n=2) # if EF is over 5.6; report 2 as sample size
      
    } else if (Simef$es < -5.65348) {
      Simef$es <- -5.65348
      cat('\t \t >>> ES was < -5.6 \n')
      Simn <- list(n=2) # if EF is below 5.6; report 2 as sample size
      
    } else { # if EF is between -5.6 and 5.6; calculate sample size
      # i.e. minimum sample size to detect that effect size at pval 0.05; power=0.8
      Simn <- pwr.t.test(n=NULL, d=Simef$es, sig.level=0.05, power=0.8)
    }
    
    # store final results
    Simfinal[i,] <- c(ko, scr, Simstats$mean[1], Simstats$sd[1], Simstats$mean[2], Simstats$sd[2], Simn$n)
    SimReturn$data <- Simfinal # first element of SimReturn
    
    # put density plots in a grid plot
    # ! will throw an error if wrong here
    Simgrid <- ggarrange(plotlist=SimPlots, nrow=ceiling(nf/3)+1, ncol=3) # that should leave enough slots
    SimReturn$grid <- Simgrid
    
    i <- i+1
  }
  
  # check last ko should be nfish
  # check last scr should be 0
  if (ko != nfish) stop('>>> Error')
  if (scr != 0) stop('>>> Error')
  
  # check that no NA left in final results
  if (sum(is.na(Simfinal)) != 0) stop('>>> Some NA left in final results')
  
  return(SimReturn)
  
}


# run multiple simulations ------------------------------------------------

# will store results in a list
# each element: one dataframe simFinal
nsim <- 10
nfish <- 100
SimsAll <- list()
SimsGrids <- list()

for (s in 1:nsim) {
  cat('\n\n\n SIMULATION ', s, '/', nsim, sep='')
  # run the simulation and store the result
  SimResult <- SimulateEF(nf=nfish)
  SimsAll[[s]] <- SimResult$data # store data of that simulation
  SimsGrids[[s]] <- SimResult$grid # store grid of that simulation
}

# extract data from simulations -------------------------------------------

Simsumm <- as.data.frame(matrix(nrow=nfish+1, ncol=2+nsim))# preallocate, will store summary of simulations
Simsumm[c(1,2)] <- SimsAll[[1]][, c(1,2)] # take metadata from first element of the list
Simsumm[3:(nsim+2)] <- sapply(SimsAll, function(x) x$minn) # extract pwr column from each simulation
colnames(Simsumm) <- c('nko', 'nwt', sprintf('sim%i', 1:nsim))

# in practice, sample size can only be a positive integer
# >> round to upper integer
Simsumm[3:(nsim+2)] <- ceiling(Simsumm[3:(nsim+2)])

# summary statistics ------------------------------------------------------
# need for plotting barplot mean +- sd

Simsumm <- tbl_df(Simsumm)

Simstats <- Simsumm %>%
  pivot_longer(-c(nko, nwt),
               names_to='sim',
               values_to='minn') %>%
  group_by(nko) %>%
  summarise_at(vars(minn),
               list(
                 mean=mean,
                 sd=sd,
                 median=median,
                 n=length
               ))

# plot --------------------------------------------------------------------

colfunc <- colorRampPalette(c('#697a87', '#f1876b'))
mycols <- colfunc(nfish+1)

Simplot <- ggplot(Simstats, aes(x=nko, y=mean, fill=as.factor(nko))) +
  geom_bar(stat='identity') +
  geom_linerange(aes(ymin=mean-sd, ymax=mean+sd)) +
  geom_hline(yintercept=16, linetype='dashed') +
  scale_fill_manual(values=mycols) +
  theme_minimal() +
  theme(
    panel.grid.minor=element_blank(),
    axis.title.x=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=7, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(size=5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size=5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    legend.position='none') +
  xlab('% biallelic KO in F0 pool') + ylab('minimum sample size') +
  coord_cartesian(ylim=c(0, 100)) +
  scale_x_continuous(breaks=c(0, 25, 50, 75, 100), labels=c('0', '25', '50', '75', '100'))
Simplot

# export ------------------------------------------------------------------

# cannot use this if running via Jobs
# setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 

# simulation takes a while to run, so will export Global Environment here
save.image(file='csnk1db_ifsomewts_SampleSize_100fish10sim.RData')

ggsave (filename='periodensity_dmso.pdf', plot=perioddensity, width=120, height=88, units='mm', useDingbats=FALSE)
ggsave (filename='periodensity_Sim.pdf', plot=perioddensitySim, width=120, height=88, units='mm', useDingbats=FALSE)
ggsave (filename='simulation_grid.pdf', plot=SimsGrids[[1]], width=500, height=2000, units='mm', limitsize=FALSE, useDingbats=FALSE)
# to keep an example of a simulation
ggsave(filename='csnk1db_nKOvsSample.pdf', plot=Simplot, width=500, height=300, units='mm', useDingbats=FALSE)
