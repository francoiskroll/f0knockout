# 26/03/2020
# from trpa1b data
# plot: % of wt in KO pool vs statistical power to detect a difference

library(ggpubr)
library(esc)
library(pwr)
library(dplyr)
library(tidyr)

# import ------------------------------------------------------------------
file <- '~/.../dpxresults.csv'
dpxtot <- read.csv(file, header=TRUE)

stats <- dpxtot %>%
  group_by(grp) %>%
  summarise_at(vars(delta),
               list(
                 mean=mean,
                 sd=sd,
                 median=median,
                 n=length
               ))

dpxtot$grp <- factor(dpxtot$grp,
                     levels=c('scr', 'trpa1b'))

# will exclude some fish --------------------------------------------------

mycols <- c('#B3BDC4', '#f1876b')

densityraw <- ggplot(dpxtot, aes(x=delta, fill=grp, colour=grp)) + # ! last iteration will be element n + 1 (eg. 21)
  geom_density(alpha=0.2) +
  geom_vline(data=stats, aes(xintercept=mean, color=grp),
             linetype='dashed') +
  scale_fill_manual(values=mycols) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  coord_cartesian(xlim=c(-1000, 10000))
densityraw

rows2remove <- c()

# trpa1b fish that may not be full KO (in wt area based on density plot above)
dpxtot[which(dpxtot$grp=='trpa1b' & dpxtot$delta > 3650),] # rows 3, 9, 27 // fish 5, 17, 53
# >> checked on the lineplot, it is the 3 trpa1b fish that increase like wt
rows2remove <- c(rows2remove, which(dpxtot$grp=='trpa1b' & dpxtot$delta > 3650))

# scr fish that does not react (actually decreases in activity)
dpxtot[which(dpxtot$grp=='scr' & dpxtot$delta < 0),] # row 6 // fish 10
rows2remove <- c(rows2remove, which(dpxtot$grp=='scr' & dpxtot$delta < 0))

rows2remove <- sort(unique(rows2remove))
dpxtotid <- dpxtot[-rows2remove,] # dpxtot 'idealised'
# >> checked on the lineplot, looks like expected


# update stats ------------------------------------------------------------

statsid <- dpxtotid %>%
  group_by(grp) %>%
  summarise_at(vars(delta),
               list(
                 mean=mean,
                 sd=sd,
                 median=median,
                 n=length
               ))
# mean & median of trpa1b are more similar now

# how density look like now? ------------------------------------------------------

# >> I am assuming this situation is 100% full KO
mycols <- c('#B3BDC4', '#f1876b')
densityid <- ggplot(dpxtotid, aes(x=delta, fill=grp, colour=grp)) + # ! last iteration will be element n + 1 (eg. 21)
  geom_density(alpha=0.2) +
  geom_vline(data=statsid, aes(xintercept=mean, color=grp),
             linetype='dashed') +
  scale_fill_manual(values=mycols) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  coord_cartesian(xlim=c(-1000, 10000))
densityid

# >> very little overlap between the two distributions

# simulate fake data ----------------------------------------------------------------

# function to generate fake scr data
simulateScr <- function(n) {
  sim <- rnorm(n, statsid$mean[1], statsid$sd[1])
  return(sim)
}

# function to generate fake ko data
# same method as above takes too many negative values; which we can assume are not likely to be observed
# function below to pick only positive from the distribution
simulateKo <- function(n) {
  sim <- rnorm(n, statsid$mean[2], statsid$sd[2]) # pick a first set
  i <- 1
  while(length(which(sim < 0)) != 0) { # then replace negative values, until there is none left
    sim[which(sim < 0)] <- rnorm(length(which(sim < 0)), statsid$mean[2], statsid$sd[2])
    # cat('\n ran', i, 'times \n')
    i <- i+1
  }
  return(sim)
}


# check if working as expected --------------------------------------------

nfish <- 1000 # a lot to see if distributions look similar

dpxtotsim <- as.data.frame(cbind(
  sprintf('f%i', 1:(nfish*2)),
  c(rep('scr', nfish), rep('trpa1b', nfish)),
  c(simulateScr(nfish), simulateKo(nfish))
))

colnames(dpxtotsim) <- c('fish', 'grp', 'delta')

# conversions
dpxtotsim$grp <- as.factor(dpxtotsim$grp)
dpxtotsim$delta <- as.numeric(as.character(dpxtotsim$delta))

dpxtotsim <- tbl_df(dpxtotsim)

# simulated distributions -------------------------------------------------

statssim <- dpxtotsim %>%
  group_by(grp) %>%
  summarise_at(vars(delta),
               list(
                 mean=mean,
                 sd=sd,
                 median=median,
                 spln=length,
                 min=min,
                 max=max
               ))

mycols <- c('#B3BDC4', '#f1876b')
densitysim <- ggplot(dpxtotsim, aes(x=delta, fill=grp, colour=grp)) + # ! last iteration will be element n + 1 (eg. 21)
  geom_density(alpha=0.2) +
  geom_vline(data=statssim, aes(xintercept=mean, color=grp),
             linetype='dashed') +
  scale_fill_manual(values=mycols) +
  scale_colour_manual(values=mycols) +
  theme_minimal() +
  coord_cartesian(xlim=c(-1000, 10000))
densitysim
# >> the simulated data was taking too many negatives; corrected above
# >> now looking pretty similar to real data


# GRADUALLY INCREASING NUMBER OF KNOCKOUTS ------------------------------
# = gradually decreasing number of wild-types
# i.e. simulating injections getting better

# will assume 20 fish in each group here
# >> effect size will get bigger and bigger as more knockouts are added to the pool
# >> more & more likely to pick it up (i.e. power increasing)

# below Function to do a full simulation
  # i.e. gradually increasing number of knockouts in group

SimulateEF <- function(nf) { # nf = number of fish in each group
  # dataframe to store each iteration of simulated data
  
  SimReturn <- list() # will return a list so can return multiple things
  
  Sim <- as.data.frame(matrix(nrow=nf*2, ncol=3))
  colnames(Sim) <- c('fish', 'grp', 'delta')
  
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
    Sim$delta <- c(simScr, simKo)
    
    # calculate summary stats
    Simstats <- Sim %>%
      group_by(grp) %>%
      summarise_at(vars(delta),
                   list(
                     mean=mean,
                     sd=sd,
                     median=median,
                     spln=length # sample n
                   ))
    
    # plot density
    mycols <- c('#B3BDC4', '#f1876b')
    SimPlots[[i]] <- ggplot(Sim, aes(x=delta, fill=grp, colour=grp)) + # ! last iteration will be element n + 1 (eg. 21)
      geom_density(alpha=0.2) +
      geom_vline(data=Simstats, aes(xintercept=mean, color=grp),
                 linetype='dashed') +
      scale_fill_manual(values=mycols) +
      scale_colour_manual(values=mycols) +
      theme_minimal() +
      coord_cartesian(xlim=c(-1000, 10000))
    
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
  geom_hline(yintercept=20, linetype='dashed') +
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

# cannot use this if run via Jobs
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

ggsave (filename='density_original.pdf', plot=densityraw, width=120, height=88, units='mm', useDingbats=FALSE)
ggsave (filename='density_ideal.pdf', plot=densityid, width=120, height=88, units='mm', useDingbats=FALSE)
ggsave (filename='density_simulated.pdf', plot=densitysim, width=120, height=88, units='mm', useDingbats=FALSE)

ggsave (filename='simulation_grid.pdf', plot=SimsGrids[[1]], width=500, height=2000, units='mm', limitsize=FALSE, useDingbats=FALSE)
# to keep an example of a simulation

# simulation takes a while to run, so will export Global Environment here
save.image(file='trpa1b_ifsomewts_SampleSize_100fish10sim.RData')

ggsave(filename='trpa1b_nKOvsSample.pdf', plot=Simplot, width=500, height=300, units='mm', useDingbats=FALSE)