# prepare per3 data for export to Biodare2

library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(openxlsx)
msgspace <- c('\n \t \t \t \t \t \t \t \t \t \t \t \t \t \t')

startdate <- '2020-02-10'

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# import ------------------------------------------------------------------
per3_path <- '~/.../220719_all.csv'
geno_path <- '~/.../220719genotype.txt'

geno <- read.table(geno_path, sep='\t', skip=1, header=TRUE)
per3_raw <- fread(per3_path, header=FALSE) # fread from data.table is fast to import

# tidy up ------------------------------------------------------------
per3_raw <- as.data.frame(per3_raw)
# per3_raw is essentially rows = wells / cols = timepoints
# but each timepoint has 3 cols: well number/CPS/time; i.e. every 3rd row is CPS

per3 <- per3_raw[, seq(2, ncol(per3_raw), 3)] # grab only the CPS
# >> per3 is now actually wells x timepoints

per3ts <- per3_raw[, seq(3, ncol(per3_raw), 3)] # same dimensions as per3 above, but contains the timestamp of each measure

per3 <- t(per3)
per3ts <- t(per3ts)
colnames(per3) <- sprintf('f%s', 1:ncol(per3))
colnames(per3ts) <- sprintf('f%s', 1:ncol(per3ts))

if (ncol(per3) != 96) 
  stop('Why do you not have 96 columns?')

per3 <- tbl_df(per3)
per3ts <- tbl_df(per3ts)

# fix dates ---------------------------------------------------------------

per3ts[] <- lapply(per3ts, function(x) paste(startdate, x, sep=' ')) # add the startdate to every timestamps (not actually accurate, will fix)
per3ts <- mutate_all(per3ts, ymd_hms) # all in good format, but dates inaccurate (only first day is ok)
# ts is for timestamps

mdsperwell <- c() # number of midnights per well, will use it to check that all the wells have same number of midnights
for (col in 1:ncol(per3ts)) {
  well <- per3ts[,col]
  tsdiff <- diff( do.call('c', well))
  mdnghts <- which(tsdiff < -1428 & tsdiff > -1432)
  cat(msgspace, '>> Well', col, '- found', length(mdnghts), 'midnights\n')
  mdsperwell <- c(mdsperwell, length(mdnghts))
  
  mdnghts <- c(mdnghts, nrow(well)) # add last row of that column to midnight indices
  mdnghts <- as.vector(mdnghts)
  mdind <- 1
  while(length(mdnghts)!=1) { # while there are still midnights (stop when 1 left = last timepoint)
    # below: add 1 day to all timepoints between midnight (value just after = first value of the new day) and the following midnight
    well[(mdnghts[1]+1) : mdnghts[1+1],] <- well[(mdnghts[1]+1) : mdnghts[1+1],] + mdind * 24*60*60
    per3ts[,col] <- well
    mdnghts <- mdnghts[-1] # remove that midnight
    mdind <- mdind + 1
  }
}

### some checks
if (length(mdsperwell) != 96) # checks it looped through 96 wells
  message(msgspace, '>> Error: did not find 96 well while looking for midnights')

if (length(unique(mdsperwell)) != 1) # checks there is the same number of midnights in each well
  message(msgspace, '>> Error: some well have different number of midnights')
###

# >>> timestamps should be correct now
# Summary: per3ts is same dimensions as per3, i.e. rows = timepoints / cols = wells (96)
# data in per3ts = timestamp of each measure // data in per3: counts per second

# add timestamps ----------------------------------------------------------

# TopCount takes around 9 minutes to scan the whole plate
# but need a single timestamp per scan, so no perfectly accurate solution
# will use the median time of each plate scan as timestamp

per3ts$zth <- 1:nrow(per3ts) # temporarily just row numbers
tstps <- apply(per3ts, 1, median) # median time of the plate scan; tsps for timestamps
zt0 <- ymd_hms(paste(startdate, '09:00:00', sep=' '))
zth <- as.numeric(difftime(tstps, zt0, unit='hours')) # Zeitgeber times, i.e. number of hours after Zeitgeber = first 9am

per3ts$zth <- zth
per3ts <- per3ts[, c(ncol(per3ts), 1:(ncol(per3ts)-1))] # and putting that column as first

per3$zth <- zth # adding that column to per3 cps
per3 <- per3[, c(ncol(per3), 1:(ncol(per3)-1))] # and putting that column as first

totalexp <- round(as.numeric(difftime(tstps[length(tstps)], tstps[1], unit='days')),2)

cat(msgspace, '>> Experiment was', totalexp, 'days\n')

# organise per3 columns for Biodare ------------------------------------------------------------

# putting genotype tags as column names
# >> allows Biodare to automatically assign genotypes

per3bd <- per3 # keep per3 as full 96 columns with fish ID as column names
# bd for Biodare

# check the genotype file
if ( sum(duplicated(sort(as.vector(unlist(geno)[!is.na(geno)])))) != 0 ) stop('Error: some fish assigned more than one genotype')

# extract empty wells for check below
allfis <- sort(as.vector(unlist(geno)[!is.na(geno)]))

emptywells <- which(! 1:96 %in% allfis) # which 1:96 is not in allfis

# assign genotypes
genocheck <- c()
genocols <- rep(NA, 96) # will store the future names of the columns

for (g in 1:ncol(geno)) {
  
  fis <- geno[g][!is.na(geno[g])] # take the fish IDs in genotype group g
  genocols[fis] <- rep(colnames(geno)[g], length(fis)) # assign the genotype for these column names
  
  genocheck <- c(genocheck, fis) # keep the columns we have assigned for check below
  
}


# CHECK
# genocheck (i.e. the columns we have assigned) should be the same as genotype file
if (!identical(sort(genocheck), sort(as.vector(unlist(geno)[!is.na(geno)])))) stop('Error during genotype assignments')
if (!identical( sort(which(is.na(genocols))) , sort(emptywells) )) stop('Error: something about empty wells')
if (ncol(per3bd)-1  != length(genocols)) stop('Error during genotype assignments')

# if OK: can assign genotypes as column names
colnames(per3bd) <- c('zth', genocols)

# and remove column which have NA as names
per3bd[, which(is.na(colnames(per3bd)))] <- NULL

# export ------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

expdate <- substr(per3_path, which(strsplit(per3_path, '')[[1]] == '/') [length(which(strsplit(per3_path, '')[[1]] == '/'))]+1, 
                  which(strsplit(per3_path, '')[[1]] == '/') [length(which(strsplit(per3_path, '')[[1]] == '/'))]+6)

write.csv(per3ts, paste(expdate, '_per3ts.csv', sep=''), row.names=FALSE)
write.csv(per3, paste(expdate, '_per3cps.csv', sep=''), row.names=FALSE)
write.csv(per3bd, paste(expdate, '_per3biodare.csv', sep=''), row.names=FALSE)