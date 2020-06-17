# genotypeGenerator 

library(tidyxl)

# import ------------------------------------------------------------------
cat('Select plate map')
file <- file.choose()
# file <- '~/Dropbox/phd/010819_galnextrapheno/010899_04genotypeMapclean.xlsx'

datebox <- substr(file, 
                  which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+1, 
                  which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))]+9) # look at the /s to get the datebox

path <- substr (file, 1, which(strsplit(file, '')[[1]] == '/') [length(which(strsplit(file, '')[[1]] == '/'))])

platemap_name <- substr(file,
                        which(strsplit(file, '')[[1]] =='/') [length(which(strsplit(file, '')[[1]] =='/'))] +1,
                        nchar(file))

experiment_name <- substr(file,
                          which(strsplit(file, '')[[1]] =='/') [length(which(strsplit(file, '')[[1]] =='/'))-1] +1,
                          which(strsplit(file, '')[[1]] =='/') [length(which(strsplit(file, '')[[1]] =='/'))] -1)

cells <- xlsx_cells(file)
cells <- cells[-(118:120) ,] # should get rid of comments if any

frmts <- cells$style_format #formats
cnts <- cells$character #contents

# clean up ----------------------------------------------------------------
# keep only wells

# each X:Y is a row; can only be 96 wells this way
frmts <- frmts[c(15:26, 28:39, 41:52, 54:65, 67:78, 80:91, 93:104, 106:117)]
cnts <- cnts[c(15:26, 28:39, 41:52, 54:65, 67:78, 80:91, 93:104, 106:117)]

if (length(frmts) != 96) stop('Not 96 wells') # checkpoint
if (length(cnts) != 96) stop('Not 96 wells') # checkpoint

# prepare list ------------------------------------------------------------

# how many genotypes are there
if ('empty' %in% frmts) {
  genos <- sort(unique(frmts)[-which(unique(frmts)=='empty')])
  geno_names <- unique(cnts[!is.na(cnts)])[-which(unique(cnts[!is.na(cnts)])=='empty')]
} else {
  genos <- sort(unique(frmts))
  geno_names <- unique(cnts[!is.na(cnts)])
}

if (length(genos) != length(geno_names)) stop('Not the same number of genotypes and genotype names.') # checkpoint

# matches genotype number to genotype name
genonames_matched <- c()

for (G in 1:length(genos)) {
  geno <- which(frmts==genos[G])
  genoname <- unique(cnts[geno])
  if (length(genoname) != 1) stop('Check the plate map. Probably same colour used for multiple genotypes.') # checkpoint
  genonames_matched <- c(genonames_matched, genoname)
}

if (setequal(geno_names, genonames_matched) == FALSE) stop('Issue with extracting the genotype names.') # checkpoint

# what is the genotype with the most number of fish
geno_lgths <- c()
for (G in 1:length(genos)) {
  fishes <- sort(which(frmts==genos[G]))
  geno_lgths <- c(geno_lgths, length(fishes))
}
maxlgth <- max(geno_lgths)


# start writing README file -----------------------------------------------
readme_header <- rbind (paste('Experiment: ', experiment_name, sep=''), 
       paste('Date: ', substr(datebox,1,6), sep=''),
       paste('Box: ', substr(datebox, 8, 10), sep=''),
       paste('Plate map: ', platemap_name, sep=''),
       '',
       paste('Total n = ', sum(geno_lgths), ' in ', length(genos), ' genotypes', sep=''),
       paste('Number of empty wells = ', length(which(cnts=='empty')), sep=''),
       '')

write(readme_header, paste(path, datebox, '_README.txt', sep=''), sep='\t') # first write the README header


# build genotype lists ----------------------------------------------------
Genocols <- vector(mode='list', length=length(genos))

Genotype <- matrix(nrow=maxlgth, ncol=length(genos))
fishcheck <- c()
cat('\n')
for (G in 1:length(genos)) {
  fishes <- sort(which(frmts==genos[G]))
  message(cat(genos[G], '=', genonames_matched[G], '// n = ', length(fishes), '// : fish ', fishes)) # append to the README file
  fishcheck <- c(fishcheck, fishes)
  length(fishes) <- maxlgth
  Genocols[[G]] <- fishes
  Genotype[,G] <- fishes
}

fishcheck <- sort(c(fishcheck, which(frmts=='empty')))
if (!identical(fishcheck, 1:96)) stop('Something wrong: not all wells are taken or a well is in multiple genotypes') # checkpoint

cat('\n')
cat('Total n =', sum(geno_lgths), '\n')
# say the number of empty wells XXX
cat('Number of empty wells =', length(which(cnts=='empty')))

# Genotype is now a dataframe of `genotypes` columns. Each element = fish ID (integer)
# missing headers

# need two header rows: `Genotype1` dummy row and genotype names
# but R does not allow elements in a same column to have different data types
# solution: first create the genotype file, then append the fish IDs to it

header <- rbind(rep('Genotype1', length(genos)), genonames_matched)
rownames(header) <-  NULL

# write genotype file -----------------------------------------------------
# first write the header
write.table(header, paste(path, datebox, 'genotype.txt', sep=''), sep='\t', na='', row.names=FALSE, col.names=FALSE)
# then append the fish IDs
write.table(Genotype, paste(path, datebox, 'genotype.txt', sep=''), sep='\t', na='', row.names=FALSE, col.names=FALSE, append=TRUE)