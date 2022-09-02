#####################################################################################################################
#####################################################################################################################
#
#     Script to derive to requested tables from data call 
#       "Data call concerning fisheries in Marine Protected Areas in the Baltic 
#        and North Seas, Atlantic EU Western Waters and EU Outermost Regions"
#
#     By Karin van der Reijden (kjova@aqua.dtu.dk)
#     August 2022
#
#     This script is a modification of the original datacall script and allows for 
#     annual analyses that are merged afterwards in the hope to save computational effort and time.
#
#####################################################################################################################
#####################################################################################################################
rm(list = ls())

#-----------------------------------------------
# Information on this script
#-----------------------------------------------
### Before running this code the following needs to be done:
##   - Cleaned tacsatEflalo files for all the years should be in the RdataPath. Each file 
##     should be saved as tacsatEflalo_"Year".rds (eg. tacsatEflalo_2019.rds)
##     These files  should be in the standard tacsatEflalo format and contain KG and EURO for all
##     species. This can be obtained using the countries own procedures or using 
##     the script here: https://github.com/ices-eg/ICES-VMS-and-Logbook-Data-Call/blob/main/2_eflalo_tacsat_analysis.R (Note to remove line 260 to keep species-specific catch information). 
##   - In addition, the cleaned tacsatEflalo should contain these columns:
##     "VE_REF", "SI_LATI", "SI_LONG", "SI_DATE", "SI_TIME", "SI_SP", "SI_HE", "SI_STATE", "SI_DATIM", "INTV", "FT_REF"
##   - Cleaned eflalo files should be in the RdataPath, saved as cleanEflalo_"Year".rds (eg. cleanEflalo_2019.rds)
##   - The polygon shapefiles should be downloaded and extracted into the folder directed to by "polyPath"

################################################
# General settings
################################################
#-----------------------------------------------
# Load libraries
#-----------------------------------------------
library(sf)
library(data.table)
library(vmstools)

#-----------------------------------------------
# Set data paths
#-----------------------------------------------
sysPath                       <- "Please provide here the link to your VMS data folder"
RdataPath                     <- paste0(sysPath, "VMS/")
polyPath                      <- paste0(sysPath, "Data/")
intPath                       <- paste0(sysPath, "Int/")
resPath                       <- paste0(sysPath, "Results/")

# If directories do not exist, create them
dir.create(sysPath,    showWarnings = FALSE)
dir.create(polyPath,   showWarnings = FALSE)
dir.create(intPath,    showWarnings = FALSE)
dir.create(resPath,    showWarnings = FALSE)
dir.create(RdataPath,  showWarnings = FALSE)

#-----------------------------------------------
# Load MPA polygon data and merge to 1 file
#-----------------------------------------------
MPAs1                         <- st_read(paste0(polyPath, "MPA_polygons1_4326.shp"))
MPAs2                         <- st_read(paste0(polyPath, "buffer_polygons_4326.shp"))
MPAs                          <- rbind(MPAs1, MPAs2)
MPAs                          <- st_transform(MPAs, 4326) # should be superfluous as the (new) MPAs shapefiles already have this projection
rm(MPAs1, MPAs2)

#-----------------------------------------------
# Set some variables
#-----------------------------------------------
Country                       <- "DK"
yearsToSubmit                 <- c(2012:2021)

#-----------------------------------------------



################################################
# Processing of tacset and eflalo files
################################################
#-----------------------------------------------
# Read in the data, per year, and process
#-----------------------------------------------
for(iYear in yearsToSubmit) {
  print(paste(Sys.time(), iYear))
  
  tacsatEflalo                <- readRDS(paste0(RdataPath,"tacsatEflalo_",iYear,".rds"))
  eflalo                      <- readRDS(paste0(RdataPath,"cleanEflalo_",iYear,".rds"))
  
  # ## If data is saved as .RData, use this instead
  # load(file=paste0(RdataPath,"tacsatEflalo_",Year,".RData"))
  # load(file=paste0(RdataPath,"cleanEflalo_",Year,".RData"))
  
  #-----------------------------------------------
  # Prepare tacsatEflalo and eflalo
  #-----------------------------------------------
  # Add relevant information to tacsatEflalo
  tacsatEflalo$LE_RECT        <- ICESrectangle(tacsatEflalo)
  tacsatEflalo$LE_GEAR        <- eflalo$LE_GEAR[match(tacsatEflalo$FT_REF, eflalo$FT_REF)]
  tacsatEflalo$LE_QUART       <- as.integer(quarter(tacsatEflalo$SI_DATIM))
  tacsatEflalo$VE_KW          <- eflalo$VE_KW[match(tacsatEflalo$FT_REF, eflalo$FT_REF)]
  tacsatEflalo$LE_MET         <- eflalo$LE_MET[match(tacsatEflalo$FT_REF, eflalo$FT_REF)] # This should be the level 6 metier definitions. Please adjust the script if needed.
  eflalo$LE_QUART             <- quarter(as.Date(eflalo$LE_CDAT, format="%d/%m/%Y"))
  eflalo$LE_YEAR              <- iYear
  tacsatEflalo$LE_YEAR        <- iYear
  
  # Make data.frames into data.tables for faster processing
  te                          <- data.table(tacsatEflalo); gc()
  e                           <- data.table(eflalo)
  
  #-----------------------------------------------
  # Calculate total landings (kg / EURO)
  #-----------------------------------------------
  # Calculate total kg and EURO if not already present in data as LE_KG_TOT and LE_EURO_TOT
  if(!"LE_KG_TOT" %in% colnames(te))
    te[,LE_KG_TOT:=rowSums(.SD, na.rm = TRUE), .SDcols = grep("KG", names(te))]
  if(!"LE_EURO_TOT" %in% colnames(te))
    te[,LE_EURO_TOT:=rowSums(.SD, na.rm = TRUE), .SDcols = grep("EURO", names(te))]
  if(!"LE_KG_TOT" %in% colnames(e))
    e[,LE_KG_TOT:=rowSums(.SD, na.rm = TRUE), .SDcols = grep("KG", names(e))]
  if(!"LE_EURO_TOT" %in% colnames(e))
    e[,LE_EURO_TOT:=rowSums(.SD, na.rm = TRUE), .SDcols = grep("EURO", names(e))]
  
  # Calculate total fishing effort in days and kwdays 
  te[, LE_EFF_DAY := INTV/1440]
  te[, LE_EFF_KWD := (LE_EFF_DAY) * VE_KW]
  
  #-----------------------------------------------
  # Save te and e for later usage
  #-----------------------------------------------
  saveRDS(te, paste0(intPath, "TACEFL_", iYear, ".rds"))
  saveRDS(e, paste0(intPath, "EFLA_", iYear, ".rds"))
} # end iYear-loop

#-----------------------------------------------



################################################
# Create Table 1 (per year and merge results)
################################################
#-----------------------------------------------
# Create format, obtain data and process
#-----------------------------------------------
Table1                       <- data.table(LE_YEAR = integer(),
                                           Csquare = character(),
                                           LE_MET  = character(),
                                           Land_kg = numeric(),
                                           Eff_day = numeric(),
                                           Eff_kwh = numeric())

for(iYear in yearsToSubmit) {
  dat                        <- readRDS(paste0(intPath, "TACEFL_", iYear, ".rds"))
  # Set NA's to 0 
  dat[is.na(dat)]            <-  0
  
  # Assign Csquare based on lon and lat
  dat$Csquare                <- CSquare(lon = dat$SI_LONG,
                                        lat = dat$SI_LATI,
                                        degrees = 0.05)
  
  #-----------------------------------------------
  # Fill Table 1 with data
  #-----------------------------------------------
  # Aggregate data and write to results folder
  tab1                       <- dat[, .(Land_kg = sum(LE_KG_TOT),
                                        Eff_day = sum(LE_EFF_DAY, na.rm=T),
                                        Eff_kwh = sum(LE_EFF_KWD, na.rm=T)),
                                    by=.(LE_YEAR, Csquare, LE_MET)]
  Table1                     <- rbind(Table1, tab1)
} # end iYear-loop

# Save as csv-file
fwrite(Table1, paste0(resPath, "Table_1_", Country, ".csv"))
rm(dat, tab1, Table1)


#-----------------------------------------------



################################################
# Create Tables 2, 3a and 3b
################################################
#-----------------------------------------------
# Create yearly MPA-aggregates
#-----------------------------------------------
for(iYear in yearsToSubmit) {
  print(paste(Sys.time(), iYear))
  
  te                         <- readRDS(file=paste0(intPath, "TACEFL_", iYear, ".rds"))
  #-----------------------------------------------
  # Convert to spatial points (sf)
  #-----------------------------------------------
  p_te                       <- te %>%
    sf::st_as_sf(coords = c("SI_LONG", "SI_LATI")) %>%
    sf::st_set_crs(4326)
  
  #-----------------------------------------------
  # Add MPAs to file and save points inside polygons
  #-----------------------------------------------
  sub                        <- st_join(p_te, MPAs, join = st_intersects)
  sub                        <- as.data.table(st_drop_geometry(sub[!is.na(sub$SITECODE),]))
  saveRDS(sub, paste0(intPath, "MPA_", iYear, ".rds"))
} # end iYear-loop
rm(te, p_te, sub)

#-----------------------------------------------
# Create table 2 (per year and merge results)
#-----------------------------------------------
Table2                       <- data.table(LE_YEAR  = integer(),
                                           LE_QUART = integer(),
                                           SITECODE = character(),
                                           LE_MET   = character(),
                                           Land_kg  = numeric(),
                                           Land_eu  = numeric(),
                                           Eff_day  = numeric(),
                                           Eff_kwh  = numeric())
for(iYear in yearsToSubmit){
  dat                        <- readRDS(paste0(intPath, "MPA_", iYear, ".rds"))
  # Aggregate data and write to results folder
  tab2                       <- dat[, .(Land_kg = sum(LE_KG_TOT),
                                        Land_eu = sum(LE_EURO_TOT),
                                        Eff_day = sum(LE_EFF_DAY, na.rm=T),
                                        Eff_kwh = sum(LE_EFF_KWD, na.rm=T)),
                                    by=.(LE_YEAR, LE_QUART, SITECODE, LE_MET)]
  Table2                     <- rbind(Table2, tab2)
} # end iYear-loop

# Save as csv-file
fwrite(Table2, paste0(resPath, "Table_2_", Country, ".csv"))


#-----------------------------------------------
# Obtain data for entire time period
#-----------------------------------------------
m_list                       <- list.files(path = intPath, pattern= "MPA", full.names = T)
dat                          <- rbindlist(lapply(m_list, readRDS), fill=T)
# Set NA's to 0 
dat[is.na(dat)]              <-  0

#-----------------------------------------------
# Create table 3a
#-----------------------------------------------
table3a                      <- data.table(SITECODE  = character(),
                                           Spec_1_kg = character(),
                                           Spec_2_kg = character(),
                                           Spec_3_kg = character(),
                                           Spec_4_kg = character(),
                                           Spec_5_kg = character(),
                                           Rest_kg   = numeric())

for(iMPA in unique(dat$SITECODE)){
  subdat                    <- dat[SITECODE == iMPA,,]
  idx                       <- grep("LE_KG", names(subdat))
  a                         <- sort(colSums(subdat[, ..idx,]), decreasing=TRUE)
  idx2                      <- grep("TOT", names(a))
  b                         <- a[-idx2]
  d                         <- b[!b==0]
  dt                        <- data.table(SITECODE  = iMPA,
                                          Spec_1_kg = substr(names(d)[1], 7, 9),
                                          Spec_2_kg = substr(names(d)[2], 7, 9),
                                          Spec_3_kg = substr(names(d)[3], 7, 9),
                                          Spec_4_kg = substr(names(d)[4], 7, 9),
                                          Spec_5_kg = substr(names(d)[5], 7, 9),
                                          Rest_kg   = round((sum(d[6:length(d)]) / a[idx2]) *100, digits=1))
  dt$Rest_kg [is.na(dt$Rest_kg)==TRUE] <- 0
  table3a                   <- rbind(table3a, dt)
} # end iMPA-loop

# Save as csv-file
fwrite(table3a, paste0(resPath, "Table_3a_", Country, ".csv"))

#-----------------------------------------------
# Create table 3b
#-----------------------------------------------
table3b                      <- data.table(SITECODE    = character(),
                                           Spec_1_euro = character(),
                                           Spec_2_euro = character(),
                                           Spec_3_euro = character(),
                                           Spec_4_euro = character(),
                                           Spec_5_euro = character(),
                                           Rest_euro   = numeric())

for(iMPA in unique(dat$SITECODE)){
  subdat                    <- dat[SITECODE == iMPA,,]
  idx                       <- grep("LE_EURO", names(subdat))
  a                         <- sort(colSums(subdat[, ..idx,]), decreasing=TRUE)
  idx2                      <- grep("TOT", names(a))
  b                         <- a[-idx2]
  d                         <- b[!b==0]
  dt                        <- data.table(SITECODE     = iMPA,
                                          Spec_1_euro = substr(names(d)[1], 9, 11),
                                          Spec_2_euro = substr(names(d)[2], 9, 11),
                                          Spec_3_euro = substr(names(d)[3], 9, 11),
                                          Spec_4_euro = substr(names(d)[4], 9, 11),
                                          Spec_5_euro = substr(names(d)[5], 9, 11),
                                          Rest_euro   = round((sum(d[6:length(d)]) / a[idx2]) *100, digits=1))
  dt$Rest_euro [is.na(dt$Rest_euro)==TRUE] <- 0
  table3b                   <- rbind(table3b, dt)
} # end iMPA-loop

# Save as csv-file
fwrite(table3b, paste0(resPath, "Table_3b_", Country, ".csv"))

rm(m_list, dat, tab2, Table2, table3a, table3b, subdat, idx, a, iMPA, b, d, dt)
#-----------------------------------------------



################################################
# Create Tables 4a and 4b (per year and merge results)
################################################
#-----------------------------------------------
# Create table 4a
#-----------------------------------------------
Table4a                      <- data.table(LE_YEAR  = integer(),
                                           LE_QUART = integer(),
                                           LE_RECT  = character(),
                                           LE_MET   = character(),
                                           Land_kg  = numeric())
# Obtain data per year
for(iYear in yearsToSubmit){
  dat                        <- readRDS(paste0(intPath, "TACEFL_", iYear, ".rds"))
  # Set NA's to 0 
  dat[is.na(dat)]            <-  0
  # Aggregate data and merge all the years
  tab4a                      <- dat[, .(Land_kg = sum(LE_KG_TOT)),
                                    by=.(LE_YEAR, LE_QUART, LE_RECT, LE_MET)]
  Table4a                    <- rbind(Table4a, tab4a)
} # end iYear-loop

# Save as csv-file
fwrite(Table4a, paste0(resPath, "Table_4a_", Country, ".csv"))
rm(dat, tab4a, Table4a)

#-----------------------------------------------
# Create table 4b
#-----------------------------------------------
Table4b                      <- data.table(LE_YEAR  = integer(),
                                           LE_QUART = integer(),
                                           LE_RECT  = character(),
                                           LE_MET   = character(),
                                           Land_kg  = numeric())


# Obtain data per year
for(iYear in yearsToSubmit){
  dat                        <- readRDS(paste0(intPath, "EFLA_", iYear, ".rds"))
  # Change NA to 0 (only numeric columns, the other one doesnt work when there is NA's in timedate columns)
  for(j in seq_along(dat)){
    set(dat, i = which(is.na(dat[[j]]) & is.numeric(dat)), j = j, value = 0)
  } 
  
  # Aggregate data and merge all the years
  tab4b                      <- dat[, .(Land_kg = sum(LE_KG_TOT)),
                                    by=.(LE_YEAR, LE_QUART, LE_RECT, LE_MET)]
  Table4b                    <- rbind(Table4b, tab4b)
} # end iYear-loop

# Save as csv-file
fwrite(Table4b, paste0(resPath, "Table_4b_", Country, ".csv"))
rm(dat, Table4b, tab4b)
#-----------------------------------------------






################################################
# New code to create tables 3a and 3b in a less computational-demanding way
#
#  Note: if you have successfully created tables 3a and 3b using the code above, there is no need to also run this part of the code.
################################################
# #-----------------------------------------------
# # Create annual tables with total catch /revenue per species
# #-----------------------------------------------
# for(iYear in yearsToSubmit) {
#   MPAdat                      <- readRDS(paste0(intPath, "MPA_", iYear, ".rds"))
#   kg                          <- c(grep("LE_KG_", colnames(MPAdat)), ncol(MPAdat))
#   kgs                         <- MPAdat[,..kg,]
#   kgs                         <- kgs[, lapply(.SD, sum, na.rm=TRUE), by="SITECODE", .SDcols=colnames(kgs)[1:ncol(kgs)-1]] 
#   kgs$Year                    <- iYear
#   saveRDS(kgs, file=paste0(intPath, "KGS_", iYear, ".rds"))
#   
#   eur                         <- c(grep("LE_EURO", colnames(MPAdat)), ncol(MPAdat))
#   euros                       <- MPAdat[,..eur]
#   euros                       <- euros[, lapply(.SD, sum, na.rm=TRUE), by="SITECODE", .SDcols=colnames(euros)[1:ncol(euros)-1]] 
#   euros$Year                  <- iYear
#   saveRDS(euros, file=paste0(intPath, "EUROS_", iYear, ".rds"))
# } # end iYear-loop
# 
# #-----------------------------------------------
# # Create Table 3a
# #-----------------------------------------------
# ## Read in data
# m_list                       <- list.files(path = intPath, pattern= "KGS_", full.names = T)
# dat                          <- rbindlist(lapply(m_list, readRDS), fill=T)
# ## Set NA's to 0 
# dat[is.na(dat)]              <-  0
# 
# ## Create base table
# table3a                      <- data.table(SITECODE  = character(),
#                                            Spec_1_kg = character(),
#                                            Spec_2_kg = character(),
#                                            Spec_3_kg = character(),
#                                            Spec_4_kg = character(),
#                                            Spec_5_kg = character(),
#                                            Rest_kg   = numeric())
# 
# for(iMPA in unique(dat$SITECODE)){
#   subdat                    <- dat[SITECODE == iMPA,,]
#   idx                       <- grep("LE_KG", names(subdat))
#   a                         <- sort(colSums(subdat[, ..idx,]), decreasing=TRUE)
#   idx2                      <- grep("TOT", names(a))
#   b                         <- a[-idx2]
#   d                         <- b[!b==0]
#   dt                        <- data.table(SITECODE  = iMPA,
#                                           Spec_1_kg = substr(names(d)[1], 7, 9),
#                                           Spec_2_kg = substr(names(d)[2], 7, 9),
#                                           Spec_3_kg = substr(names(d)[3], 7, 9),
#                                           Spec_4_kg = substr(names(d)[4], 7, 9),
#                                           Spec_5_kg = substr(names(d)[5], 7, 9),
#                                           Rest_kg   = round((sum(d[6:length(d)]) / a[idx2]) *100, digits=1))
#   dt$Rest_kg [is.na(dt$Rest_kg)==TRUE] <- 0
#   table3a                   <- rbind(table3a, dt)
# } # end iMPA-loop
# 
# # Save as csv-file
# fwrite(table3a, paste0(resPath, "Table_3a_", Country, ".csv"))
# 
# 
# #-----------------------------------------------
# # Create Table 3b
# #-----------------------------------------------
# ## Read in data
# m_list                       <- list.files(path = intPath, pattern= "EUROS_", full.names = T)
# dat                          <- rbindlist(lapply(m_list, readRDS), fill=T)
# ## Set NA's to 0 
# dat[is.na(dat)]              <-  0
# 
# ## Create base table
# table3b                      <- data.table(SITECODE    = character(),
#                                            Spec_1_euro = character(),
#                                            Spec_2_euro = character(),
#                                            Spec_3_euro = character(),
#                                            Spec_4_euro = character(),
#                                            Spec_5_euro = character(),
#                                            Rest_euro   = numeric())
# 
# for(iMPA in unique(dat$SITECODE)){
#   subdat                    <- dat[SITECODE == iMPA,,]
#   idx                       <- grep("LE_EURO", names(subdat))
#   a                         <- sort(colSums(subdat[, ..idx,]), decreasing=TRUE)
#   idx2                      <- grep("TOT", names(a))
#   b                         <- a[-idx2]
#   d                         <- b[!b==0]
#   dt                        <- data.table(SITECODE     = iMPA,
#                                           Spec_1_euro = substr(names(d)[1], 9, 11),
#                                           Spec_2_euro = substr(names(d)[2], 9, 11),
#                                           Spec_3_euro = substr(names(d)[3], 9, 11),
#                                           Spec_4_euro = substr(names(d)[4], 9, 11),
#                                           Spec_5_euro = substr(names(d)[5], 9, 11),
#                                           Rest_euro   = round((sum(d[6:length(d)]) / a[idx2]) *100, digits=1))
#   dt$Rest_euro [is.na(dt$Rest_euro)==TRUE] <- 0
#   table3b                   <- rbind(table3b, dt)
# } # end iMPA-loop
# 
# # Save as csv-file
# fwrite(table3b, paste0(resPath, "Table_3b_", Country, ".csv"))
