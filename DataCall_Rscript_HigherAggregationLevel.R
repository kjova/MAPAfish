#####################################################################################################################
#####################################################################################################################
#
#     Script to derive to more aggregated tables from data call 
#       "Data call concerning fisheries in Marine Protected Areas in the Baltic 
#        and North Seas, Atlantic EU Western Waters and EU Outermost Regions"
#
#     By Karin van der Reijden (kjova@aqua.dtu.dk)
#     September 2022
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
library(stringr)

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
# Processing of tacset and eflalo files           ## This is similar as in the 'normal' script
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
# Create Table 1 (adjusted to higher aggregation level)
################################################
#-----------------------------------------------
# Obtain data and process
#-----------------------------------------------
m_list                       <- list.files(path = intPath, pattern= "TACEFL", full.names = T)
dat                          <- rbindlist(lapply(m_list, readRDS), fill=T)
# Set NA's to 0 
dat[is.na(dat)]              <-  0

# Assign Csquare based on lon and lat
dat$Csquare                  <- CSquare(lon = dat$SI_LONG,
                                        lat = dat$SI_LATI,
                                        degrees = 0.05)

# Redefine metier (basically down to gear-level, e.g. TBB, OTB, PS, SDN, ...)
Metiers                      <- str_split_fixed(dat$LE_MET, pattern="_", n=6)
dat$LE_MET                   <- Metiers[,1] 

#-----------------------------------------------
# Create table 1
#-----------------------------------------------
# Aggregate data and write to results folder
table1                       <- dat[, .(Land_kg = sum(LE_KG_TOT),
                                        Eff_day = sum(LE_EFF_DAY, na.rm=T),
                                        Eff_kwh = sum(LE_EFF_KWD, na.rm=T)),
                                    by=.(LE_YEAR, Csquare, LE_MET)]

# Save as csv-file
fwrite(table1, paste0(resPath, "Table_1HAL_", Country, ".csv"))
rm(m_list, dat, table1)

#-----------------------------------------------





################################################
# Create Table 2 (adjusted to higher aggregation level)
################################################
#-----------------------------------------------
# Create yearly MPA-aggregates                    ## This is similar as in the 'normal' script
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
# Obtain data for entire time period
#-----------------------------------------------
m_list                       <- list.files(path = intPath, pattern= "MPA", full.names = T)
dat                          <- rbindlist(lapply(m_list, readRDS), fill=T)
# Set NA's to 0 
dat[is.na(dat)]              <-  0

# Redefine metier (basically down to gear-level, e.g. TBB, OTB, PS, SDN, ...)
Metiers                      <- str_split_fixed(dat$LE_MET, pattern="_", n=6)
dat$LE_MET                   <- Metiers[,1] 

#-----------------------------------------------
# Create table 2
#-----------------------------------------------
# Aggregate data and write to results folder
table2                       <- dat[, .(Land_kg = sum(LE_KG_TOT),
                                        Land_eu = sum(LE_EURO_TOT),
                                        Eff_day = sum(LE_EFF_DAY, na.rm=T),
                                        Eff_kwh = sum(LE_EFF_KWD, na.rm=T)),
                                    by=.(LE_YEAR, SITECODE, LE_MET)]

# Save as csv-file
fwrite(table2, paste0(resPath, "Table_2HAL_", Country, ".csv"))

#-----------------------------------------------




