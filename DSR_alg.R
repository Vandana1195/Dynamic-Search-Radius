# Dynamic Search Radius Algorithm

# Libraries Required ------------------------------------------------------
# for distance calculation using Haversine (distHaversine)
library(geosphere)
library(sp)
library(rgdal)
library(sf)
# for fun.zero.omit function
library(GLDEX)
library(lubridate)
library(stlplus)


# Data -------------------------------------------------------------------------
# TG data filtered and subsetted 
psmsl <- read.csv("Bible_PSMSL1994_2021May_1_July.csv")
# Filtered 
sub_psmsl <- psmsl[psmsl$id !=0,]

# Reference time array from January 1994 to December 2021
time_array =seq(as.Date("1/1/1994", format= "%m/%d/%Y"), by = "month", length.out = 336)

# All XTRACK observation points; XTRACK key file 
Xtrack <- read.csv("xtrack_dataPass.csv", header=TRUE)

# going along the coastline and finding distance to the shelf break
# inputID is the TG station ID, TargetID is the coastal point, 
# DistanceKM is the distance from the TG
along_coastTG <- read.csv("coast_tg.csv",header = TRUE)

# the distance of each of those coastline points(including station points) from the shelf break
coast_to_shelf <- read.csv("coast_to_shelf.csv",header = TRUE)

# target id is the points on the coastline and hubdist is the distance to the shelf break
names(coast_to_shelf)[names(coast_to_shelf) == "ID"] <- "TargetID"

# minor islands list
isleP1 <- read.csv("isl/isl1.csv",header = TRUE)
isleP2 <- read.csv("isl/isl2.csv",header = TRUE)
isleP3 <- read.csv("isl/isl3.csv",header = TRUE)
isleP4 <- read.csv("isl/isl4.csv",header = TRUE)
isleP5 <- read.csv("isl/isl5.csv",header = TRUE)
islandTG <- rbind(isleP1, isleP2, isleP3, isleP4, isleP5)

# Filtering Function for XTRACK-------------------------------------------------------------------------

# inputs are: Xtrack, sub_psmsl, filtering radius= 350 km); 
# the function returns 
filtering_funct <-  function(al, tg, dist_m){
  
  # altimetry observation latitude, longitude
  al_lat <- al$lat
  al_lon <- al$lon
  
  # TG latitude, longitude
  tg_lat <- tg$latitude
  tg_lon <- tg$longitude
  
  # XTRACK file path
  al_fpath <- al$file_name
  
  # extract altimetry observation points
  al_points <- al$nb_points
  
  # code id (with file name and point number)
  al_cid = paste(al$file_name,al$nb_points, sep = "_") 
  
  # distance to nearest GSHHS 1.3 coastline  
  al_dist = al$dist_to_coast
  
  # column binding lon,lat
  coordmat_al = cbind(al_lon, al_lat)
  coordmat_tg= cbind(tg_lon,tg_lat)
  
  # a table with the following columns: code id, latitude, longitude, path address, distance to coast
  finalData_al = cbind(al_cid, al_lat, al_lon, al_fpath, al_points, al_dist)
  
  # declaring an empty matrix with columns = tg stations, rows = altimetry observation point
  f_al = matrix(NA, nrow = nrow(finalData_al), ncol = nrow(tg) )
  
  # logical matrix showing 1 if the point is inside the radius defined by dist_m else 0
  for (i in 1:ncol(f_al)) {
    pdist = distHaversine( coordmat_al, coordmat_tg[i,] )
    ptrue = pdist <= dist_m
    f_al[,i] = ptrue
  }
  
  # sum of all the logical ones for each altimetry observation point; 
  verdict = apply(f_al, 1 , sum)
  # logical vector to flag those points that are not near any tide gauges
  v_gt0 = verdict > 0
  
  # table where each row represents an altimetry point and the column "verdict" lists the 
  # count of TGs whose defined proximity (dist_m) includes that altimetry point
  finalData_al= cbind(finalData_al,verdict)
  
  # for filtering out only those altimetry points which are near TGs
  final_data_req_ind = seq(1, nrow(finalData_al) ) * v_gt0
  final_data_req_ind = fun.zero.omit(final_data_req_ind)
  t_final_data = finalData_al[final_data_req_ind, ]
  t_tg_true = f_al[final_data_req_ind, ]
  
  # trial out matrix has the last column "final_tg_id" listing the TG IDs whose defined 
  # proximity (dist_m) includes that altimetry point
  final_tg_id = apply( t_tg_true, 1, FUN = Tru_id1)
  trial_out = cbind(t_final_data, final_tg_id)
  
  filt_dum = as.data.frame(trial_out)
  return(filt_dum)
}

# function to list out the TG IDs
Tru_id1 = function(x) {
  y = which( x == TRUE )
  
  output = paste( sub_psmsl$id[y], sep = "," )
  
}

# Applying the filtering function; now this is subsetted xtrack key file
xtrack_subupd <- filtering_funct(Xtrack, tg= sub_psmsl, dist_m = 350000)


# Main code -------------------------------------------------------------------------

# to get the unique ids
id_U = xtrack_subupd$final_tg_id
id_unlist_U= as.numeric(unlist(id_U))
id_uniqueUpd= unique(id_unlist_U)

# declaring an empty matrix to store sla time series from altimetry;
# columns: TG station; rows: time
x_sla = matrix(NA, nrow = length(time_array), ncol = nrow(sub_psmsl))

# declaring an empty matrix to store auxilliary information about the altimetry time series
TG_extraInfo = matrix(NA, nrow = nrow(sub_psmsl), ncol = 4)

# stations omitted due to not meeting algorithm criteria
notAt_all_there=c()

# parameter setting
# search radius = 250km
SR_randu = 250
# threshold = 15
Th_s  = 15
# coastline extent = 15km
cst_s = 15

# looping through each TG station to get its corresponding altimetry SLA time series from DSR
for (tg in 1:nrow(sub_psmsl)){
  # checking if any altimetry observations are available near the station
  if (sub_psmsl$id[tg] %in% id_uniqueUpd){
    
    #tg station index; getting the TG station ID
    id_stn <- sub_psmsl$id[tg]
    
    #pattern matching with the TG station ID
    pat_match = gsub(" ", "", paste("\\b",as.character(id_stn), "\\b"))
    
    #subsetting the subsetted xtrack key to find the alt nbpoints having observations around 
    #350km2 from this particular TG(result used is our filtering result)
    sub_pt = xtrack_subupd %>% dplyr::filter(if_any(final_tg_id, ~grepl(pat_match,.)))
    
    # points along the coastline near the station (nearest 500 points)
    tg_sub <- along_coastTG %>% dplyr::filter(InputID==id_stn)
    
    # case 1: if the tg station is on an island
    if(id_stn %in% islandTG$id){
      
      # if the TG station is on an island; it will be seeing/capturing the open ocean signals since it is on an island
      shelf_break=200
    }else{
      # shelf break is either 200km (for islands) or max shelf break (for stations in the continental land mass)
      # parameter 1: coastaline extent "cst_s"; we take points along the coastline extending 15km in
      # in either direction of the TG
      tg_sub2 <- tg_sub[tg_sub$Distancekm<=cst_s,]
      com_id <- merge(coast_to_shelf,tg_sub2, by="TargetID")
      # maximum shelf width (distance to the shelf break) is chosen as SR1 (denoted as shelf_break)
      shelf_break= max(com_id$HubDist)
      
    }
    # auxiliary data storing; column 2: SR1
    TG_extraInfo[tg,2]= shelf_break
    
    # no relevant data; station omitted
    if(nrow(com_id)==0){
      paste_thg <- paste("stn_not_taken",id_stn)
      # prints the omitted station id 
      print(paste_thg)
      notAt_all_there = c(notAt_all_there, id_stn)
      TG_extraInfo[tg, 1] = id_stn
      TG_extraInfo[tg, 2] = NA
      TG_extraInfo[tg, 3] = NA
      TG_extraInfo[tg, 4] = NA
      next
    }
    
    # list of latitudes of the relevant altimetry measurements for this particular TG 
    alLat <- as.double(unlist(sub_pt$al_lat))
    # list of longitudes of the relevant altimetry measurements for this particular TG
    alLon <- as.double(unlist(sub_pt$al_lon))
    
    #al lat lon combined vector; contains all pairs of lat lon corresponding to the alt observations
    coordmatAL <- cbind(alLon,alLat)
    #tg lat lon combined vector
    coormatTG <- cbind(sub_psmsl$longitude[tg], sub_psmsl$latitude[tg])
    
    f_al <- matrix(NA, nrow=nrow(sub_pt), ncol=1)
    
    # filtering observations 
    for(i in 1: nrow(sub_pt)){
      pdist= distHaversine(coordmatAL, coormatTG[1,1:2])
      ptrue <- (pdist/1000) <=shelf_break
      f_al[,1] <- ptrue
    }
    sub_pt2 <- cbind(sub_pt,f_al)
    
    #final subset of relevent altimetry observation for the validation/comparison
    sub_stage2_1 <- sub_pt2 %>% dplyr::filter(f_al==TRUE)
    
    # parameter 2: threshold
    # we set a lower limit of 15 altimetry points to be present within SR1
    # if the condition is not met; we use SR2 (second search radius)
    if(nrow(sub_stage2_1)<Th_s && !(id_stn %in% islandTG)){
      
      SR = SR_randu
      # auxiliary data; column 3: SR2
      TG_extraInfo[tg, 3] = SR_randu
      
      alLat <- as.double(unlist(sub_pt$al_lat))
      alLon <- as.double(unlist(sub_pt$al_lon))
      
      coordmatAL <- cbind(alLon,alLat)
      coormatTG <- cbind(sub_psmsl$longitude[tg], sub_psmsl$latitude[tg])
      
      f_al <- matrix(NA, nrow=nrow(sub_pt), ncol=1)
      
      # filtering observations with SR2
      for(i in 1: nrow(sub_pt)){
        pdist= distHaversine(coordmatAL, coormatTG[1,1:2])
        ptrue <- (pdist/1000) <=SR_randu
        f_al[,1] <- ptrue
      }
      sub_pt2 <- cbind(sub_pt,f_al)
      # final subset of the relevent altimetry observation after filtering
      sub_stage2_inter <- sub_pt2 %>% dplyr::filter(f_al == TRUE)
      # altimetry points only within the shelf width inside SR2 will be filtered
      sub_stage2_1 <- sub_stage2_inter %>% dplyr::filter(al_dist<=shelf_break)
      
    }
    
    # if threshold is not met, the station is omitted
    if(nrow(sub_stage2_1)<Th_s&& !(id_stn %in% islandTG)){
      notAt_all_there <- c(notAt_all_there, id_stn)
      next
    }
    
    # auxilliary information: column 3: maximum shelf width over the extent
    TG_extraInfo[tg, 3] = shelf_break
    # auxilliary information: column 1: TG station ID
    TG_extraInfo[tg,1]= id_stn
    
    # getting the unique files (each xtrack file corresponds to specific satellite mission,
    # region and track number. note: we use TOPEX/Poseidon+Jason1+Jason2+Jason3 data)
    path_list_1 <- unique(sub_stage2_1$al_fpath)
    
    # total number files having relevant altimetry points 
    l=length(path_list_1)
    # initialising an empty list
    x_tg_1= list()
    
    # initialising an empty matrix to store altimetry SLA time series at all the selected points
    inter_alt <- matrix(NA, nrow=336, ncol=l)
    x_tgP= list()
    x_df= data.frame()
    
    #going through each file one by one  
    for(p in 1:l){
      
      # auxilliary data; column 4: total number of satellite passes contributing to the final
      # altimetry sla time series
      TG_extraInfo[tg, 4] = p
      
      x_tg= list()
      var_2 <- strsplit(path_list_1[[p]], 'J3')
      var_3 = var_2[[1]][3]
      
      # subsetting for a particular track; all the relevant altimetry points 
      sub_2 = sub_stage2_1 %>% dplyr::filter(if_any(al_fpath, ~grepl(var_3,.)))
      
      xy_name = gsub('[.]','',var_3)
      
      csv_f1 = paste(substring(xy_name,1,nchar(xy_name)-5), substring(xy_name,nchar(xy_name)-4, nchar(xy_name)-2), sep='_')
      # getting the altimetry path .nc file
      file_name1 = paste("xtrack_sla_350_94_21_month/",csv_f1,".csv")
      
      file_name1 = gsub(" ","",file_name1)
      
      # altimetry path file; rows: time, columns: lat,lon point;
      # 1st column: sl.no; 2nd column: date column; data from 3rd column
      alt_matrix = read.csv(file_name1, header=TRUE)
      alty = alt_matrix
      
      #fixing the time array
      start_date = alty$xtrack_dummymonth.month[1]
      
      time_period = length(alty$xtrack_dummymonth.month)
      
      end_date = alty$xtrack_dummymonth.month[time_period]
      
      date_array = seq(ymd(start_date),ymd(end_date), by = '1 month')
      
      ind_start <- which(date_array==time_array[1])
      ind_end <- which(date_array==time_array[length(time_array)])
      
      for (f in 1:nrow(sub_2)){
        col_id <- as.numeric(sub_2$al_points[f])
        col_f <- alty[ind_start:ind_end, col_id+2]
        
        x_tg[[f]] = col_f
      }
      
      x_tgP[[p]] = x_tg
      
    }
    x_df <- as.data.frame(x_tgP)
    
    # averaging over all relevant altimetry points associated with that TG station
    x_sla[,tg] <- rowMeans(x_df,na.rm=TRUE)
    
  }
}

# saving 
# altimetry SLA time series: columns: time; rows: altimetry SLA associated to each TG station 
write.csv(x_sla,"XTRACK_sla_250km_15th_15km.csv",row.names = FALSE)
# auxilliary data: column 1: TG ID; column 2: SR1; column 3: maximum shelf width; column 4: total number of satellite passes used
write.csv(TG_extraInfo, "auxilliary_data.csv", row.names = FALSE)
  
    
    
    

