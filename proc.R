#install.packages(c("rsconnect"),repos='http://cran.r-project.org',lib="/home/kcm/LIBS/R")
#.libPaths("/home/kcm/LIBS/R")
.libPaths("/usr/local/lib/R/site-library")

#rmarkdown::draft("dashboard.Rmd",template="flex_dashboard",package="flexdashboard")


#Rrsconnect::setAccountInfo(name='kkobchul', token='6C05E2EEB1D428475D34F5F2C472D342', secret='lJXTvy53IuQOzUl+uvi0b8a9rAUeME4iqsUELIrd')
#rsconnect::deployApp('path/to/your/app')

#rmarkdown::run("index.Rmd")
rmarkdown::render_site()

#df<-df0 %>% filter((lon>=(min(nNewLon)-1) & lon<=(max(nNewLon)+1) & lat>=(min(nNewLat)-1) & lat<=(max(nNewLat)+1)) & int==tt)
#dfTmp<-akima::interp(df$lon,df$lat,df$Wspd,xo=nNewLon,yo=nNewLat)
#
#spGeoData = expand.grid(nNewLon, nNewLat)
#coordinates(spGeoData) = ~ Var1 + Var2
#gridded(spGeoData) = TRUE
#
#dfDataL2<-data.frame(crtn_tm = df$crtn_tm[1],fcst_tm = df$fcst_tm[1],int = df$int[1],nLon = spGeoData$Var1, nLat = spGeoData$Var2, nVal = c(dfTmp$z)) %>% mutate(nVal=round(nVal,1))

###### Power curve MODULE
#curve1<-data.frame()
#for (i in c(1:length(SheetList))){
#tmp<-SheetList[[i]]
#tmp1<-str_split(names(tmp),pattern=' ',simplify = TRUE)
#curve1<-rbind(curve1,data.frame(manu=tmp1[2,1],model=tmp1[1,1],check.names=FALSE))
#curve<-data.frame(manu=tmp1[2,1],model=tmp1[1,1],wspd=as.numeric(tmp[,1]),wp=as.numeric(tmp[,2])) %>% mutate(wspd=ifelse(is.na(wspd),0,wspd),wp=ifelse(is.na(wp),0,wp))
#print(curve)
##varname<-paste0(tmp1[2,1],"_",tmp1[1,1])
##dfDataL2<-dfDataL2 %>% mutate(!!varname := round(pracma::interp1(curve$wspd,curve$wp,xi=nVal,method="linear"),1))
#}
#wp1<-read_csv("/home/kcm/PROJ/VPPLAB/data/wp_manu_info_data.csv")
#wp2<-read_csv("/home/kcm/PROJ/VPPLAB/data/wp_perf_info_data.csv")
#wp1<-wp1 %>% dplyr::select(MANU_ID,MANU_NM)
#wp2<-wp2 %>% dplyr::select(MODEL_ID,MODEL_NM,RTD_CAPA,RTR_DMTR,MANU_ID,END_WS,AIR_DENSITY,POWER_CURVE_JSON) 
#wp3<-left_join(wp2,wp1,by=c("MANU_ID")) %>% filter(str_detect(MANU_NM,"Doosan|GAMESA|Unison|Siemens|VESTAS")) 
#curve<-data.frame()
#for (i in c(1:nrow(wp3))){
#wp4<-wp3[i,]
#tmp<-cbind(wp4 %>% dplyr::select(-POWER_CURVE_JSON),fromJSON(wp4$POWER_CURVE_JSON)) %>% dplyr::select(MANU_NM,MANU_ID,MODEL_NM,MODEL_ID,RTD_CAPA,RTR_DMTR,END_WS,AIR_DENSITY,WSPD,QGEN)
#curve<-rbind(curve,tmp)
#}

#dfDataL2<-dfDataL2 %>% mutate(lon=nLon,lat=nLat) %>% st_as_sf(coords=c("lon","lat"),crs=4326) %>% st_intersection(shp1) %>% st_set_geometry(NULL) %>% dplyr::select(-area)
####### plot module
#
#cbMatlab =  colorRamps::matlab.like(11)
#mapKor = read_sf("gadm36_KOR_shp/gadm36_KOR_1.shp")
#
#ggplot() + 
#coord_fixed(ratio = 1.1) +
#theme_bw() +
#geom_tile(data = dfDataL2, aes(x = nLon, y = nLat, fill = nVal)) +
##geom_tile(data = dfDataL2, aes(x = nLon, y = nLat, fill = Doosan_WinDS5560)) + 
#scale_fill_gradientn(colours = cbMatlab, limits=c(0, 10), breaks = seq(0, 10, 1), na.value = cbMatlab[length(cbMatlab)]) +
##scale_fill_gradientn(colours = cbMatlab, limits=c(0, 1000), breaks = seq(0, 1000, 100), na.value = cbMatlab[length(cbMatlab)]) +
#geom_sf(data = mapKor, color = "black", fill = NA) +
#metR::scale_x_longitude(expand = c(0, 0), breaks = seq(125.5, 127.5, 0.5), limits = c(125.5, 127.5)) +
#metR::scale_y_latitude(expand = c(0, 0), breaks = seq(32.5, 34, 0.5), limits = c(32.5, 34)) +
#labs(x="",y="",fill="",colour="",title="") +
#theme(
#    plot.title = element_text(face = "bold", size = 16, color = "black")
#    , axis.title.x = element_text(face = "bold", size = 16, colour = "black")
#    , axis.title.y = element_text(face = "bold", size=16, colour = "black", angle = 90)
#    , axis.text.x = element_text(face = "bold", size=16, colour = "black")
#    , axis.text.y = element_text(face = "bold", size=16, colour = "black")
#    , legend.position = c(1, 1.03)
#    , legend.justification = c(1, 1)
#    , legend.key = element_blank()
#    , legend.text = element_text(size = 14, face = "bold")
#    , legend.title = element_text(face = "bold", size = 14, colour = "black")
#    , legend.background=element_blank()
#    , plot.margin = unit(c(0, 8, 0, 0), "mm")
#) 
#df_new1<-rbind(df_new1,dfDataL2)
#}
#
#head(df_new1)
#df_new2<-df_new1 %>% group_by(nLon,nLat) %>% summarize(wnd=paste0(as.character(nVal),collapse=","),wp1=paste0(as.character(Doosan_WinDS5560),collapse=","),wp2=paste0(as.character(GE_GE.5500.158),collapse=","),wp3=paste0(as.character(Hyosung_HS139.5.5MW),collapse=","),wp4=paste0(as.character(Vestas_V150.5.6MW),collapse=","),wp5=paste0(as.character(Unison_U136.4.2MW),collapse=","),wp6=paste0(as.character(Unison_U151.4.3MW),collapse=","),wp7=paste0(as.character(Vestas_V90.3MW),collapse=","),wp8=paste0(as.character(GE_Haliade.150.6MW),collapse=","),wp9=paste0(as.character(Siemens_SWT.8.0.167),collapse=","),wp10=paste0(as.character(Enercon_E.82.E3),collapse=","),wp11=paste0(as.character(Vestas_V80.2MW),collapse=",")) %>% ungroup() %>% dplyr::select(-nLon,-nLat)

####### DB create & write module(sudo service mysql start)
#con <- dbConnect(MySQL(),username="kcm",password="rhcjfals",host="localhost",port=3306,dbname="wind_power",client.flag=CLIENT_MULTI_RESULTS)
##dbSendQuery(con,"CREATE TABLE grid_info(gid INT NOT NULL AUTO_INCREMENT, 
##                                        lon DOUBLE,
##                                        lat DOUBLE,
##                                        PRIMARY KEY (gid));")
##dbSendQuery(con,"CREATE TABLE time_info(tid INT NOT NULL AUTO_INCREMENT,
##                                        crtn_time DATETIME,
##                                        PRIMARY KEY (tid));")
##dbSendQuery(con,"CREATE TABLE turb_info(pid INT NOT NULL AUTO_INCREMENT,
##                                        manu VARCHAR (20), 
##                                        modl VARCHAR (20),
##                                        PRIMARY KEY (pid));")
##dbSendQuery(con,"CREATE TABLE powr_info(tid INT NOT NULL,
##                                        gid INT NOT NULL,
##                                        wnd TEXT, wp1 TEXT, wp2 TEXT, wp3 TEXT, wp4 TEXT, wp5 TEXT, wp6 TEXT, wp7 TEXT, wp8 TEXT, wp9 TEXT, wp10 TEXT, wp11 TEXT);")
#
##dbWriteTable(con,"grid_info",data.frame(lon=dfDataL2$nLon,lat=dfDataL2$nLat),overwrite=F,append=T,row.names=F)
##dbWriteTable(con,"turb_info",data.frame(manu=curve1$manu,modl=curve1$model),overwrite=F,append=T,row.names=F)
#
#dbWriteTable(con,"time_info",data.frame(crtn_time=crtn_time),overwrite=F,append=T,row.names=F)
#tid<-dbGetQuery(con, "SELECT LAST_INSERT_ID() FROM time_info;")[1,1]
#dbWriteTable(con,"powr_info",data.frame(tid=tid,gid=c(1:nrow(dfDataL2)),df_new2),overwrite=F,append=T,row.names=F)
#
#dbDisconnect(con)

















#e <- extent(123.3102, 132.7750, 31.6518, 43.3935)
#p <- as(e, "SpatialPolygons")
#crs(p) <- "+proj=longlat +datum=WGS84" 
#p<-p %>% st_as_sf
#
#pp <- p %>% st_transform(lcc)
#pp
#r <- raster(pp)
#res(r)<-5000
#r

#grid<-raster(extent(shp))
#res(grid)<-5
#proj4string(grid)<-proj4string(shp)
#gridpolygon<-rasterToPolygons(grid)
#grid_new<-intersect(shp,gridpolygon)
#plot(grid_new)


#coords_sf <- st_as_sf(grd, grd = c("lon", "lat"), crs = 4326)
## Then you can transform them to the system you want
#coords_shp <- coords_sf %>% st_transform(crs = st_crs(shp))
#coords_shp
#shp %>% st_contains(coords_shp)




###################################################
#
#get_wether_forecast_nwp <- function(type, base, lead, var, service_key){
#
#if (type == "LDAPS") {
#    url <- "http://apis.data.go.kr/1360000/NwpModelInfoService/getLdapsUnisAll"
#    res <- GET(url = url,
#               query = list(ServiceKey=service_key,
#                            numOfRows="10",
#                            dataType="JSON",
#                            baseTime=paste0(substr(base,1,4),substr(base,6,7),substr(base,9,10),substr(base,12,13),substr(base,15,16)),
#                            leadHour=lead,
#                            dataTypeCd=var))
#
#} else if (type == "RDAPS") {
#    url <- "http://apis.data.go.kr/1360000/NwpModelInfoService/getRdapsUnisAll"
#    res <- GET(url = url,
#               query = list(ServiceKey=service_key,
#                            numOfRows="10",
#                            dataType="JSON",
#                            baseTime=paste0(substr(base,1,4),substr(base,6,7),substr(base,9,10),substr(base,12,13),substr(base,15,16)),
#                            leadHour=lead,
#                            dataTypeCd=var))
#}
#
#    res_json <- res %>% content(as="text") %>% fromJSON()
#    fcst_df <- res_json$response$body$items$item
#
#if (type == "LDAPS") {    
#    nwp_grd<-read.csv("nwp_ldps_grd.csv",header=T)
#} else if (type == "RDAPS") {
#    nwp_grd<-read.csv("nwp_rdps_grd.csv",header=T)
#}
#    cbind(nwp_grd,val=as.numeric(unlist(str_split(fcst_df$value,pattern=","))))
#
#}
#
#############################################################################
#
#get_wether_forecast_dfs <- function(type, ix, iy, datetime, service_key){
#
#if (type == "SHRT") {
#    url <- "http://apis.data.go.kr/1360000/VilageFcstInfoService_2.0/getVilageFcst"
#} else if (type == "VSRT") {
#    url <- "http://apis.data.go.kr/1360000/VilageFcstInfoService_2.0/getUltraSrtFcst"
#} else if (type == "ODAM") {
#    url <- "http://apis.data.go.kr/1360000/VilageFcstInfoService_2.0/getUltraSrtNcst"
#}
#
#    res <- GET(url = url,
#               query = list(ServiceKey=service_key,
#                            numOfRows="1000",
#                            dataType="JSON",
#                            base_date=paste0(substr(datetime,1,4),substr(datetime,6,7),substr(datetime,9,10)),
#                            base_time=paste0(substr(datetime,12,13),substr(datetime,15,16)),
#                            nx=ix,
#                            ny=iy))
#
#
#    res_json <- res %>% content(as="text") %>% fromJSON()
#    fcst_df <- res_json$response$body$items$item
#
#if (type == "SHRT") {
#    category <- c("TMP", "REH", "WSD", "VEC", "SKY", "PCP")
#    fcst_df %>%
#      filter(category %in% !!category) %>%
#      pivot_wider(id_cols=c("fcstDate", "fcstTime"),names_from=category,values_from=fcstValue) %>%
#      mutate(fcst_time=as.POSIXct(paste0(fcstDate,fcstTime),format="%Y%m%d%H%M")) %>% select(-fcstTime, -fcstDate) %>%
#      select(fcst_time, everything()) %>% mutate(PCP=ifelse(PCP=="강수없음","0",PCP)) 
#
#} else if (type == "VSRT") {
#    category <- c("T1H", "REH", "WSD", "VEC", "SKY", "RN1")
#    fcst_df %>%
#      filter(category %in% !!category) %>%
#      pivot_wider(id_cols=c("fcstDate", "fcstTime"),names_from=category,values_from=fcstValue) %>%
#      mutate(fcst_time=as.POSIXct(paste0(fcstDate,fcstTime),format="%Y%m%d%H%M")) %>% select(-fcstTime, -fcstDate) %>%
#      select(fcst_time, everything()) %>% rename(TMP=T1H,PCP=RN1) %>% mutate(PCP=ifelse(PCP=="강수없음","0",PCP))
#
#} else if (type == "ODAM") {
#    category <- c("T1H", "REH", "WSD", "VEC", "RN1")
#    fcst_df %>%
#      filter(category %in% !!category) %>%
#      pivot_wider(id_cols=c("baseDate", "baseTime"),names_from=category,values_from=obsrValue) %>%
#      mutate(fcst_time=as.POSIXct(paste0(baseDate,baseTime),format="%Y%m%d%H%M")) %>% select(-baseTime, -baseDate) %>%
#      select(fcst_time, everything()) %>% rename(TMP=T1H,PCP=RN1) 
#}
#
#}
##########################################################
#
#service_key1 <- readLines("./service_key1.txt")
#service_key2 <- readLines("./service_key2.txt")
#
##ctrn_time<-as.POSIXct("2022-03-11 11:00:00")
###df1<-get_wether_forecast_dfs("ODAM",55,127,ctrn_time,service_key2)
###df2<-get_wether_forecast_dfs("VSRT",55,127,ctrn_time,service_key2)
##
##df<-get_wether_forecast_dfs("SHRT",55,127,ctrn_time,service_key2)
##cbind(ctrn_time,df) %>% mutate(leadTime=as.numeric(fcst_time-ctrn_time,units="hours")) %>% as.data.frame
#
#df<-get_wether_forecast_nwp("LDAPS","2022-03-11 03:00:00",0,"Temp",service_key1)
#df


#df1<-df %>% mutate(val=ifelse(val==-9999,NA,val)) %>% na.omit
#fld <- with(df1, interp(x = lon, y = lat, z = val, ))
#filled.contour(x = fld$x,
#               y = fld$y,
#               z = fld$z,
#               color.palette =
#                 colorRampPalette(c("white", "blue")),
#               xlab = "Longitude",
#               ylab = "Latitude",
#               main = "Rwandan rainfall",
#               key.title = title(main = "Rain (mm)", cex.main = 1))

#df <- melt(fld$z, na.rm = TRUE)
#names(df) <- c("x", "y", "Rain")
#df$Lon <- fld$x[df$x]
#df$Lat <- fld$y[df$y]
#
#ggplot(data = df, aes(x = Lon, y = Lat, z = Rain)) +
#  geom_tile(aes(fill = Rain)) +
#  stat_contour() +
#  ggtitle("Rwandan rainfall") +
#  xlab("Longitude") +
#  ylab("Latitude") +
#  scale_fill_continuous(name = "Rain (mm)",
#                        low = "white", high = "blue") +
#  theme(plot.title = element_text(size = 25, face = "bold"),
#        legend.title = element_text(size = 15),
#        axis.text = element_text(size = 15),
#        axis.title.x = element_text(size = 20, vjust = -0.5),
#        axis.title.y = element_text(size = 20, vjust = 0.2),
#        legend.text = element_text(size = 10))














#get_wether_forecast_vsrt <- function(ix, iy, base, service_key){
#
#    url <- "http://apis.data.go.kr/1360000/VilageFcstInfoService_2.0/getUltraSrtFcst"
#
##    base_date <- strftime(today(), format="%Y%m%d") # 오늘날짜 자동 입력 ex) '20210509'
##    region_xy <- list(dangjin=c("53", "114"),
##                      ulsan=c("102", "83"),
##                      gujwa=c("56", "39"),
##                      pyoseon=c("56","35"))
#    res <- GET(url = url,
#               query = list(ServiceKey=service_key,
#                            numOfRows="100",
#                            dataType="JSON",
#                            base_date=paste0(substr(base,1,4),substr(base,6,7),substr(base,9,10)),
#                            base_time=paste0(substr(base,12,13),substr(base,15,16)),
##                            base_date=base_date,
##                            base_time=base_time,            # 14시 예보
#                            nx=ix,
#                            ny=iy))
##                            nx=region_xy[[region]][1],
##                            ny=region_xy[[region]][2]))
#
#
#    res_json <- res %>% content(as="text") %>% fromJSON()
#    fcst_df <- res_json$response$body$items$item
#
#    category <- c("REH", "T1H", "VEC", "WSD", "SKY") # 사용할 예보 종류
#
#    fcst_df %>%
#      filter(category %in% !!category) %>% #,
#             #fcstDate != base_date) %>%
#      pivot_wider(id_cols=c("fcstDate", "fcstTime"),
#                  names_from=category,
#                  values_from=fcstValue) %>%
#      rename(humidity=REH,
#             cloud=SKY,
#             temperature=T1H,
#             winddirection=VEC,
#             windspeed=WSD) %>%
#      mutate(forecast_time=as_datetime(glue("{fcstDate} {fcstTime}"), format="%Y%m%d %H%M")) %>%
#      select(-fcstTime, -fcstDate) %>%
#      select(forecast_time, everything()) %>% as.data.frame
#
#}
#
##############################################################
#
#get_wether_forecast_shrt <- function(ix, iy, base, service_key){
#    
#    url <- "http://apis.data.go.kr/1360000/VilageFcstInfoService_2.0/getVilageFcst"
#
##    base_date <- strftime(today(), format="%Y%m%d") # 오늘날짜 자동 입력 ex) '20210509'
##    region_xy <- list(dangjin=c("53", "114"),
##                      ulsan=c("102", "83"),
##                      gujwa=c("56", "39"),
##                      pyoseon=c("56","35"))
#
#    res <- GET(url = url,
#               query = list(ServiceKey=service_key,
#                            numOfRows="1000",
#                            dataType="JSON",
#                            base_date=paste0(substr(base,1,4),substr(base,6,7),substr(base,9,10)),
#                            base_time=paste0(substr(base,12,13),substr(base,15,16)),
##                            base_date=base_date,
##                            base_time=base_time,            # 14시 예보 
#                            nx=ix,
#                            ny=iy))
##                            nx=region_xy[[region]][1],
##                            ny=region_xy[[region]][2]))
#    
#
#    res_json <- res %>% content(as="text") %>% fromJSON()
#    fcst_df <- res_json$response$body$items$item 
#
#    category <- c("REH", "TMP", "VEC", "WSD", "SKY") # 사용할 예보 종류
#
#    fcst_df %>%
#      filter(category %in% !!category) %>% #,
#             #fcstDate != base_date) %>%
#      pivot_wider(id_cols=c("fcstDate", "fcstTime"),
#                  names_from=category,
#                  values_from=fcstValue) %>%
#      rename(humidity=REH,
#             cloud=SKY,
#             temperature=TMP,
#             winddirection=VEC,
#             windspeed=WSD) %>%
#      mutate(forecast_time=as_datetime(glue("{fcstDate} {fcstTime}"), format="%Y%m%d %H%M")) %>%
#      select(-fcstTime, -fcstDate) %>%
#      select(forecast_time, everything()) %>% as.data.frame
#
#}
#
##################################################################
#get_wind_speed_by_altitude <- function(pid, tid, ix, iy, rl, base){
#
#df1<-get_wether_forecast_vsrt(ix, iy, format(strptime(base,"%Y-%m-%d %H:%M:%S")-3600*0.5,"%Y-%m-%d %H:%M:%S"), service_key)
#df2<-get_wether_forecast_shrt(ix, iy, base, service_key)
#     
#     bind_rows(df1,df2) %>% 
#     distinct(.,forecast_time,.keep_all=TRUE) %>% 
#     select(windspeed,winddirection) %>% 
#     rename(ws10=windspeed,wd10=winddirection) %>% 
#     mutate(ws10=as.numeric(ws10),wd10=as.numeric(wd10),ws10=ifelse(ws10<=0,0.1,ws10)) %>% 
#     mutate(ws30=ws10*(30/10)^((1/(log(sqrt(30*10)/rl))+0.088/(1-0.088*log(10/10)))+((-0.088/(1-0.088*log(10/10)))*log(ws10))), # deacon eq.
#            ws50=ws10*(50/10)^((1/(log(sqrt(50*10)/rl))+0.088/(1-0.088*log(10/10)))+((-0.088/(1-0.088*log(10/10)))*log(ws10))),
#            ws80=ws10*(80/10)^((1/(log(sqrt(80*10)/rl))+0.088/(1-0.088*log(10/10)))+((-0.088/(1-0.088*log(10/10)))*log(ws10)))) %>%
##        mutate(we10=(1/2*(1.23)*(pi*(52)^2)*(ws10^3)*0.4)/10^3, # units: kwh
##               we30=(1/2*(1.23)*(pi*(52)^2)*(ws30^3)*0.4)/10^3,
##               we50=(1/2*(1.23)*(pi*(52)^2)*(ws50^3)*0.4)/10^3,
##               we80=(1/2*(1.23)*(pi*(52)^2)*(ws80^3)*0.4)/10^3) %>%
#     mutate_if(is.numeric,round,1) %>%
#     summarize(ws10=paste(ws10,collapse=","),ws30=paste(ws30,collapse=","),ws50=paste(ws50,collapse=","),ws80=paste(ws80,collapse=","),
#               wd10=paste(wd10,collapse=","),pid=pid,tid=tid)
##                 we10=paste(we10,collapse=","),we30=paste(we30,collapse=","),we50=paste(we50,collapse=","),we80=paste(we80,collapse=","),
#}
#
#################################################################
#
#if (!file.exists("stn_info.csv")){
#stn<-read.xlsx("/mnt/d/랩이앤이/회사사업/발전소_포인트위치/211201_발전량 산출포인트_200지점.xlsx",sheet=1) %>% select(1:5) %>% rename(pid=X1,sido=도,sigungu=시군,lat=N,lon=E) %>% select(lon,lat,sido,sigungu) %>% na.omit %>% mutate(pid=c(1:length(lon))) %>% select(pid,lon,lat,sido,sigungu)
#
#xy<-list()
#for (i in c(1:nrow(stn))){
#xy<-append(xy,dfs_xy_conv("toXY",stn$lat[i],stn$lon[i]))
#}
#
#stn<-bind_cols(stn,data.frame(matrix(unlist(xy),ncol=2,byrow=TRUE))) %>% rename(ix=X1,iy=X2) %>% select(pid,lon,lat,ix,iy,sido,sigungu)
#
#nc<-nc_open("geo_em.d01.nc")
#lon<-ncvar_get(nc,"XLONG_M")
#lat<-ncvar_get(nc,"XLAT_M")
#lu_index<-ncvar_get(nc,"LU_INDEX")
#geo<-data.frame(lon1=c(lon),lat1=c(lat),lu_index=c(lu_index))
#landuse<-read_csv("LANDUSE.TBL",col_names=FALSE,show_col_types = FALSE) %>% select(X1,X5) %>% rename(id=X1,rl=X5) %>% mutate(id=as.numeric(id),rl=rl/100,rl=ifelse(rl<=0.1,0.1,rl))
#
#dist<-geodist(stn[,c("lon","lat")],geo[,c("lon1","lat1")],paired=FALSE,sequential=FALSE,pad=FALSE,measure="haversine")
#stn1<-data.frame(stn,lu_index=geo[apply(dist,1,which.min),c("lu_index")],dist=round(apply(dist,1,min)/1000,3))
#
#stn1<-left_join(stn1,landuse,by=c("lu_index"="id"))
#dbWriteTable(con,"wt_region",stn,row.names=F,append=FALSE,overwrite=TRUE)
#write_csv(stn1,"stn_info.csv")
#} else {
#stn1<-read_csv("stn_info.csv",col_names=TRUE,show_col_types = FALSE)
#}
#
#####################################################################
#
#dt0<-"2021-12-29 05:00:00"
#dt<-data.frame(crtn_time=dt0)
#dbWriteTable(con,"wt_date",dt,row.names=F,append=TRUE,overwrite=FALSE)
#tid<-dbGetQuery(con, "SELECT LAST_INSERT_ID() FROM wt_date;")[1,1]
#
#for (i in c(1:nrow(stn1))){
#df<-get_wind_speed_by_altitude(i,tid,stn1$ix[i],stn1$iy[i],stn1$rl[i],dt$crtn_time)
#print(df)
#dbWriteTable(con,"wt_energy",df,row.names=F,append=TRUE,overwrite=FALSE)
#}
#
#
#dbDisconnect(con)






















#con <- dbConnect(MySQL(),username="kcm",password="rhcjfals",host="localhost",port=3306,dbname="testdb",client.flag=CLIENT_MULTI_RESULTS)


#format_future_fcst <- function(fcst_data) {
#    data_dummy <- data.frame(forecast_time=seq(min(fcst_data$forecast_time),
#                                       max(fcst_data$forecast_time),
#                                       by="1 hour"))
#
#data_dummy %>%
#left_join(fcst_data) %>%
#arrange(forecast_time) %>%
#    mutate_if(is.character, type.convert) %>%
#    mutate_if(is.numeric, na.approx) %>%
#    slice(-1)
#}

#
##postscriptFonts()
##pdf.options(family = "Korea1deb")
#
#df<-read.xlsx('KIMST 전문가명부 현황_업데이트_(주)지니_취합본(211206)_최종본.xlsx',sheet=1,startRow=3) 
#df<-df[,c(1,2,4,5,8:15)]
#colnames(df)<-c("num","name","deg","div","main_cat_2018","mid_cat_2018","sub_cat_2018","etc_cat_2018","main_cat_2017","mid_cat_2017","sub_cat_2017","etc_cat_2017")
#
#df[df=="-" | df==" "]<-NA
#
#dir0<-c("/mnt/d/기술품질평가센터/회사사업/해양용역/전문가명부_분석/")
#
############################################
#
##df1<-df %>% select(num,name,deg) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(deg))) 1 else which(!is.na(deg))) %>% distinct() %>% as.data.frame
##df2<-df1 %>% group_by(deg) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
##df20<-df2 %>% mutate(percent=round(total/sum(total)*100,2)) %>% mutate(deg=ifelse(is.na(deg),"확인불가",deg)) %>% as_tibble 
##write.xlsx(df20,paste0(dir0,"학위비율_공란포함.xlsx"))
##df21<-df2 %>% na.omit %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
##write.xlsx(df21,paste0(dir0,"학위비율_공란미포함.xlsx"))
#
##df1<-df %>% select(num,name,div) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(div))) 1 else which(!is.na(div))) %>% distinct() %>% as.data.frame
##df2<-df1 %>% group_by(div) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
##df20<-df2 %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
##write.xlsx(df20,paste0(dir0,"소속기관비율_기타포함.xlsx"))
#
##############################################################
#df1<-df %>% mutate(main_cat_2018=str_replace_all(main_cat_2018,'[\r\n]| ',''),
#                   main_cat_2018=str_replace_all(main_cat_2018,'EB]재료','EB]'),
#                   main_cat_2018=str_replace_all(main_cat_2018,'EB]','EB]재료'),
#                   main_cat_2018=str_replace_all(main_cat_2018,'OC]과학기술과인물사회','OC]과학기술과인문사회'),
#                   main_cat_2018=str_replace_all(main_cat_2018,'해양안전/교통','[EI]건설/교통'),
#
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'[\r\n]| ',''),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'EA14]재난안전장비','EA14]'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'EA14]','EA14]재난안전장비'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'EE08]홈네크워크','EE08]홈네트워크'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'EE12ITS/]텔레매틱스','EE12]ITS/텔레매틱스'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'EI06]해양환경','EI06]철도교통기술'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'LC15]독성/안전성관리기반기반기술','LC15]독성/안전성관리기반기술'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'NC0]유기화학','NC02]유기화학'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'ND08해양과학','ND08]해양과학'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'SB12]분야별/유형별/행정/정책','SB12]분야별/유형별행정/정책'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'SC11]경영정보/E-비즈니스','SC11]경영정보/e-비즈니스'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'LA03]발생/신경생물학','[LA03]발생/신경생물학'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'NB99]기타물리학','[NB99]기타물리학'),
#                   mid_cat_2018=str_replace_all(mid_cat_2018,'\\[\\[','\\[')#,
#
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'[\r\n]| ',''),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'EA14]재난안전장비','[EA1408]위험감지/모니터링장비'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'EB103]복합재료','EB0103]복합재료'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'EE0202]SW솔루션','EE0202]S/W솔루션'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'EE0599]달리분류되지안흔위성/전파','EE0599]달리분류되지않는위성/전파'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'EE0802.유/무선홈네트워킹기술','EE0802]유/무선홈네트워킹기술'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'EH1405]환경친화적제품설계기술','EH1405]환경친화적제품설계기'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'EH1405]환경친화적제품설계기','EH1405]환경친화적제품설계기술'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'LC0102]생화학','LC0103]생화학'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'ND0103.광상/자원지질학]','ND0103.광상/자원지질학'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'ND0103.광상/자원지질학','ND0103]광상/자원지질학'),
##                   sub_cat_2018=str_replace_all(sub_cat_2018,'ND1002]해양수산생물자원유전현상규명','ND1002]해양생물자원유전현상규명')
#
#)
#
#df2<-df1 %>% select(num,name,main_cat_2018) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(main_cat_2018))) 1 else which(!is.na(main_cat_2018))) %>% distinct() %>% as.data.frame
#head(df2,n=30)
#df2<-df2 %>% group_by(main_cat_2018) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
#head(df2,n=30)
#df20<-df2 %>% mutate(main1_cat_2018=ifelse(is.na(main_cat_2018),"DBX","DBO")) %>% group_by(main1_cat_2018) %>% summarize(total=sum(total)) %>% mutate(percent=round(total/sum(total)*100,2))
#write.xlsx(df20,paste0(dir0,"국가과학기술표준분류체계(2018)_대분류_DB유무비율.xlsx"))
#df21<-df2 %>% na.omit %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
#write.xlsx(df21,paste0(dir0,"국가과학기술표준분류체계(2018)_대분류_유형별.xlsx"))

#df3<-df1 %>% select(num,name,mid_cat_2018) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(mid_cat_2018))) 1 else which(!is.na(mid_cat_2018))) %>% distinct() %>% mutate(mid_cat_2018=ifelse(name %in% c('박혜신','양월수','고명국') & mid_cat_2018 %in% c('[OC03]과학기술정책/사회'),'[OC99]기타과학기술과인문사회',mid_cat_2018)) %>% as.data.frame
#df3<-df3 %>% group_by(mid_cat_2018) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
#df30<-df3 %>% mutate(mid1_cat_2018=ifelse(is.na(mid_cat_2018),"DBX","DBO")) %>% group_by(mid1_cat_2018) %>% summarize(total=sum(total)) %>% mutate(percent=round(total/sum(total)*100,2))
#write.xlsx(df30,paste0(dir0,"국가과학기술표준분류체계(2018)_중분류_DB유무비율.xlsx"))
#df31<-df3 %>% na.omit %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
#write.xlsx(df31,paste0(dir0,"국가과학기술표준분류체계(2018)_중분류_유형별.xlsx"))

#df4<-df1 %>% fill(name,.direction=c("down")) %>% mutate(sub_cat_2018=ifelse(name %in% c('박혜신','양월수','고명국') & sub_cat_2018 %in% c('[OC99]기타과학기술과인문사회'),'[OC9999]달리분류되지않는과학기술과인문사회',sub_cat_2018)) %>% filter(sub_cat_2018!="-") %>% group_by(name) %>% distinct(sub_cat_2018) %>% as.data.frame
#df4<-df4 %>% group_by(sub_cat_2018) %>% summarize(total=n()) %>% mutate(percent=round(total/sum(total)*100,2)) %>% arrange(sub_cat_2018) %>% as.data.frame #%>% arrange(desc(percent)) %>% as.data.frame
#df4

##########################
#
#df1<-df %>% mutate(main_cat_2017=str_replace_all(main_cat_2017,'[\r\n]| ',''),
#                   main_cat_2017=str_replace_all(main_cat_2017,'HLG]해안/항만물류','HLG]해안/항만'),
#                   main_cat_2017=str_replace_all(main_cat_2017,'HLG]해안/항만','HLG]해안/항만물류'),
#
#                   mid_cat_2017=str_replace_all(mid_cat_2017,'[\r\n]| ',''),
#                   mid_cat_2017=str_replace_all(mid_cat_2017,'MBT02]해양생명현상규명','MBT02]해양수산생명현상규명'),
#                   mid_cat_2017=str_replace_all(mid_cat_2017,'MEG02]공학','MEG02]선박공학'),
#                   mid_cat_2017=str_replace_all(mid_cat_2017,'MEV03]해양환경위해성평가관리','MEV03]해양환경위해성평가․관리'),
#                   mid_cat_2017=str_replace_all(mid_cat_2017,'OOF0101]해양관측및감시','OOF01]해양관측및감시'),
#
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'[\r\n]| ',''),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'FSP0202]수산저장/포장기술','FSP0202]수산물저장/포장기술'),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'HLG0301]해안항만구조물설계기술','HLG0301]해안·항만구조물설계기술'),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'MBT0303]환경․에너지소재개발기술','MBT0303]환경.에너지소재개발기술'),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'MDP0103]해양이상현상예측및대비기술','MDP0103]해양이상현상예측및대책기술'),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'MRS0101]쇄실성광물자원개발기술','MRS0101]쇄설성광물자원개발기술'),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'MSI0101]연구개발정보수집/보관/가공/분석/유통/활요을위한연구와응용기술','MSI0101]연구개발정보수집/보관/가공/분석/유통/활용을위한연구와응용기술'),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'MEV0303]달리분류되지않는해양환경위해성평가․관리기술','달리분류되지않는해양환경위해성평가·관리기술'),
#                   sub_cat_2017=str_replace_all(sub_cat_2017,'달리분류되지않는해양환경위해성평가·관리기술','[MEV0303]달리분류되지않는해양환경위해성평가·관리기술'),
#
#                   etc_cat_2017=str_replace_all(etc_cat_2017,'[\r\n]| ','')                   
#)

#df2<-df1 %>% select(num,name,main_cat_2017) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(main_cat_2017))) 1 else which(!is.na(main_cat_2017))) %>% distinct() %>% as.data.frame
#df2<-df2 %>% group_by(main_cat_2017) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
#df20<-df2 %>% mutate(main1_cat_2017=ifelse(is.na(main_cat_2017),"DBX","DBO")) %>% group_by(main1_cat_2017) %>% summarize(total=sum(total)) %>% mutate(percent=round(total/sum(total)*100,2))
#write.xlsx(df20,paste0(dir0,"해양수산과학기술분류체계(2017)_대분류_DB유무비율.xlsx"))
#df21<-df2 %>% na.omit %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
#write.xlsx(df21,paste0(dir0,"해양수산과학기술분류체계(2017)_대분류_유형별.xlsx"))

###

#df3<-df1 %>% select(num,name,mid_cat_2017) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(mid_cat_2017))) 1 else which(!is.na(mid_cat_2017))) %>% distinct() %>% as.data.frame
#df3<-df3 %>% group_by(mid_cat_2017) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
#df30<-df3 %>% mutate(mid1_cat_2017=ifelse(is.na(mid_cat_2017),"DBX","DBO")) %>% group_by(mid1_cat_2017) %>% summarize(total=sum(total)) %>% mutate(percent=round(total/sum(total)*100,2))
#write.xlsx(df30,paste0(dir0,"해양수산과학기술분류체계(2017)_중분류_DB유무비율.xlsx"))
#df31<-df3 %>% na.omit %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
#write.xlsx(df31,paste0(dir0,"해양수산과학기술분류체계(2017)_중분류_유형별.xlsx"))

#df4<-df1 %>% select(num,name,sub_cat_2017) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(sub_cat_2017))) 1 else which(!is.na(sub_cat_2017))) %>% distinct() %>% as.data.frame
#df4<-df4 %>% group_by(sub_cat_2017) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
#df40<-df4 %>% mutate(sub1_cat_2017=ifelse(is.na(sub_cat_2017),"DBX","DBO")) %>% group_by(sub1_cat_2017) %>% summarize(total=sum(total)) %>% mutate(percent=round(total/sum(total)*100,2))
#write.xlsx(df40,paste0(dir0,"해양수산과학기술분류체계(2017)_소분류_DB유무비율.xlsx"))
#df41<-df4 %>% na.omit %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
#write.xlsx(df41,paste0(dir0,"해양수산과학기술분류체계(2017)_소분류_유형별.xlsx"))

#df5<-df1 %>% select(num,name,etc_cat_2017) %>% fill(c(num,name),.direction=c("down")) %>% group_by(num,name) %>% slice(if(all(is.na(etc_cat_2017))) 1 else which(!is.na(etc_cat_2017))) %>% distinct() %>% mutate(etc_cat_2017=ifelse(etc_cat_2017=="기계지능화/로봇","기계지능화·로봇",ifelse(etc_cat_2017=="수질정화복원기술","수질정화·복원기술",etc_cat_2017) ))%>%as.data.frame
#df5<-df5 %>% group_by(etc_cat_2017) %>% summarize(total=n()) %>% arrange(desc(total)) %>% as.data.frame
#df50<-df5 %>% mutate(etc1_cat_2017=ifelse(is.na(etc_cat_2017),"DBX","DBO")) %>% group_by(etc1_cat_2017) %>% summarize(total=sum(total)) %>% mutate(percent=round(total/sum(total)*100,2))
#write.xlsx(df50,paste0(dir0,"해양수산과학기술분류체계(2017)_기타_DB유무비율.xlsx"))
#df51<-df5 %>% na.omit %>% mutate(percent=round(total/sum(total)*100,2)) %>% as_tibble
#write.xlsx(df51,paste0(dir0,"해양수산과학기술분류체계(2017)_기타_유형별.xlsx"))


#######################################################################

#png(width = 1600, height = 800, filename = "test.png")
##pdf(file="test.pdf",width = 1600, height = 800)
#ggplot(data=df2,aes(x="", y=percent, fill=reorder(main_cat_2017,total))) +
#  geom_bar(stat="identity", width=1, color="white") +
#  coord_polar("y", start=0) +
#  theme_void() + 
#  geom_text(aes(label=paste0(main_cat_2017,'\n',round(percent,1), '%')),
#            position=position_stack(vjust=0.5), size=6) +
#  #+ theme(legend.position = "none")
#  guides(fill=guide_legend(title="main_cat_2017",reverse=TRUE)) 
#dev.off()


#  theme(legend.position="none") #+
#  geom_text(aes(y = ypos, label = group), color = "white", size=6) +
#  scale_fill_brewer(palette="Set1") 
