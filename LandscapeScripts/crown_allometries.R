library(pacman)
p_load(sf,dplyr,ggplot2,segmented, minipack.lm)

path="D://BCI_50ha_timeseries//crownmap//BCI_50ha_2022_2023_crownmap_raw.shp"

#shp to dataframe
crownmap <- st_read(path)
crownmap <- crownmap %>%
    st_drop_geometry()

colnames(crownmap)

crownmap<-crownmap %>%
    filter(crownmap$crown_area>0)%>%
    mutate(crown_area=as.numeric(crown_area),
        dbh=as.numeric(dbh),
        crown=as.numeric(crown))%>%
    filter(crown_area>0)%>%
    filter(dbh>0)%>%
    filter(!is.na(iluminatio))%>%
    filter(iluminatio > 2)%>%
    filter(crown>3)

nrow(crownmap)

windows()
ggplot(crownmap,aes(x=log(crown_area)))+
    geom_histogram(binwidth=0.1)+
    labs(title="Distribution of Crown Area",
         x="Crown Area (m^2)",
         y="Frequency")+
    theme_minimal()


windows()
ggplot(crownmap,aes(x=log(dbh),y=log(crown_area),color=iluminatio))+
    geom_point()+
    labs(title="Log-Log Plot of Crown Area vs DBH colored by Illumination Index",
         x="Log(DBH) (cm)",
         y="Log(Crown Area) (m^2)")

#crowding index, 
#sum of the distances to neighboring wighted by dbh
# circularity
#for all crowns surrounding the focal crown

windows()
ggplot(crownmap,aes(x=dbh,y=crown_area,color=iluminatio))+
    geom_point()+
    labs(title="Crown area vs DBH colored by illumination",
         x="DBH (cm)",
         y="Crown Area (m^2)")
        labs(title="Crown area vs DBH faceted by illumination",
             x="DBH (cm)",
             y="Crown Area (m^2)")+
        theme(legend.position="none")


#fit a linear model as crown_area ~ dbh

crownmap$log_dbh <- log(crownmap$dbh)
crownmap$log_crown_area <- log(crownmap$crown_area)

log_model <- lm(log_crown_area ~ log_dbh, data = crownmap)
summary(log_model)
pred_area <- exp(predict(log_model))

windows()
ggplot(crownmap,aes(x=log_dbh,y=log_crown_area,color=iluminatio))+
    geom_point()+
    geom_line(aes(y=log(pred_area)), color="blue")+
    labs(title="Crown area vs DBH with predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")+
    labs(title="Crown area vs DBH faceted by illumination with predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")


release_power_model <- nlsLM(
  log_crown_area ~ log_alpha +
    (beta0 + beta1 / (1 + exp(-(log_dbh - c) / s))) * log_dbh,
  data = crownmap,
  start = list(
    log_alpha = coef(log_model)[1],
    beta0 = 0.6,
    beta1 = 0.6,
    c = median(crownmap$log_dbh),
    s = 0.5
  ),
  lower = c(-Inf, 0, 0, min(crownmap$log_dbh), 0.05),
  upper = c(Inf, Inf, Inf, max(crownmap$log_dbh), 2)
)

summary(crownmap$log_dbh)
summary(release_power_model)

pred_release_power <- predict(release_power_model)

windows()
ggplot(crownmap,aes(x=log_dbh,y=log_crown_area,color=iluminatio))+
    geom_point()+
    geom_line(aes(y=pred_release_power), color="red")+
    labs(title="Crown area vs DBH with release power predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")+
#    facet_wrap(~iluminatio,scales="free")+
    labs(title="Crown area vs DBH faceted by illumination with release power predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")



p_load(MASS)

mod_robust <- rlm(log_crown_area ~ log_dbh, data = crownmap)
summary(mod_robust)
pred_robust <- predict(mod_robust)

windows()
ggplot(crownmap,aes(x=log_dbh,y=log_crown_area,color=iluminatio))+
    geom_point()+
    geom_line(aes(y=pred_robust), color="green")+
    labs(title="Crown area vs DBH with robust predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")+
#    facet_wrap(~iluminatio,scales="free")+
    labs(title="Crown area vs DBH faceted by illumination with robust predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")


windows()
ggplot(crownmap, aes(x=dbh,y=crown_area,color=iluminatio))+
    geom_point()+
    geom_line(aes(y=exp(pred_release_power)), color="red")+
    labs(title="Crown area vs DBH colored by species",
         x="DBH (cm)",
         y="Crown Area (m^2)")+
        labs(title="Crown area vs DBH faceted by illumination",
             x="DBH (cm)",
             y="Crown Area (m^2)")+
        theme(legend.position="none")


p_load(nlme)

mod_gls <- gls(
  log_crown_area ~ log_dbh,
  data = crownmap,
  weights = varPower(form = ~ log_dbh)
)
summary(mod_gls)
pred_gls <- predict(mod_gls)

windows()
ggplot(crownmap,aes(x=dbh,y=crown_area,color=iluminatio))+
    geom_point()+
    geom_line(aes(y=exp(pred_gls)), color="blue")+
    labs(title="Crown area vs DBH with GLS predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")+
#    facet_wrap(~iluminatio,scales="free")+
    labs(title="Crown area vs DBH faceted by illumination with GLS predictions",
         x="DBH (cm)",
         y="Crown Area (m^2)")



#lets go with the mass model

path_stems="C:\\Users\\vasquezvicente\\repo\\stem_summ.gpkg"

data<-st_read(path_stems) %>%
    filter(plot_id=="BCI 50 ha plot")

data$diam_mm<- as.numeric(data$diam_cm)*10
data$long<- st_coordinates(data)[,1]
data$lat<- st_coordinates(data)[,2]
pred_log_crown <- predict(mod_robust, newdata = data.frame(log_dbh = log(data$diam_mm)))
data$crown_area <- exp(pred_log_crown)
data$radius <- sqrt(data$crown_area / pi)

pred_csv_path <- "crown_area_predictions.csv"
pred_gpkg_path <- "crown_area_predictions.gpkg"

write.csv(st_drop_geometry(data), pred_csv_path, row.names = FALSE)
