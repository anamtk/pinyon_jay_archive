################################################################################
## This script creates annual BA data for pinyon pine from 2010 to 2022, and modifies
# pinyon cone predictions for analysis

## Code by Kyle C. Rodman, Ecological Restoration Institute. 
# 8/25/2026

################################################################################
### Bring in necessary packages
package.list <- c("here", "tidyverse", "sf", "terra", "rnaturalearth",
                  "tidyterra", "patchwork")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

################################################################################
### Bring in data

## Cone predictions
  # Historical
histCones <- rast(here("Data", "ClimateChange", "Cones", "historic_cone.tif"))
  # Futures
csiroCones <- rast(here("Data", "ClimateChange", "Cones", "cone_predictions_csiro_rcp45.tif"))
gfdlCones <- rast(here("Data", "ClimateChange", "Cones", "cone_predictions_gfdl_rcp45.tif"))
ipslCones <- rast(here("Data", "ClimateChange", "Cones", "cone_predictions_ipsl_rcp45.tif"))

## Future suitability for PIED
piedSuit <- rast(here("Data", "ClimateChange", "PIED_Suit", "pinus_edulis_proj_suits.tif"))

## North american political boundaries
pol_bounds <- ne_countries(continent = c("north america", "south america"), returnclass = "sf") %>%
  st_transform(crs(histCones)) %>%
  filter(!geounit %in% c("Canada", "United States of America"))
studyBounds <- ne_states(country = "United States of America", returnclass = "sf") %>%
  st_transform(crs(histCones))

################################################################################
### Prep layers for plotting

## Get cone production predictions for historical, mid-century, and late-century
  # Historical
histMeanCones_proj <- app(histCones[[2:21]], mean, na.rm = TRUE) ## Mean from 1981-2000. Noel did 1970-2000, but this is more consistent with projection windows
  # Mid-century
csiro_mc_MeanCones <- app(csiroCones[[36:55]], mean, na.rm = TRUE)
gfdl_mc_MeanCones <- app(gfdlCones[[36:55]], mean, na.rm = TRUE)
ipsl_mc_MeanCones <- app(ipslCones[[36:55]], mean, na.rm = TRUE)
mcMeanCones_proj <- project(app(rast(list(csiro_mc_MeanCones,gfdl_mc_MeanCones,ipsl_mc_MeanCones)), mean, na.rm = TRUE),
                            histMeanCones_proj, method = "average")
mcMeanCones_final <- mask(mcMeanCones_proj, histMeanCones_proj)
  # End-century
csiro_ec_MeanCones <- app(csiroCones[[75:94]], mean, na.rm = TRUE)
gfdl_ec_MeanCones <- app(gfdlCones[[75:94]], mean, na.rm = TRUE)
ipsl_ec_MeanCones <- app(ipslCones[[75:94]], mean, na.rm = TRUE)
ecMeanCones_proj <- project(app(rast(list(csiro_ec_MeanCones,gfdl_ec_MeanCones,ipsl_ec_MeanCones)), mean, na.rm = TRUE),
                            histMeanCones_proj, method = "average")
ecMeanCones_final <- mask(ecMeanCones_proj, histMeanCones_proj)

## Climate suitability for PIED
  # Get only the layers we care about
names(piedSuit)
piedSuit_sub <- piedSuit[[c(1, 3, 9)]] ## These three bands represent historical, mid-century, and end-of-century suitability for two-needle Pinyon pine
piedSuit_proj <- project(piedSuit_sub, histMeanCones_proj, method = "average")
piedSuit_final <- mask(piedSuit_proj, histMeanCones_proj)

################################################################################
## Make plots

## Historical cones
(a <- ggplot() + geom_spatraster(data = histMeanCones_proj, aes(fill = mean), maxcell = 1000000) + 
   geom_sf(data = studyBounds, fill = NA) +
   geom_sf(data = pol_bounds, fill = "grey98") + 
   coord_sf(
     xlim = c(ext(histMeanCones_proj)[1:2]),
     ylim = c(ext(histMeanCones_proj)[3:4]),
     expand = FALSE
   ) + theme_bw() + 
   scale_fill_viridis_c(na.value = NA, option = "viridis", direction = 1, 
                        limits = c(0.25,0.45), labels = c(0.25, 0.35, 0.45),
                        breaks = c(0.25, 0.35, 0.45), name = "Mean Cone\nAvailability") +
   theme_bw() + ggtitle("Historical") +  
   theme(plot.title = element_text(hjust = 0.5)))

## Mid-Century cones
(b <- ggplot() + geom_spatraster(data = mcMeanCones_final, aes(fill = mean), maxcell = 1000000) + 
    geom_sf(data = studyBounds, fill = NA) +
    geom_sf(data = pol_bounds, fill = "grey98") + 
    coord_sf(
      xlim = c(ext(histMeanCones_proj)[1:2]),
      ylim = c(ext(histMeanCones_proj)[3:4]),
      expand = FALSE
    ) + theme_bw() + 
    scale_fill_viridis_c(na.value = NA, option = "viridis", direction = 1, 
                         limits = c(0.25,0.45), labels = c(0.25, 0.35, 0.45),
                         breaks = c(0.25, 0.35, 0.45), name = "Mean Cone\nAvailability") +
    theme_bw() + ggtitle("Mid-Century") + 
    theme(plot.title = element_text(hjust = 0.5)))

## End-of-Century cones
(c <- ggplot() + geom_spatraster(data = ecMeanCones_final, aes(fill = mean), maxcell = 1000000) + 
    geom_sf(data = studyBounds, fill = NA) +
    geom_sf(data = pol_bounds, fill = "grey98") + 
    coord_sf(
      xlim = c(ext(histMeanCones_proj)[1:2]),
      ylim = c(ext(histMeanCones_proj)[3:4]),
      expand = FALSE
    ) + theme_bw() + 
    scale_fill_viridis_c(na.value = NA, option = "viridis", direction = 1, 
                         limits = c(0.25,0.45), labels = c(0.25, 0.35, 0.45),
                         breaks = c(0.25, 0.35, 0.45), name = "Mean Cone\nAvailability") +
    theme_bw() + ggtitle("End-of-Century") +
    theme(plot.title = element_text(hjust = 0.5)))

## Historical pied suit
(d <- ggplot() + geom_spatraster(data = piedSuit_final[[1]], aes(fill = current_suit), maxcell = 1000000) + 
    geom_sf(data = studyBounds, fill = NA) +
    geom_sf(data = pol_bounds, fill = "grey98") + 
    coord_sf(
      xlim = c(ext(histMeanCones_proj)[1:2]),
      ylim = c(ext(histMeanCones_proj)[3:4]),
      expand = FALSE
    ) + theme_bw() + 
    scale_fill_viridis_c(na.value = NA, option = "viridis", direction = 1, 
                         limits = c(0,1), labels = c(0.0, 0.5, 1.0),
                         breaks = c(0, 0.5, 1), name = "Climate\nSuitability") +
    theme_bw() + ggtitle("Historical") +  
    theme(plot.title = element_text(hjust = 0.5)))

## Mid-Century pinyon suit
(e <- ggplot() + geom_spatraster(data = piedSuit_final[[2]], aes(fill = median_suit_245_mid), maxcell = 1000000) + 
    geom_sf(data = studyBounds, fill = NA) +
    geom_sf(data = pol_bounds, fill = "grey98") + 
    coord_sf(
      xlim = c(ext(histMeanCones_proj)[1:2]),
      ylim = c(ext(histMeanCones_proj)[3:4]),
      expand = FALSE
    ) + theme_bw() + 
    scale_fill_viridis_c(na.value = NA, option = "viridis", direction = 1, 
                         limits = c(0,1), labels = c(0.0, 0.5, 1.0),
                         breaks = c(0, 0.5, 1), name = "Climate\nSuitability") +
    theme_bw() + ggtitle("Mid-Century") + 
    theme(plot.title = element_text(hjust = 0.5)))

## End-of-Century pinyon suit
(f <- ggplot() + geom_spatraster(data = piedSuit_final[[3]], aes(fill = median_suit_245_end), maxcell = 1000000) + 
    geom_sf(data = studyBounds, fill = NA) +
    geom_sf(data = pol_bounds, fill = "grey98") + 
    coord_sf(
      xlim = c(ext(histMeanCones_proj)[1:2]),
      ylim = c(ext(histMeanCones_proj)[3:4]),
      expand = FALSE
    ) + theme_bw() + 
    scale_fill_viridis_c(na.value = NA, option = "viridis", direction = 1, 
                         limits = c(0,1), labels = c(0.0, 0.5, 1.0),
                         breaks = c(0, 0.5, 1), name = "Climate\nSuitability") +
    theme_bw() + ggtitle("End-of-Century") +
    theme(plot.title = element_text(hjust = 0.5)))

## Boxplots of values - cones
cones_combined <- c(histMeanCones_proj, mcMeanCones_final, ecMeanCones_final)
names(cones_combined) <- c("Historical", "Mid-Century", "End-of-Century")
coneDF <- as.data.frame(cones_combined, na.rm = TRUE) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "TimePeriod", 
    values_to = "ConeAvailability"
  )
coneDF$TimePeriod <- factor(coneDF$TimePeriod, levels = c("Historical", "Mid-Century", "End-of-Century"))
(g <- ggplot(coneDF, aes(x = TimePeriod, y = ConeAvailability, fill = TimePeriod)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    ylim(c(0.25,0.45)) +
    labs(
      x = "Time Period",
      y = "Mean Cone Availability"
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "BrBG", direction = -1))

## Boxplots of values - suitability
suitCombined <- piedSuit_final
names(suitCombined) <- c("Historical", "Mid-Century", "End-of-Century")
suitsDF <- as.data.frame(suitCombined, na.rm = TRUE) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "TimePeriod", 
    values_to = "ClimateSuitability"
  )
suitsDF$TimePeriod <- factor(suitsDF$TimePeriod, levels = c("Historical", "Mid-Century", "End-of-Century"))
(h <- ggplot(suitsDF, aes(x = TimePeriod, y = ClimateSuitability, fill = TimePeriod)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    ylim(c(0,1)) +
    labs(
      x = "Time Period",
      y = "Climate Suitability"
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "BrBG", direction = -1))

### Merge and Export in two draft figure files
(a | b | c)/(d | e | f)/(g | h) +  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'a')
ggsave(here("Figures", "FutureChange.pdf"), width = 11, height = 10,
       device = cairo_pdf())
