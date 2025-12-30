#===============================================================================
#===============Pituffik Stream Isotope Script==================================
#===============================================================================
# This script is focused on the stream water isotopes in the Pituffik and 
# Sioraq Rivers of the Pituffik region, Greenland. These analyses are based
# on samples taken in 2018 and 2019. This follows the broader surface water and
# lake water isotopic analyses completed in a prior project.
# Written Pete D Akers, July 2024-December 2025

# Loading required packages
library(Rmisc)
library(cowplot)
library(tidyverse)
library(simmr)
library(beepr)
library(ggnewscale)

####======Reading in data========####
stream_group_index <- c("Pituffik", "Sioraq")

minor_stream_group_index <- c("Pituffik Mouth", "Pituffik Snoutwash", "Pituffik Ice Wall", "Pituffik Ice Wall Head",
                              "Sioraq Mouth", "Sioraq Pingorsuit", "Other Sioraq", "Sioraq Blackmoss", "Sioraq Tuto", "Sioraq Tuto Head")

freq_stream_group_index <- c("Pituffik Mouth", "Pituffik Snoutwash", "Pituffik Ice Wall", # The sites frequently sampled
                             "Sioraq Mouth", "Sioraq Tuto", "Sioraq Pingorsuit")

stream_iso_full <- as_tibble(read.csv("pituffik_stream_iso_2018_2019.csv", header=TRUE, sep=","))
stream_iso_full$date <- as.Date(stream_iso_full$date, format= "%d/%m/%Y")
stream_iso_full$yr <- year(stream_iso_full$date) # Extracting year
stream_iso_full$doy <- yday(stream_iso_full$date)
stream_iso_full$stream_group <- factor(stream_iso_full$stream_group, levels=stream_group_index)  # Making so the order is always correct
stream_iso_full$minor_stream_group <- factor(stream_iso_full$minor_stream_group, levels=minor_stream_group_index)  # Making so the order is always correct
stream_iso <- stream_iso_full[stream_iso_full$qc_flag != 1, ] # Removing data flagged for quality

# Subsetting samples only from frequently sampled sites and streams
freq_stream_iso <- stream_iso %>%
  filter(minor_stream_group %in% freq_stream_group_index)

# Setting plotting and analysis variables
iso_clm <- c("d18O", "d2H", "dxs") # columns with plotting data
iso_clm_index <- NA
for (i in 1:length(iso_clm)) {
  iso_clm_index[i] <- which(colnames(stream_iso) == iso_clm[i])
}

# Creating sourcetype variable sets
sourcetype_index <- c("glacial", "snowpack", "prcp_act", "lentic")
sourcetype_detailed_index <- c("direct_ice", "firn", "runoff", "premelt", "freshet", 
                               "prcp", "active", "lake", "pool")

# Stream water source isotopic data
sources_iso_full <- as_tibble(read.csv("pituffik_stream_sources_iso_2018_2019.csv", header=TRUE, sep=","))
sources_iso_full$date <- as.Date(sources_iso_full$date, format= "%d/%m/%Y")
sources_iso_full$yr <- year(sources_iso_full$date)
sources_iso_full$doy <- yday(sources_iso_full$date)
sources_iso_full$mo <- month(sources_iso_full$date)

# Further source data structuring and QC
sources_iso_full$sourcetype <- factor(sources_iso_full$sourcetype, levels=sourcetype_index)  # Making so the order is always correct
sources_iso_full$sourcetype_detailed <- factor(sources_iso_full$sourcetype_detailed, levels=sourcetype_detailed_index)  # Making so the order is always correct
sources_iso <- sources_iso_full %>%
  filter(qc_flag != 1) # Removing data flagged for quality

# Setting sources summary data
sources_iso_mean <- sources_iso %>% # Stream source waters mean isotopic values
  group_by(sourcetype) %>%
  summarize(across(all_of(iso_clm),  \(x) mean(x, na.rm = TRUE)))
sources_iso_sd <- sources_iso %>% # Stream source waters standard deviation of isotopic values
  group_by(sourcetype) %>%
  summarize(across(all_of(iso_clm),  \(x) sd(x, na.rm = TRUE)))

# Creating color-iso variable table
iso_color_names <- c("turquoise3", "firebrick", "mediumorchid4")
iso_color_index <- data.frame(colnames(stream_iso[iso_clm_index]), iso_color_names[seq(1,length(iso_clm_index))])
colnames(iso_color_index) <- c("iso", "color")

# Other plotting color schemes
colorset_catchment <- setNames(c("deepskyblue3", "gold2"), nm = stream_group_index)
type_labels <- stream_group_index

colorset_sourcetype <- setNames(c("olivedrab", "aquamarine2", "slateblue4", "#FF99FF"), nm = sourcetype_index)


#==Climate and weather data
# Pituffik weather data from Giovanni Muscari and USAF and water vapor isotopes from Akers 2020
pfk_isowx_day <- as_tibble(read.csv("pituffik_iso_wx_day.csv", header=TRUE))
pfk_isowx_day$daybreak <- as.Date(pfk_isowx_day$daybreak, format= "%d/%m/%Y")
pfk_isowx_day <- pfk_isowx_day %>% rename(date = daybreak)
pfk_isowx_day$doy <- yday(pfk_isowx_day$date)
pfk_isowx_day$yr <- year(pfk_isowx_day$date)

# Loading PET data extracted for Pituffik from Singer et al., 2021
# Data extraction point coordinates (76.5 N 68.5 W), bilinear interpolation 
pfk_pet_1981_2023 <- as_tibble(read.csv("pituffik_pet_1981_2023.csv", header=TRUE)) # PET in mm/day
pfk_pet_1981_2023$date <- as.Date(pfk_pet_1981_2023$date, format="%Y-%m-%d")
pfk_isowx_day <- pfk_isowx_day %>%
  left_join(pfk_pet_1981_2023 %>% select(pet, date), by = "date")

# Loading and joining snow height and ice surface data from GEUS THU-L station
thu_L_melt <- as_tibble(read.csv("THU_L_melt_data.csv", header=TRUE)) %>% # Snow height in m
  mutate(date = as.Date(date, format= "%d/%m/%Y"))
pfk_isowx_day <- pfk_isowx_day %>%
  left_join(thu_L_melt %>% select(snow_height, z_ice_surf_yr_anom, date), by = "date")

####========END DATA READ======####


####============================####
####========BASIC STATS=========####
####============================####

# Total stream isotope means and CI
stream_iso %>%
  reframe(across(all_of(iso_clm), CI))

# Total isotope 1*SD
stream_iso %>%
  summarize(across(all_of(iso_clm), sd))

# Total isotope 2*SD
stream_iso %>%
  summarize(across(all_of(iso_clm), sd))*2

####========END BASIC STATS======####


####============================####
####========LWL CALCUL.=========####
####============================####

# Stream LWL by stream group for all sampling sites
lwl_stream_bystreamgroup <- stream_iso %>% # Running regression by type
  group_by(stream_group) %>%
  nest() %>%
  mutate(lwl = map(data, ~lm(d2H~d18O, data=.)))

lwl_slope <- 
  lwl_stream_bystreamgroup %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "d18O") %>%
  dplyr::select(stream_group, estimate, std.error)
colnames(lwl_slope) <- c("stream_group", "slope", "se_slope")
lwl_slope$conf_int_slope <- lwl_slope$se_slope*qnorm(0.975)

lwl_intercept <- 
  lwl_stream_bystreamgroup %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(stream_group, estimate, std.error)
colnames(lwl_intercept) <- c("stream_group", "intercept", "se_intercept")
lwl_intercept$conf_int_intercept <- lwl_intercept$se_intercept*qnorm(0.975)

lwl_r2 <- 
  lwl_stream_bystreamgroup %>%
  mutate(glanceit = map(lwl, broom::glance)) %>%
  unnest(glanceit) %>%
  dplyr::select(stream_group, r.squared, p.value, nobs)

lwl_stream_bystreamgroup_params <- lwl_slope %>% # Combining all data
  inner_join(lwl_intercept, by="stream_group") %>%
  inner_join(lwl_r2, by="stream_group")

print(lwl_stream_bystreamgroup_params)


# LWLs for the frequently sampled stream sites
lwl_freq_stream_bystream <- freq_stream_iso %>% # Running regression by type
  group_by(minor_stream_group) %>%
  nest() %>%
  mutate(lwl = map(data, ~lm(d2H~d18O, data=.)))

lwl_freq_slope <- 
  lwl_freq_stream_bystream %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "d18O") %>%
  dplyr::select(minor_stream_group, estimate, std.error)
colnames(lwl_freq_slope) <- c("minor_stream_group", "slope", "se_slope")
lwl_freq_slope$conf_int_slope <- lwl_freq_slope$se_slope*qnorm(0.975)

lwl_freq_intercept <- 
  lwl_freq_stream_bystream %>%
  mutate(tidyit = map(lwl, broom::tidy)) %>%
  unnest(tidyit) %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(minor_stream_group, estimate, std.error)
colnames(lwl_freq_intercept) <- c("minor_stream_group", "intercept", "se_intercept")
lwl_freq_intercept$conf_int_intercept <- lwl_freq_intercept$se_intercept*qnorm(0.975)

lwl_freq_r2 <- 
  lwl_freq_stream_bystream %>%
  mutate(glanceit = map(lwl, broom::glance)) %>%
  unnest(glanceit) %>%
  dplyr::select(minor_stream_group, r.squared, p.value, nobs)

lwl_freq_stream_bystream_params <- lwl_freq_slope %>% # Combining all data
  inner_join(lwl_freq_intercept, by="minor_stream_group") %>%
  inner_join(lwl_freq_r2, by="minor_stream_group")

print(lwl_freq_stream_bystream_params)

####========END LWL======####


####=============================####
####=====Source Mixing Model=====####
####=============================####

#### NOTE: Because the model produces slightly different output each run, we do
####  not run the model new each time. Instead, we ran 100 iterations of the
####  model (code below) and then calculated the means of those iterations. For
####  the later analyses, we load in and use these mean values. Note also that
####  running the code straight through with 100 iterations takes a long time
####  and is prone to crashing out. We broke up the loop into 20 iterations at
####  a time, run five times, to build up our 100 iterations.

# catchment_iter <- list()
# sample_iter <- list()
# for (i in 1:100) {
#   #### Getting source fractions estimated for whole stream groups for use as priors
#   # Setting up sample input data
#   catchment_sources_simmr_in <- list()
#   catchment_sources_simmr_out <- list()
#   pfk_catchment_sources_mixmodel_round <- list()
#   sourcetype_index_out <- c("glacial_frac", "glacial_frac_sd", "snowpack_frac", "snowpack_frac_sd",
#                                           "prcp_act_frac", "prcp_act_frac_sd", "lentic_frac", "lentic_frac_sd")
# 
#   # Setting up input data
#   catchment_sources_simmr_in <- simmr_load( # Loading in input data
#     mixtures = stream_iso %>% select(all_of(iso_clm)),
#     source_names = sourcetype_index,
#     source_means = sources_iso_mean %>% filter(sourcetype %in% sourcetype_index) %>%
#       mutate(sourcetype = factor(sourcetype, levels = sourcetype_index)) %>% arrange(sourcetype) %>% # Placing this sort to ensure that the sourcetype names are assigned properly in the mixing model
#       select(all_of(iso_clm)),
#     source_sds = sources_iso_sd %>% filter(sourcetype %in% sourcetype_index) %>%
#       mutate(sourcetype = factor(sourcetype, levels = sourcetype_index)) %>% arrange(sourcetype) %>% # Placing this sort to ensure that the sourcetype names are assigned properly in the mixing model
#       select(all_of(iso_clm)),
#     group = stream_iso %>% pull(stream_group)
#   )
# 
#   # Running mixing model
#   catchment_sources_simmr_out <- simmr_mcmc(catchment_sources_simmr_in)
# 
#   # Organizing output data
#   pfk_catchment_sources_mixmodel_iter_full <- list()
#   pfk_catchment_sources_mixmodel_iter <- as_tibble(matrix(ncol=length(sourcetype_index)*2, nrow=length(catchment_sources_simmr_out$output)))
#   for (k in 1:length(catchment_sources_simmr_out$output)){
#     pfk_catchment_sources_mixmodel_iter_full[[k]] <- catchment_sources_simmr_out$output[[k]]$BUGSoutput$summary
#     for (j in 1:length(sourcetype_index)) {
#       pfk_catchment_sources_mixmodel_iter[k,(j*2-1)] <- pfk_catchment_sources_mixmodel_iter_full[[k]][(j+1),1]
#       pfk_catchment_sources_mixmodel_iter[k,(j*2)] <- pfk_catchment_sources_mixmodel_iter_full[[k]][(j+1),2]
#     } # End j
#   } # End k
#   names(pfk_catchment_sources_mixmodel_iter_full) <- names(catchment_sources_simmr_out$output)
#   colnames(pfk_catchment_sources_mixmodel_iter) <- sourcetype_index_out
#   pfk_catchment_sources_mixmodel_iter <- pfk_catchment_sources_mixmodel_iter %>%
#     add_column(stream_group = names(catchment_sources_simmr_out$output), .before = 1)
#   pfk_catchment_sources_mixmodel <- tibble(pfk_catchment_sources_mixmodel_iter%>%select(stream_group), #Rounding to 3 digits for output csv
#                                                         round(pfk_catchment_sources_mixmodel_iter%>%select(-c(stream_group)),3))
#   beep(1)
# 
#   # Setting priors based on stream group fractional means
#   prior_means_catchment <- pfk_catchment_sources_mixmodel %>% select(-sourcetype_index_out[grepl("_sd", sourcetype_index_out)])
#   prior_sd_catchment <- pfk_catchment_sources_mixmodel %>% select(-sourcetype_index_out[!grepl("_sd", sourcetype_index_out)])
#   prior_sd_catchment_backup <- prior_sd_catchment
#   prior_sd_catchment[,-1] <- 0.33 # Using output SD restricts possible variability to near zero. SD=0.33 roughly gives 3SD = 1, providing a full theoretical full range
# 
#   sample_sources_groupprior_simmr_out <- list()
#   pfk_sample_sources_mixmodel_bycatchment <- list()
#   for(j in 1:nrow(prior_means_catchment)) { # Looping over stream groups to include different priors per stream group
#     # Calculating prior input values per simmr functions
#     prior_iter <- simmr_elicit(n_sources=length(sourcetype_index), proportion_means=as.vector(unlist(prior_means_catchment[j,-1])), proportion_sd=as.vector(unlist(prior_sd_catchment[j,-1])))
# 
#     # Running simmr mix model with priors for individual stream samples
#     sample_sources_simmr_in_iter <- simmr_load( # Loading in input data
#       mixtures = stream_iso %>% filter(stream_group == prior_means_catchment$stream_group[j]) %>% select(all_of(iso_clm)),
#       source_names = sourcetype_index,
#       source_means = sources_iso_mean %>% filter(sourcetype %in% sourcetype_index) %>%
#         mutate(sourcetype = factor(sourcetype, levels = sourcetype_index)) %>% arrange(sourcetype) %>% # Placing this sort to ensure that the sourcetype names are assigned properly in the mixing model
#         select(all_of(iso_clm)),
#       source_sds = sources_iso_sd %>% filter(sourcetype %in% sourcetype_index) %>%
#         mutate(sourcetype = factor(sourcetype, levels = sourcetype_index)) %>% arrange(sourcetype) %>% # Placing this sort to ensure that the sourcetype names are assigned properly in the mixing model
#         select(all_of(iso_clm)),
#       group = stream_iso %>% filter(stream_group == prior_means_catchment$stream_group[j]) %>% pull(sample_id)
#     )
# 
#     sample_sources_groupprior_simmr_out[[j]] <- simmr_mcmc(sample_sources_simmr_in_iter, prior_control = list(means=prior_iter$mean, sd=prior_iter$sd,
#                                                                                                               sigma_shape = rep(3, sample_sources_simmr_in_iter$n_tracers),
#                                                                                                               sigma_rate = rep(3/50,sample_sources_simmr_in_iter$n_tracers)))
#     # Organizing output data
#     pfk_sample_sources_mixmodel_iter_full <- list()
#     pfk_sample_sources_mixmodel_iter <- as_tibble(matrix(ncol=length(sourcetype_index)*2, nrow=length(sample_sources_groupprior_simmr_out[[j]]$output)))
#     for (k in 1:length(sample_sources_groupprior_simmr_out[[j]]$output)){
#       pfk_sample_sources_mixmodel_iter_full[[k]] <- sample_sources_groupprior_simmr_out[[j]]$output[[k]]$BUGSoutput$summary
#       for (m in 1:length(sourcetype_index)) {
#         pfk_sample_sources_mixmodel_iter[k,(m*2-1)] <- pfk_sample_sources_mixmodel_iter_full[[k]][(m+1),1]
#         pfk_sample_sources_mixmodel_iter[k,(m*2)] <- pfk_sample_sources_mixmodel_iter_full[[k]][(m+1),2]
#       } # End m
#     } # End k
#     names(pfk_sample_sources_mixmodel_iter_full) <- names(sample_sources_groupprior_simmr_out[[j]]$output)
#     colnames(pfk_sample_sources_mixmodel_iter) <- sourcetype_index_out
#     pfk_sample_sources_mixmodel_iter <- pfk_sample_sources_mixmodel_iter %>%
#       add_column(sample_id = names(sample_sources_groupprior_simmr_out[[j]]$output), .before = 1)
#     pfk_sample_sources_mixmodel_iter <- pfk_sample_sources_mixmodel_iter %>%
#       add_column(stream_group =  left_join(pfk_sample_sources_mixmodel_iter, stream_iso, by="sample_id") %>% pull(stream_group),
#                  minor_stream_group = left_join(pfk_sample_sources_mixmodel_iter, stream_iso, by="sample_id") %>% pull(minor_stream_group),
#                  date = left_join(pfk_sample_sources_mixmodel_iter, stream_iso, by="sample_id") %>% pull(date),
#                  .before = 2)
#     pfk_sample_sources_mixmodel_bycatchment[[j]]<- tibble(pfk_sample_sources_mixmodel_iter%>%select(sample_id, stream_group, minor_stream_group, date), #Rounding to 3 digits for output csv
#                                                      round(pfk_sample_sources_mixmodel_iter%>%select(-c(sample_id, stream_group, minor_stream_group, date)),3))
#   } # End j
#   names(pfk_sample_sources_mixmodel_bycatchment) <- prior_means_catchment$stream_group
#  # print(pfk_sample_sources_mixmodel_bycatchment)
#   pfk_sample_sources_mixmodel <- pfk_sample_sources_mixmodel_bycatchment %>% bind_rows()
# 
#   catchment_iter[[i]] <- pfk_catchment_sources_mixmodel
#   sample_iter[[i]] <- pfk_sample_sources_mixmodel
#   print(i)
#   beep(2)
# } # End i
# beep(3)
# 
# catchment_iter_bound <- catchment_iter %>% bind_rows()
# sample_iter_bound <- sample_iter %>% bind_rows()
# write.csv(catchment_iter_bound, "pituffik_catchment_sources_iterations.csv")
# write.csv(sample_iter_bound, "pituffik_sample_sources_iterations.csv")
# beep(4)
  

# Reading in mixing model output data from the code above
catchment_iterations <- as_tibble(read.csv("pituffik_catchment_sources_iterations.csv", header=TRUE))
sample_iterations <- as_tibble(read.csv("pituffik_sample_sources_iterations.csv", header=TRUE))

# NOTE: Catchment data are where all values in a stream group are combined together
#  for the mixing model run, producing one fractional value per stream catchment.
#  Sample data are where the mixing model ran individually for each stream water sample.

# Calculating overall source fraction means and sds from the 100 iterations of the mixing model.
#  Note that each iteration produces a source fraction mean and an sd of the
#  iteration's distribution of possible fractions. This means that our summarizing
#  produces four columns per source: mean (the mean of the 100 means), sd (the sd of the 100 means),
#  sd_mean (the mean of the 100 sds), and sd_sd (the sd of the 100 sds).
catchment_iterations_summary <- catchment_iterations %>% # Display mean standard deviations by stream group
  group_by(stream_group) %>%
  summarize(across(where(is.numeric), list(
    mean = ~mean(.x, na.rm = TRUE),
    sd = ~sd(.x, na.rm = TRUE)),
    .names = "{.col}_{.fn}"), .groups = "drop")

sample_iterations_summary <- sample_iterations %>% # Display mean standard deviations by stream group
  group_by(sample_id) %>%
  summarize(across(where(is.numeric), list(
    mean = ~mean(.x, na.rm = TRUE),
    sd = ~sd(.x, na.rm = TRUE)),
    .names = "{.col}_{.fn}"), .groups = "drop") %>%
  left_join(stream_iso %>% select(sample_id, stream_group, minor_stream_group, date), by = "sample_id")

# Renaming variables for later analyses and only keeping data columns of means
pfk_catchment_sources_mixmodel <- catchment_iterations_summary %>%
  select(-ends_with("_sd")) %>%
  rename_with(~ sub("_mean$", "", .x), ends_with("_mean"))
pfk_sample_sources_mixmodel <- sample_iterations_summary %>%
  select(-ends_with("_sd")) %>%
  rename_with(~ sub("_mean$", "", .x), ends_with("_mean"))
pfk_sample_sources_mixmodel$date <- as.Date(pfk_sample_sources_mixmodel$date, format="%d/%m/%Y")

####========END MIXING MODELS======####



####============================####
####=====Spatial Analyses ======####
####============================####

# Mean isotopic values by stream group
mean_iso_streams <- stream_iso %>%
  group_by(stream_group) %>%
  summarize(count=n(),mean_d18O=mean(d18O, na.rm=TRUE)*1000, conf_int_d18O=sd(d18O, na.rm=TRUE)/sqrt(sum(!is.na(d18O)))*qnorm(0.975)*1000,
            mean_d2H=mean(d2H, na.rm=TRUE)*1000, conf_int_d2H=sd(d2H, na.rm=TRUE)/sqrt(sum(!is.na(d2H)))*qnorm(0.975)*1000,
            mean_dxs=mean(dxs, na.rm=TRUE)*1000, conf_int_dxs=sd(dxs, na.rm=TRUE)/sqrt(sum(!is.na(dxs)))*qnorm(0.975)*1000)
print(mean_iso_streams)

# Mean isotopic values by stream group for frequently sampled catchments
mean_freq_iso_streams <- freq_stream_iso %>%
  group_by(minor_stream_group) %>%
  summarize(count=n(),mean_d18O=mean(d18O, na.rm=TRUE)*1000, conf_int_d18O=sd(d18O, na.rm=TRUE)/sqrt(sum(!is.na(d18O)))*qnorm(0.975)*1000,
            mean_d2H=mean(d2H, na.rm=TRUE)*1000, conf_int_d2H=sd(d2H, na.rm=TRUE)/sqrt(sum(!is.na(d2H)))*qnorm(0.975)*1000,
            mean_dxs=mean(dxs, na.rm=TRUE)*1000, conf_int_dxs=sd(dxs, na.rm=TRUE)/sqrt(sum(!is.na(dxs)))*qnorm(0.975)*1000)
print(mean_freq_iso_streams)


# Mean source proportions by stream group
pfk_sample_sources_mixmodel %>%
  replace(is.na(.), 0) %>% # This is a safeguard so that excluded sources count as zero
  group_by(stream_group) %>%
  summarize(n=n())
pfk_sample_sources_mixmodel %>% # Display mean contributions by stream group
  replace(is.na(.), 0) %>% # This is a safeguard so that excluded sources count as zero
  group_by(stream_group) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  select(c(stream_group, glacial_frac, snowpack_frac, prcp_act_frac, lentic_frac))
pfk_sample_sources_mixmodel %>% # Display mean standard deviations by stream group
  replace(is.na(.), 0) %>% # This is a safeguard so that excluded sources count as zero
  group_by(stream_group) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  select(c(stream_group, glacial_frac_sd, snowpack_frac_sd, prcp_act_frac_sd, lentic_frac_sd))

# Mean isotopic values by minor stream group
mean_iso_streams_minor <- stream_iso %>% # Calculating the mean isotopic values of all stream groups
  group_by(minor_stream_group) %>%
  summarize(count=n(),mean_d18O=mean(d18O, na.rm=TRUE)*1000, conf_int_d18O=sd(d18O, na.rm=TRUE)/sqrt(sum(!is.na(d18O)))*qnorm(0.975)*1000,
            mean_d2H=mean(d2H, na.rm=TRUE)*1000, conf_int_d2H=sd(d2H, na.rm=TRUE)/sqrt(sum(!is.na(d2H)))*qnorm(0.975)*1000,
            mean_dxs=mean(dxs, na.rm=TRUE)*1000, conf_int_dxs=sd(dxs, na.rm=TRUE)/sqrt(sum(!is.na(dxs)))*qnorm(0.975)*1000)
print(mean_iso_streams_minor)

####========END SPATIAL ANALYSIS======####


#=====================================================================#
####==================Plotting Figures=============================####
#=====================================================================#
# Defining common plot elements
common_x <- scale_x_date(
  name = NULL,
  position = "bottom",
  limits = as.Date(c("2018-06-10", "2018-08-25")),
  date_breaks = "2 weeks",
  date_labels = "%d %b"
)

common_theme <- theme(
  legend.position = "none",
  axis.title.x = element_text(size = 16, color = "gray40"),
  axis.text.x  = element_text(size = 14, color = "gray40", angle = 45, hjust = 1),
  axis.ticks.x = element_line(color = "gray40"),
  axis.line.x  = element_line(color = "gray40")
)

y_axis_theme <- function(y_axis_col) {
  theme(
    axis.title.y = element_text(size = 16, color = y_axis_col),
    axis.text.y  = element_text(size = 14, color = y_axis_col),
    axis.ticks.y = element_line(color = y_axis_col),
    axis.line.y  = element_line(color = y_axis_col)
  )
}

####==================================####
####====Spatial data results plots====####
####==================================####

### Stream sample isotopic values versus source isotopic means
sample_source_iso_d2H_plot <- ggplot() +
  theme_classic() +
  geom_jitter(aes(x=d18O*1000, y=d2H*1000, color=stream_group), data=stream_iso, alpha=1) +
  scale_color_manual(values=colorset_catchment) +
  new_scale_color() +
  geom_point(aes(x=d18O*1000, y=d2H*1000, color=sourcetype), data=sources_iso_mean) +
  geom_errorbar(aes(x=sources_iso_mean$d18O*1000, y=sources_iso_mean$d2H*1000, 
                    ymin=sources_iso_mean$d2H*1000-sources_iso_sd$d2H*1000,
                    ymax=sources_iso_mean$d2H*1000+sources_iso_sd$d2H*1000,
                    color=sources_iso_mean$sourcetype)) +
  geom_errorbar(aes(x=sources_iso_mean$d18O*1000, y=sources_iso_mean$d2H*1000, 
                     xmin=sources_iso_mean$d18O*1000-sources_iso_sd$d18O*1000,
                     xmax=sources_iso_mean$d18O*1000+sources_iso_sd$d18O*1000,
                     color=sources_iso_mean$sourcetype)) +
  geom_text(aes(x=d18O*1000, y=d2H*1000, label=sourcetype, color=sourcetype), data=sources_iso_mean) +
  scale_color_manual(values=colorset_sourcetype) +
  scale_x_continuous(name=paste("d18O","(‰)")) +
  scale_y_continuous(name=paste("d2H","(‰)")) +
  common_theme + y_axis_theme("gray40")

sample_source_iso_dxs_plot <- ggplot() +
  theme_classic() +
  geom_jitter(aes(x=d18O*1000, y=dxs*1000, color=stream_group), data=stream_iso, alpha=1) +
  scale_color_manual(values=colorset_catchment) +
  new_scale_color() +
  geom_point(aes(x=d18O*1000, y=dxs*1000, color=sourcetype), data=sources_iso_mean) +
  geom_errorbar(aes(x=sources_iso_mean$d18O*1000, y=sources_iso_mean$dxs*1000, 
                    ymin=sources_iso_mean$dxs*1000-sources_iso_sd$dxs*1000,
                    ymax=sources_iso_mean$dxs*1000+sources_iso_sd$dxs*1000,
                    color=sources_iso_mean$sourcetype)) +
  geom_errorbar(aes(x=sources_iso_mean$d18O*1000, y=sources_iso_mean$dxs*1000, 
                     xmin=sources_iso_mean$d18O*1000-sources_iso_sd$d18O*1000,
                     xmax=sources_iso_mean$d18O*1000+sources_iso_sd$d18O*1000,
                     color=sources_iso_mean$sourcetype)) +
  geom_text(aes(x=d18O*1000, y=dxs*1000, label=sourcetype, color=sourcetype), data=sources_iso_mean) +
  scale_color_manual(values=colorset_sourcetype) +
  scale_x_continuous(name=paste("d18O","(‰)")) +
  scale_y_continuous(name=paste("dxs","(‰)")) +
  common_theme + y_axis_theme("gray40")

windows(height=8, width=14)
#pdf("Figures/sample_source_compare_plot.pdf", height=8, width=14)
plot_grid(plotlist = list(sample_source_iso_d2H_plot, sample_source_iso_dxs_plot),
          ncol=2, align = "hv", axis = "bt", byrow=FALSE)   
#ggsave("Figures/sample_source_compare_plot.png", height=8, width=14)
#dev.off()

#### END SPATIAL DATA PLOTS ####

####===================================================####
####============THU-L 2018 vs. 2019 plots==============####
####===================================================####

# Setting up tibble to contain sampling period dates for plotting
sample_period <- tibble(yr = c(2018,2019),
                        start = as.Date(c("2018-06-14", "2019-07-17")),
                        end = as.Date(c("2018-08-22", "2019-08-01"))) %>%
  mutate(doy_start = yday(start), doy_end = yday(end), y = c(0.53, 0.55), yend = c(0.53, 0.55))

# THU-L snow height plot
snowht_ts_plot <- ggplot(pfk_isowx_day) +
  theme_classic() +
  geom_line(aes(x = doy, y = snow_height, group=yr, color=factor(yr))) +
  geom_ribbon(aes(x = doy, ymin=0, ymax=snow_height, fill=factor(yr)), alpha=0.15) +
  geom_segment(data=sample_period, aes(x=doy_start, xend=doy_end, y=y, yend=yend, group=yr, color=factor(yr)), inherit.aes = FALSE) +
  scale_color_manual(values=c("2018" = "mediumpurple1", "2019" = "violetred3")) +
  scale_fill_manual(values=c("2018" = "mediumpurple1", "2019" = "violetred3")) +
  scale_x_continuous(name = NULL, position = "bottom", 
                     breaks = cumsum(days_in_month(1:12))-31,labels = month.abb) +
  scale_y_continuous(name = "Snow Height (m)", position = "left", limits = c(0, 0.55)) +
  common_theme + y_axis_theme("gray40")

# Finding plotting values for ice surface loss
ice_surf <- pfk_isowx_day %>%
  filter(yr == 2018 | yr == 2019) %>%
  group_by(yr) %>%
  summarise(ice_surf_anom = min(z_ice_surf_yr_anom, na.rm = TRUE), .groups = "drop")

# THU-L surface ice change plot
surf_ice_change_plot <- ggplot(data=ice_surf, aes(x = factor(yr), y = ice_surf_anom, fill = factor(yr))) +
  theme_classic() +
  geom_col(width = 0.5) +
  geom_hline(yintercept=0, color="gray40", linewidth=1.0) +
  scale_fill_manual(values = c("2018" = "mediumpurple1", "2019" = "violetred3")) +
  labs(x = NULL, y = "Ice Surface Change (m)") +
  common_theme + y_axis_theme("gray40")

windows(height=4, width=14)
#pdf("Figures/thuL_plots.pdf", height=4, width=14)
plot_grid(plotlist = list(snowht_ts_plot, surf_ice_change_plot),
          ncol=2, align = "hv", rel_widths = c(0.8, 0.2))   
#ggsave("Figures/thuL_plots.png", height=4, width=14)
#dev.off()

#### END THU-L PLOTS ####

####===================================================####
####=====Times Series of Lake and Stream Isotopes======####
####===================================================####

#======Stream temporal plots======#
stream_ts_pituffik_index <- c("Pituffik Mouth", "Pituffik Snoutwash", "Pituffik Ice Wall")
stream_ts_sioraq_index <- c("Sioraq Mouth", "Sioraq Tuto", "Sioraq Pingorsuit")
stream_ts_index <- list(stream_ts_pituffik_index, stream_ts_sioraq_index)
stream_ts_label_index <- c("Pituffik River", "Sioraq River")
yr_index <- c(2018, 2019)
stream_ts_plot <- list()
stream_ts_plot_byyear <- list()
title_iter <- list()
ylims <- data.frame(c(-25,-185,3),c(-19,-140,14))
colnames(ylims) <- c("ymin", "ymax")

for (j in seq_along(stream_ts_index)) {
  for (k in seq_along(yr_index)) {
    data_iter <- stream_iso %>%
      filter(yr == yr_index[k],
             minor_stream_group %in% stream_ts_index[[j]]) %>%
      mutate(minor_stream_group = factor(minor_stream_group,
                                         levels = stream_ts_index[[j]]))
    title_iter <- ggdraw() + 
      draw_label(stream_ts_label_index[j], x = 0, hjust = 0, size = 14) +
      theme(plot.margin = margin(0, 0, 0, 4))
    
    # Build isotope plots as a list
    iso_ts_plot_stream <- map(seq_along(iso_clm_index), function(i) {
      
      iso_name <- colnames(data_iter)[iso_clm_index[i]]
      color_select <- iso_color_index %>%
        filter(iso == iso_name) %>%
        pull(2) %>% as.character()
      
      data_iter %>%
        ggplot(aes(
          x = as.Date(doy, origin = "2017-12-31"),
          y = .data[[iso_clm[i]]] * 1000,
          group = factor(minor_stream_group)
        )) +
        theme_classic() +
        geom_line(aes(linetype = factor(minor_stream_group)), color = color_select) +
        geom_point(aes(shape = factor(minor_stream_group)), color = color_select) +
        scale_y_continuous(
          name = paste0(iso_name, " (‰)"),
          limits = c(ylims[i, 1], ylims[i, 2])
        ) +
        common_x + common_theme + y_axis_theme(color_select)
    })
    
    # Combine title and isotope plots
    stream_ts_plot[[(j - 1) * length(yr_index) + k]] <- plot_grid(
      title_iter,
      plotlist = iso_ts_plot_stream,
      ncol = 1, align = "v",
      rel_heights = c(0.2, rep(1, length(iso_ts_plot_stream)))
    )
  }
}

# Source water isotope boxplots (to go at right side for reference)
source_stream_temporal_boxplot <- list()
title_iter <- ggdraw() + 
  draw_label("Sources", x = 0, hjust = 0, size=14) +
  theme(plot.margin = margin(0, 0, 0, 4))

source_stream_temporal_boxplot <- map(seq_along(iso_clm), function(i) {
  
  iso_name <- iso_clm[i]
  
  color_select <- iso_color_index %>%
    filter(iso == colnames(data_iter[iso_clm_index[i]])) %>%
    pull(2) %>% as.character()
  
  ggplot(data = sources_iso, aes(x = sourcetype, y = .data[[iso_name]] * 1000)) +
    theme_classic() +
    geom_boxplot(aes(color = sourcetype, fill = sourcetype)) +
    scale_color_manual(values = colorset_sourcetype) +
    scale_fill_manual(values = alpha(colorset_sourcetype, 0.3)) +
    scale_y_continuous(name = paste0(iso_name, " (‰)")) +
    scale_x_discrete(labels = toupper(sourcetype_index), position = "bottom", name = NULL) +
    coord_cartesian(ylim = c(ylims[i, 1], ylims[i, 2])) +
    common_theme + y_axis_theme("gray40")
})

p_source_boxplots <- plot_grid(title_iter,
                               plotlist=source_stream_temporal_boxplot,
                               ncol=1, align = "v",
                               rel_heights = c(0.2, seq(1,1,length=length(source_stream_temporal_boxplot))))

# Adding boxplot to stream iso time series plots
stream_ts_plot[[length(stream_ts_plot)+1]] <- p_source_boxplots

windows(height=8, width=12)
#pdf("Figures/stream_temporal_plots.pdf", height=8, width=12)
plot_grid(plotlist = stream_ts_plot, ncol=5, align = "hv",
          rel_widths = c(seq(1,1,length=length(stream_ts_label_index)*length(yr_index)), 0.5))
#ggsave("Figures/stream_temporal_plots.png", height=8, width=12, dpi=600)
#dev.off()

#======Stream source mixing model temporal plots======#
stream_ts_pituffik_index <- c("Pituffik Mouth", "Pituffik Snoutwash", "Pituffik Ice Wall")
stream_ts_sioraq_index <- c("Sioraq Mouth", "Sioraq Tuto", "Sioraq Pingorsuit")
stream_ts_index <- list(stream_ts_pituffik_index, stream_ts_sioraq_index)
stream_ts_label_index <- c("Pituffik River", "Sioraq River")
yr_index <- c(2018, 2019)

# Creating color-source mixing model variable table
mixmodel_clm <- c("glacial_frac", "snowpack_frac", "prcp_act_frac", "lentic_frac")
mixmodel_color_index <- tibble(source_frac = mixmodel_clm, color= colorset_sourcetype[seq_along(mixmodel_clm)])

max_y <- pfk_sample_sources_mixmodel %>%
  filter(stream_group %in% c("Pituffik", "Sioraq")) %>%
  summarize(max_y = max(.[mixmodel_clm], na.rm=TRUE)) %>%
  pull(max_y)
source_frac_ylims <- data.frame(ymin = rep(0,length(mixmodel_clm)), ymax = rep(max_y,length(mixmodel_clm)))
prior_means_lineplot <- pfk_catchment_sources_mixmodel %>% select(stream_group, all_of(mixmodel_clm)) 

mixmodel_ts_plot <- list()

for (j in seq_along(stream_ts_index)) {
  for (k in seq_along(yr_index)) {
    data_mixmodel_iter <- pfk_sample_sources_mixmodel %>%
      mutate(yr = year(date), doy = yday(date)) %>%
      filter(yr == yr_index[k],
             minor_stream_group %in% stream_ts_index[[j]]) %>%
      mutate(minor_stream_group = factor(minor_stream_group,
                                         levels = stream_ts_index[[j]]))
    title_iter <- ggdraw() + 
      draw_label(stream_ts_label_index[j], x = 0, hjust = 0, size = 14) +
      theme(plot.margin = margin(0, 0, 0, 4))
    
    # Mixing model plots per column
    mixmodel_ts_plot_stream <- map(seq_along(mixmodel_clm), function(i) {
      
      color_select <- mixmodel_color_index %>%
        filter(source_frac == mixmodel_clm[i]) %>%
        pull(color) %>% as.character()
      
      ggplot(data_mixmodel_iter, aes(
        x = as.Date(doy, origin = "2017-12-31"),
        y = .data[[mixmodel_clm[i]]],
        group = factor(minor_stream_group)
      )) +
        theme_classic() +
        geom_hline(yintercept = prior_means_lineplot[[j, i + 1]],
                   linetype = "dashed", color = "gray60") +
        geom_line(aes(linetype = factor(minor_stream_group)), color = color_select) +
        geom_point(aes(shape = factor(minor_stream_group)), color = color_select) +
        scale_y_continuous(mixmodel_clm[i], position = "left",
                           limits = c(source_frac_ylims[i, 1], source_frac_ylims[i, 2])) +
        common_x + common_theme + y_axis_theme(color_select)
    })
    
    # Weather data selection
    wxdata_select <- pfk_isowx_day %>% filter(yr == yr_index[k])
    
    # Temperature plot
    tavg_ts_plot <- ggplot(wxdata_select) +
      theme_classic() +
      geom_line(aes(x = as.Date(doy, origin = "2017-12-31"), y = tavg), color = "tomato4") +
      scale_y_continuous(name = "Air Temperature (°C)", position = "left", limits = c(-2, 15)) +
      common_x + common_theme + y_axis_theme("tomato4")
    
    # Precipitation plot
    prcp_ts_plot <- ggplot(wxdata_select) +
      theme_classic() +
      geom_col(aes(x = as.Date(doy, origin = "2017-12-31"), y = prcp_H2O), fill = "slateblue4") +
      scale_y_continuous(name = "Precipitation (mm)", position = "left", limits = c(0, 10.5)) +
      common_x + common_theme + y_axis_theme("slateblue4")
    
    # PET plot
    pet_ts_plot <- ggplot(wxdata_select) +
      theme_classic() +
      geom_line(aes(x = as.Date(doy, origin = "2017-12-31"), y = pet), color = "darkorange") +
      scale_y_continuous(name = "PET (mm d-1)", position = "left", limits = c(0, 4)) +
      common_x + common_theme + y_axis_theme("darkorange")
    
    # Combine title, mixing model plots, and weather plots
    mixmodel_ts_plot[[(j - 1) * length(yr_index) + k]] <- plot_grid(
      title_iter,
      plotlist = append(mixmodel_ts_plot_stream, list(tavg_ts_plot, prcp_ts_plot, pet_ts_plot)),
      ncol = 1, align = "v",
      rel_heights = c(0.2, rep(1, length(mixmodel_ts_plot_stream)), 0.5, 0.5, 0.5)
    )
  }
}
windows(height=12.5, width=10)
#pdf("Figures/sample_sources_temporal_plots.pdf", height=12.5, width=10)
plot_grid(plotlist = mixmodel_ts_plot, ncol=4, align = "hv")
#ggsave("Figures/sample_sources_temporal_plots.png", height=12.5, width=10, dpi=600)
#dev.off()

#### END TIME SERIES PLOTS ####

