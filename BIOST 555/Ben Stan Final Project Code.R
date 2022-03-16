# Load libraries and data
library(dplyr)
library(ggplot2)
library(haven)
library(sf)
library(INLA)
library(SUMMER)
library(survey)
library(tinytex)
library(knitr)
library(rgdal)
library(sp)
library(spdep)
library(SpatialEpi)
library(RColorBrewer)
library(ggplot2)
library(ggridges)
rm(list = ls())
options(digits = 3) 
setwd("/Users/bstan/Documents/UW/Courses/BIOST 555")
cols_needed = c("v001","v023","v024","v005","v190","v301","v302a")
# v302a is answer to "Ever used anything or tried to delay or avoid getting pregnant"
individual = read_dta("cameroon/cameroon_dhs/cameroon_dhs.DTA",
                      col_select = cols_needed)
hiv = read_dta("cameroon/cameroon_hiv/cameroon_hiv.DTA")

# Calculate Covariate values
individual$v302a = as.numeric(individual$v302a)
individual_cov = individual %>% filter(!is.na(v302a))
table(individual_cov$v302a)/nrow(individual_cov)
cluster_cov = individual_cov %>% group_by(v001) %>% 
  dplyr::summarize(n_prevent = sum(v302a),
                   n_total = n(),
                   p_naive = sum(v302a)/n(),
                   weighted_numerator = sum(v005*v302a),
                   weighted_denominator = sum(v005),
                   p_weighted = sum(v005*v302a)/sum(v005))

# Join tables to get strata for individuals
indiv_trim = individual %>% distinct(v001,v023,v024) %>% 
  dplyr::rename(hivclust=v001,
                strata=v023,
                region=v024)
full_indiv = hiv %>% inner_join(indiv_trim,by="hivclust")
table(full_indiv$hiv03)

# Load maps and plot admin regions along with cluster locations
regions_dhs = st_read("cameroon/dhs_shapefiles/CMGE71FL.shp")
regions_gadm1 = st_read("cameroon/gadm_shapefiles/gadm40_CMR_1.shp")
regions_gadm2 = st_read("cameroon/gadm_shapefiles/gadm40_CMR_2.shp")
regions_gadm2 = regions_gadm2 %>% arrange(NAME_2)
ggplot(data=regions_gadm1) + geom_sf(aes(fill=VARNAME_1))
ggplot(data=regions_gadm2) + geom_sf(aes(fill=NAME_2)) + 
  theme(legend.text=element_text(size=5),legend.key.size=unit(.45,'cm'))
ggplot(data=regions_gadm2) + geom_sf()
ggplot() + geom_sf(data=regions_gadm2) + 
  geom_sf(data=regions_dhs, alpha=0.3, size=0.8, color="blue")
ggplot(data=regions_dhs) + geom_sf(aes(col=ADM1NAME))

# Join cluster table and admin 2 table
regions_dhs_trim = regions_dhs %>% select(DHSCLUST,ADM1DHS,ADM1NAME,URBAN_RURA,geometry)
regions_gadm2_trim = regions_gadm2 %>% select(NAME_1,NAME_2,geometry)
cluster_geo = st_join(regions_gadm2_trim,regions_dhs_trim,left=TRUE)

# Create smaller version of table and count number of clusters per admin2
cluster_geo_trim = cluster_geo %>% select(NAME_1,NAME_2,DHSCLUST,URBAN_RURA)
cluster_geo_grouped = cluster_geo_trim %>% 
  mutate(cluster_string = as.character(DHSCLUST)) %>% 
  group_by(NAME_2) %>% 
  dplyr::summarize(clusters = paste0(cluster_string, collapse = ","),
                   num_clusters = n())
# Plot map of clusters per region
ggplot(data=cluster_geo_grouped) + 
  geom_sf(aes(fill=num_clusters)) +
  scale_fill_viridis_c(option = "D",direction = -1)
# Plot histogram of clusters per region
ggplot(data=cluster_geo_grouped) + 
  geom_histogram(aes(x=num_clusters,y=..density..)) +
  labs(x = "Number of Clusters in Admin 2 Region", y = "Proportion")
# Join with individual HIV data
st_geometry(cluster_geo_trim) = NULL
#print(nrow(cluster_geo_trim %>% filter(is.na(DHSCLUST)))) # 7 regions without a DHS cluster
df = inner_join(full_indiv,cluster_geo_trim,by=c("hivclust"="DHSCLUST"))
df = df %>% filter(hiv03 != 9)
head(df)

# Join with covariate data to plot covariate values
admin2_cov = cluster_cov %>% left_join(cluster_geo_trim,by=c("v001"="DHSCLUST")) %>%
  group_by(NAME_2) %>% 
  dplyr::summarize(n_prevent = sum(n_prevent),
                   n_total = sum(n_total),
                   p_naive = sum(n_prevent)/sum(n_total),
                   weighted_numerator = sum(weighted_numerator),
                   weighted_denominator = sum(weighted_denominator),
                   p_weighted = sum(weighted_numerator)/sum(weighted_denominator))
admin2_cov_geo = regions_gadm2 %>% left_join(admin2_cov,by="NAME_2")
ggplot(data=admin2_cov_geo) + 
  geom_sf(aes(fill=p_naive)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=admin2_cov_geo) + 
  geom_sf(aes(fill=p_weighted)) +
  scale_fill_viridis_c(option = "D",direction = -1)

# Calculate unweighted estimates
design_unweighted = svydesign(ids = ~hivclust,
                              strata = ~strata,
                              data = df)
results_unweighted = svyby(~hiv03, ~NAME_2, design_unweighted, svymean) %>% 
  dplyr::rename(unweighted_est=hiv03,
                unweighted_se=se) 
# Calculate weighted estimates
design_weighted = svydesign(ids = ~hivclust, 
                            weights = ~hiv05, 
                            strata = ~strata,
                            data = df)
results_weighted = svyby(~hiv03, ~NAME_2, design_weighted, svymean) %>% 
  dplyr::rename(weighted_est=hiv03,
                weighted_se=se) 
# Calculate smoothed estimates
nb.r = poly2nb(regions_gadm2, queen = F, row.names = regions_gadm2$NAME_2) 
mat = nb2mat(nb.r, style = "B", zero.policy = TRUE)
colnames(mat) = sort(unique(regions_gadm2$NAME_2))
rownames(mat) = sort(unique(regions_gadm2$NAME_2))
df$NAME_2 = factor(df$NAME_2, levels = regions_gadm2$NAME_2)
smoothed = smoothSurvey(data = as.data.frame(df), geo = regions_gadm2, Amat = mat, 
                        responseType = "binary", responseVar = "hiv03", 
                        strataVar = NULL, weightVar = NULL, regionVar = "NAME_2", 
                        clusterVar = NULL, CI = 0.95)
smooth_post = smoothed$smooth %>% 
  mutate(smoothed_sd=sqrt(var)) %>% 
  dplyr::rename(smoothed_median=median)
smooth_trim = smooth_post %>% dplyr::select(region,smoothed_median,smoothed_sd)
# Calculate FH estimates
smoothed_weight = smoothSurvey(data = as.data.frame(df), geo = regions_gadm2, Amat = mat, 
                               responseType = "binary", responseVar = "hiv03", 
                               strataVar = "strata", weightVar = "hiv05", regionVar = "NAME_2", 
                               clusterVar = "~hivclust", CI = 0.95)
smooth_weight_post = smoothed_weight$smooth %>% 
  mutate(fh_sd=sqrt(var)) %>% 
  dplyr::rename(fh_median=median)
smooth_weight_trim = smooth_weight_post %>% dplyr::select(region,fh_median,fh_sd)
# Calculate FH estimates with covariate
all_admin2 = as.data.frame(rownames(mat)) %>% dplyr::rename("NAME_2" = "rownames(mat)")
cov_df = all_admin2 %>% left_join(admin2_cov %>% dplyr::select(NAME_2,p_weighted), by="NAME_2")
smoothed_weight_cov = smoothSurvey(data = as.data.frame(df), geo = regions_gadm2, Amat = mat,
                                   X = as.data.frame(cov_df), responseType = "binary", 
                                   responseVar = "hiv03",strataVar = "strata", weightVar = "hiv05", 
                                   regionVar = "NAME_2", clusterVar = "~hivclust", CI = 0.95)
smooth_weight_cov_post = smoothed_weight_cov$smooth %>% 
  mutate(fh_cov_sd=sqrt(var)) %>% 
  dplyr::rename(fh_cov_median=median)
smooth_weight_cov_trim = smooth_weight_cov_post %>% 
  dplyr::select(region,fh_cov_median,fh_cov_sd)
# Compare phi for FH models with and without covariate
smoothed_weight$fit$summary.hyperpar
smoothed_weight_cov$fit$summary.hyperpar

# Create unified table with all estimates and plot
all_est = df %>% 
  group_by(NAME_2) %>% 
  dplyr::summarize(count = n(),
                   prop_hiv_naive = sum(hiv03)*1.0/n()) %>% 
  left_join(results_unweighted,by="NAME_2") %>% 
  left_join(results_weighted,by="NAME_2") %>% 
  left_join(smooth_trim,by=c("NAME_2"="region")) %>% 
  left_join(smooth_weight_trim,by=c("NAME_2"="region")) %>% 
  left_join(smooth_weight_cov_trim,by=c("NAME_2"="region"))
head(all_est)
full_results = regions_gadm2_trim %>% 
  left_join(all_est,by="NAME_2") %>% 
  left_join(admin2_cov,by="NAME_2")
ggplot(data=full_results) + geom_sf(aes(fill=unweighted_est)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=weighted_est)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=smoothed_median)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=fh_median)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=fh_cov_median)) +
  scale_fill_viridis_c(option = "D",direction = -1)

# Plot standard errors/deviations
ggplot(data=full_results) + geom_sf(aes(fill=unweighted_se)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=weighted_se)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=smoothed_sd)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=fh_sd)) +
  scale_fill_viridis_c(option = "D",direction = -1)
ggplot(data=full_results) + geom_sf(aes(fill=fh_cov_sd)) +
  scale_fill_viridis_c(option = "D",direction = -1)

# Compare various estimates to one another
ggplot(full_results, aes(x = unweighted_est, y = weighted_est)) + 
  geom_point() +
  labs(x = "Simple (Unweighted) Estimate", y = "Direct (Weighted) Estimate") +
  geom_abline(intercept = 0,slope = 1, color = "red")
ggplot(full_results, aes(x = unweighted_est, y = smoothed_median)) + 
  geom_point() +
  labs(x = "Simple (Unweighted) Estimate", y = "Smoothed Estimate") +
  geom_abline(intercept = 0,slope = 1, color = "red") +
  xlim(0,0.1) + 
  ylim(0,0.1)
ggplot(full_results, aes(x = weighted_est, y = fh_median)) + 
  geom_point() +
  labs(x = "Direct (Weighted) Estimate", y = "FH Estimate") +
  geom_abline(intercept = 0,slope = 1, color = "red") +
  xlim(0,0.15) + 
  ylim(0,0.15)
ggplot(full_results, aes(x = weighted_se, y = fh_sd)) + 
  geom_point() +
  labs(x = "Direct (Weighted) SE", y = "FH Posterior SD") +
  geom_abline(intercept = 0,slope = 1, color = "red") +
  xlim(0,0.05) + 
  ylim(0,0.05)
ggplot(full_results, aes(x = p_weighted, y = weighted_est)) + 
  geom_point() +
  geom_smooth() +
  labs(x = "Proportion with Pregnancy Prevention", y = "Direct (Weighted) Estimate")
ggplot(full_results, aes(x = fh_median, y = fh_cov_median)) + 
  geom_point() +
  geom_abline(intercept = 0,slope = 1, color = "red") +
  labs(x = "FH Estimate", y = "FH Estimate with Covariate")