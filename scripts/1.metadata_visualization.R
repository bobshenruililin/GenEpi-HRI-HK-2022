# Ensure required packages are installed
options(repos = c(CRAN = "https://cloud.r-project.org"))
required_packages <- c(
  "tidyverse", "jsonlite", "dplyr", "tidyr", "lubridate", "ggplot2",
  "ggforce", "ggrepel", "MetBrewer", "ggbump", "scales", "sp",
  "gridExtra", "RColorBrewer", "writexl", "readxl"
)
to_install <- setdiff(required_packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

# Load libraries
library(tidyverse)
library(jsonlite)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(MetBrewer)
library(ggbump)
library(scales) # For the rescale function
library(sp)
library(gridExtra)
library(RColorBrewer)
library(grid)      # for unit()
# readxl and writexl are used via :: but ensure installed above

# Panel A: Aggregated confirmed case classification and sequenced samples in Hong Kong by date (2022).
data <- fromJSON("data/case_curve/HK_case_data.json")
## Note: The JSON fields for the starting date are zero-indexed for the month, so add 1.
start_date <- as.Date(sprintf("%d-%02d-%02d",
                                data$hkcaseaggr_v202209_confirm_first_date_utc_y,
                                data$hkcaseaggr_v202209_confirm_first_date_utc_m + 1,
                                data$hkcaseaggr_v202209_confirm_first_date_utc_d))

## Create a vector of dates for the time series
n <- length(data$hkcaseaggr_v202209_confirm_classification_imported_pcr)
dates <- start_date + 0:(n - 1)

## Build a data frame with the four series.
df <- tibble(
  date = dates,
  `PCR (imported)` = data$hkcaseaggr_v202209_confirm_classification_imported_pcr,
  `PCR (local)`    = data$hkcaseaggr_v202209_confirm_classification_local_pcr,
  `RAT (imported)` = data$hkcaseaggr_v202209_confirm_classification_imported_rat,
  `RAT (local)`    = data$hkcaseaggr_v202209_confirm_classification_local_rat
)
writexl::write_xlsx(df, "data/case_curve/hk_epicurve.xlsx")

## Limit the data to 2022
df_2022 <- df %>%
  filter(date >= as.Date("2022-01-01") & date <= as.Date("2022-12-31"))

## Convert from wide to long format.
df_long <- df_2022 %>%
  pivot_longer(cols = c(`PCR (imported)`, `PCR (local)`, `RAT (imported)`, `RAT (local)`),
               names_to = "case_type", values_to = "count")

## Updated helper function to navigate the complex grob structure
ggforce_hack_facet_zoom_x <- function(p) {
  # Instead of trying to modify the grob structure, let's use a direct approach
  # Create two versions of the plot with different text angles
  p_main_rotated <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p_all_horizontal <- p + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  # Convert to grobs
  g_rotated <- ggplotGrob(p_main_rotated)
  g_horizontal <- ggplotGrob(p_all_horizontal)
  
  # Identify axis components
  main_axis_id <- which(g_rotated$layout$name == "axis-b-1")
  zoom_axis_id <- which(g_rotated$layout$name == "axis-b-4")
  
  # Take the rotated main axis and horizontal zoom axis
  if (length(main_axis_id) > 0 && length(zoom_axis_id) > 0) {
    g_horizontal$grobs[[main_axis_id]] <- g_rotated$grobs[[main_axis_id]]
    return(g_horizontal)
  } else {
    print("Could not identify axis elements, returning original plot")
    # Try a simpler direct theme approach
    p_direct <- p + 
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      theme(strip.placement = "outside") +
      coord_cartesian(clip = "off")
    return(ggplotGrob(p_direct))
  }
}

## Panel A2: GISAID sequences (removing known hospital-acquired infections), representing community genomic surveillance

metadata_gisaid <- readxl::read_xlsx("data/GISAID_data/filter.xlsx")
### only keep the lineages achieving >= 2% frequency, relabel the rest as "Other"
check <- table(metadata_gisaid$pangolin_lineage)/nrow(metadata_gisaid)>0.02
vec_lineages_keep <- names(check)[check]
metadata_gisaid$Lineages <- ifelse(metadata_gisaid$pangolin_lineage %in% vec_lineages_keep, metadata_gisaid$pangolin_lineage, "Others")
df_lin_freq <- metadata_gisaid %>% transmute(date = as_date(date), Lineages) %>% filter(!is.na(date)) %>% group_by(date, Lineages) %>% summarise(`No. of sequences` = n())

## Simplified approach for the final figure
fig_1a2 <- ggplot() +
  geom_col(aes(x = date, y = count, fill = case_type), data=df_long) +
  geom_label(x=ymd("2022-04-03"), y=60000, label="HK Fifth wave (phase 1)", size=2.5, color="black", alpha=0.8, data=df_long %>% filter(date=="2022-03-14" & case_type=="PCR (local)")) +
  geom_label(x=ymd("2022-08-01"), y=20000, label="HK Fifth wave (phase 2)", size=2.5, color="black", alpha=0.8, data=df_long %>% filter(date=="2022-08-01" & case_type=="PCR (local)")) +
  geom_vline(xintercept = as.Date("2022-05-06"), linetype = "dashed", color = "#053C29", size = 0.8, alpha = 0.8) + # Relaxed visit scheduel first round
  geom_vline(xintercept = as.Date("2022-05-31"), linetype = "dashed", color = "#053C29", size = 0.8, alpha = 0.8) + # Relaxed visit scheduel second round
  geom_label(x=ymd("2022-05-31"), y=60000, label="Relaxed visit schedule\n(2nd round, remaining major hospitals)", size=2.5, color="#053C29", alpha=0.8, data=df_long %>% filter(date=="2022-05-31" & case_type=="PCR (local)")) +
  geom_label(x=ymd("2022-05-06"), y=60000, label="Relaxed visit schedule\n(1st round, 26 public hospitals)", size=2.5, color="#053C29", alpha=0.8, data=df_long %>% filter(date=="2022-05-06" & case_type=="PCR (local)")) +
  scale_x_date(date_labels = "%b-%d", date_breaks = "15 days", expand = c(0, 0)) +
  labs(
    x = "Date of confirmation",
    y = "No. of cases",
    fill = "Case Type"
  ) +
  scale_fill_manual(values = c(`PCR (imported)` = "blue", `PCR (local)` = "darkblue",
                               `RAT (imported)` = "red", `RAT (local)` = "darkred")) +
  theme(
    legend.position = "bottom", 
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey", size = 0.5),
    panel.grid.minor = element_line(color = "lightgrey", size = 0.25),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "lines")
  ) +
  ggforce::facet_zoom(xy = date >="2022-05-28"& date<= "2022-08-18", horizontal = FALSE) +
  geom_line(data = df_lin_freq, aes(x = date, y = `No. of sequences`*100, group = Lineages), color="white", linewidth = 1.1, alpha=0.7) +
  geom_line(data = df_lin_freq, aes(x = date, y = `No. of sequences`*100, color = Lineages), linewidth = 0.6, alpha=1) +
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "No. of sequences")) +
  scale_color_brewer(palette = "Dark2") +
  guides(
    color = guide_legend(ncol = 4, override.aes = list(linewidth = 1.5)),
    fill = guide_legend(ncol = 2, override.aes = list(alpha = 1))
  )

## Save the final figure with the hack applied
grid_fig_1a2 <- ggforce_hack_facet_zoom_x(fig_1a2)
ggsave("results/fig_1a.png", grid_fig_1a2, width = 10, height = 8, dpi = 300)
ggsave("results/fig_1a.pdf", grid_fig_1a2, width = 10, height = 8)


# Panel B: Hospital-acquired infections (HAIs) by date, each ward having one horizontal line, each case represented by a dot (stratified by staff and patient).
metadata_hai <- readxl::read_xlsx("data/hospital_data/HAI-28MAY2022-16AUG2022.xlsx")
metadata_cases <- readxl::read_xlsx("data/hospital_data/HK_nosocomial_metadata_2025-09-13_GS_031025.xlsx")
table(metadata_cases$case_type_updated)

# metadata_cases$change_of_type <- paste0(metadata_cases$case_type_old, " -> ", metadata_cases$case_type_updated)
# change_counts <- table(metadata_cases$change_of_type)
# paste0(names(change_counts), ": ", as.integer(change_counts))
# metadata_cases %>% filter(!change_of_type %in% c("patient -> Inpatient", "staff -> Staff")) %>% group_by(HOSPITAL, `WARD/CLUSTERS`, change_of_type) %>% summarise(N = n()) %>% arrange(desc(N)) %>% print(n=100)

metadata_cases$Ward_ID <- metadata_cases$`WARD/CLUSTERS`
metadata_cases$Ward_ID[!metadata_cases$Ward_ID %in% metadata_hai$`WARD/CLUSTERS`[metadata_hai$`Confirmed cluster`]] <- paste0(metadata_cases$Ward_ID[!metadata_cases$Ward_ID %in% metadata_hai$`WARD/CLUSTERS`[metadata_hai$`Confirmed cluster`]], "*") # add "*" to the ward ID of the clusters not reported as nosocomial infection

metadata_cases$sample_id <- sapply(metadata_cases$FASTA, function(x) as.character(strsplit(x, "/")[[1]][1]))
sort(table(metadata_cases$sample_id))
metadata_cases$date <- sapply(metadata_cases$FASTA, function(x) as.character(ymd(strsplit(x, "/")[[1]][3])))
metadata_cases$date_imputed <- is.na(metadata_cases$date)
metadata_cases$date[metadata_cases$date_imputed] <- metadata_cases %>% filter(is.na(date)) %>% transmute(date = (LAST_COL_DATE - `1ST_COL_DATE`)/2+`1ST_COL_DATE`) %>% mutate(date = as.character(as_date(date))) %>% pull(date)

metadata_cases$HOSPITAL <- factor(metadata_cases$HOSPITAL, metadata_cases %>% arrange(`1ST_COL_DATE`) %>% pull(HOSPITAL) %>% unique())
metadata_cases$Ward_ID <- factor(metadata_cases$Ward_ID, metadata_cases %>% arrange(desc(`1ST_COL_DATE`)) %>% pull(Ward_ID) %>% unique())

metadata_cases$standardized_datetime <- as.Date(as.numeric(metadata_cases$standardized_datetime), origin = "1899-12-30")
metadata_cases$API <- as.numeric(ymd(metadata_cases$date)-metadata_cases$standardized_datetime) %>% ceiling()
quantile(metadata_cases$API, na.rm=TRUE, probs=seq(0, 1, by=0.01))

metadata_cases$Is_psychiatric  <- metadata_cases$HOSPITAL %in% c("CPH", "TKP")

metadata_cases$case_type_new <- metadata_cases$case_type_updated
metadata_cases$case_type_new[metadata_cases$case_type_updated == "Inpatient" & (metadata_cases$API < 3)] <- "Inpatient (API < 3)"
metadata_cases$case_type_new[metadata_cases$case_type_updated == "Inpatient" & (metadata_cases$API >= 3)] <- "Inpatient (API >= 3)"
metadata_cases$case_type_new[metadata_cases$case_type_updated == "Inpatient" & is.na(metadata_cases$API)] <- "Inpatient (API unknown)"
metadata_cases$case_type_new[metadata_cases$case_type_new == "Staff" | metadata_cases$case_type_new == "Outpatient"] <- "Outpatient/Staff"

metadata_cases_plot <- metadata_cases %>% transmute(Date=ymd(date), Wards=Ward_ID, Hospital = HOSPITAL, `Case type`=case_type_new, `Hospital clusters` = HA_CLUSTERS, Lineages = LINEAGE)
metadata_cases_plot <- metadata_cases_plot %>% group_by(Date, `Hospital clusters`, Hospital, Wards, `Case type`, Lineages) %>% summarise(N = n())
metadata_cases_plot$`Case type` <- factor(metadata_cases_plot$`Case type`, levels = c("Inpatient (API < 3)", "Inpatient (API >= 3)", "Inpatient (API unknown)", "Outpatient/Staff"))

df_Hospital_anonymized <- tibble(Hospital=levels(metadata_cases_plot$Hospital), Hospital_anonymized = LETTERS[seq_along(unique(metadata_cases_plot$Hospital))])
write_csv(df_Hospital_anonymized, "data/hospital_data/Hospital_anonymized.csv")
metadata_cases_plot <- left_join(metadata_cases_plot, df_Hospital_anonymized, by = "Hospital")
write_csv(metadata_cases_plot %>% ungroup() %>% select(-Hospital), "data/hospital_data/metadata_cases_plot.csv")

dark2_colors <- brewer.pal(7, "Dark2")

# Create a dataset that identifies which date/ward combinations have multiple case types
overlap_data <- metadata_cases_plot %>%
  group_by(Date, Wards, Hospital_anonymized) %>%
  summarize(case_type_count = n_distinct(`Case type`), .groups = "drop") %>%
  filter(case_type_count > 1)

# Join this info back to original data
metadata_cases_plot_with_overlap <- metadata_cases_plot %>%
  left_join(overlap_data %>% select(Date, Wards, Hospital_anonymized, case_type_count), 
            by = c("Date", "Wards", "Hospital_anonymized"))

# Now create the plot with conditional jittering
fig_1b <- ggplot() +
  # Plot points without overlap (no jittering needed)
  geom_point(data = metadata_cases_plot_with_overlap %>% filter(is.na(case_type_count)),
             aes(x = Date, y = Wards, color = Lineages, shape = `Case type`, size = N), 
             alpha = 0.8) +
  # Plot points with overlap (apply jittering)
  geom_point(data = metadata_cases_plot_with_overlap %>% 
               filter(!is.na(case_type_count) & `Case type` == "Inpatient (API < 3)"),
             aes(x = Date, y = Wards, color = Lineages, shape = `Case type`, size = N),
             alpha = 0.8, position = position_nudge(y = 0.3)) +
  geom_point(data = metadata_cases_plot_with_overlap %>% 
               filter(!is.na(case_type_count) & `Case type` == "Inpatient (API >= 3)"),
             aes(x = Date, y = Wards, color = Lineages, shape = `Case type`, size = N),
             alpha = 0.8, position = position_nudge(y = 0.1)) +
  geom_point(data = metadata_cases_plot_with_overlap %>% 
               filter(!is.na(case_type_count) & `Case type` == "Inpatient (API unknown)"),
             aes(x = Date, y = Wards, color = Lineages, shape = `Case type`, size = N),
             alpha = 0.8, position = position_nudge(y = -0.1)) +
  geom_point(data = metadata_cases_plot_with_overlap %>% 
               filter(!is.na(case_type_count) & `Case type` == "Outpatient/Staff"),
             aes(x = Date, y = Wards, color = Lineages, shape = `Case type`, size = N),
             alpha = 0.8, position = position_nudge(y = -0.3)) +
  geom_vline(xintercept = as.Date("2022-05-31"), linetype = "dashed", 
             color = "#053C29", size = 0.8, alpha = 0.8) +
  geom_label_repel(x = as.Date("2022-05-31"), y = "AE*", 
                  label = "Relaxed visit schedule, 2022-05-31\n(2nd round, remaining major hospitals)", 
                  size = 2.5, alpha = 0.8, color = "#053C29", nudge_x = 0.5, 
                  show.legend = FALSE, 
                  data = metadata_cases_plot %>% filter(Wards == "AE*") %>% .[1,]) +
  scale_x_date(date_labels = "%b-%d", date_breaks = "15 days", expand = c(0, 0), 
               limits = c(as.Date("2022-05-28"), as.Date("2022-08-18"))) +
  scale_size_continuous(name = "No. of cases", breaks = c(1, 3, 5), 
                       labels = c(1, 3, 5), range = c(1, 6)) +
  scale_color_manual(values = c("BA.2.12.1" = dark2_colors[2], 
                               "BA.2.2" = dark2_colors[3], 
                               "BA.5.6" = dark2_colors[6])) +
  scale_shape_manual(values = c("Outpatient/Staff" = 5, 
                                "Inpatient (API unknown)" = 0,
                                "Inpatient (API < 3)" = 17, 
                                "Inpatient (API >= 3)" = 16)) +
  facet_grid(rows = vars(Hospital_anonymized), scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.background = element_rect(color = "black", fill = "white", size = 0.5)
  ) +
  ggtitle("b") +
  NULL

# Save the combined figure
combined_grid <- gridExtra::grid.arrange(grid_fig_1a2, ggplotGrob(fig_1b), ncol = 1, heights = c(1.5, 2))
ggsave("results/fig_1ab.png", combined_grid, width = 10, height = 16, dpi = 300)
ggsave("results/fig_1ab.pdf", combined_grid, width = 10, height = 16)

## We want to plot the distribution of API (admission-to-positive interval) for the cases in figure 1b
plot_API_all <- ggplot(metadata_cases %>% filter(!is.na(API))) +
  geom_histogram(aes(x = API, fill = Is_psychiatric), binwidth = 50, color = "black", alpha = 0.7) +
  scale_x_continuous(
    breaks = seq(0, 8000, by = 1000),
    labels = scales::comma
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = "Admission-to-positive interval (days)",
    y = "No. of inpatient cases"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey", size = 0.5),
    panel.grid.minor = element_line(color = "lightgrey", size = 0.25),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1, "lines"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "top"
  )
plot_API_all

plot_API_0_100 <- ggplot(metadata_cases %>% filter(!is.na(API) & API <= 100)) +
  geom_histogram(aes(x = API, fill = Is_psychiatric), binwidth = 1, color = "black", alpha = 0.7) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 10),
    labels = scales::comma
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = "Admission-to-positive interval (days)",
    y = "No. of inpatient cases"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey", size = 0.5),
    panel.grid.minor = element_line(color = "lightgrey", size = 0.25),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1, "lines"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "top"
  )
plot_API_0_100

ggsave("results/fig_S_API_all.pdf", plot_API_all, width = 7, height = 5)
ggsave("results/fig_S_API_0_100.pdf", plot_API_0_100, width = 6, height = 5)


## link fig_1b to a map
hkmap = readRDS("data/hospital_data/HKG_adm1.rds") 

# df_plot_1b_left_link <- metadata_cases_plot %>% ungroup() %>% group_by(Hospital, `Hospital clusters`) %>%
#   summarise(N = sum(N)) %>%
#   arrange(Hospital) %>% 
#   ungroup() %>%
#   mutate(
#     # bump_y_start = -90,
#     rank = row_number(),
#     bump_y_start = scales::rescale(-rank, to = c(22.13, 22.6)),
#     # bump_x_start = normalize(-rank, range = c(-180, 180), method = "range")
#     bump_x_start = 113.8
#   )

df_pos <- read_csv("data/hospital_data/hong_kong_hospitals_with_districts.csv")
# df_pos <- left_join(df_pos, df_Hospital_anonymized, "Hospital")
# df_plot_1b_left_link <- left_join(df_plot_1b_left_link, df_pos, "Hospital")
# df_plot_1b_left_link <- df_plot_1b_left_link %>% filter(Hospital %in% metadata_hai$HOSPITAL[metadata_hai$`Confirmed cluster`])

hkmapdf = fortify(hkmap)
hkmapmeta <- read_csv("data/hospital_data/Hong_Kong_Districts_with_HA_Clusters.csv")
hkmapdf = merge(hkmapdf, hkmapmeta, by.x="id", by.y="mapid")

p_fig_1b_right <- ggplot(hkmapdf) +
  geom_polygon(aes(long, lat, group=group, fill=`HA Cluster`)) + 
  # geom_sigmoid(data = df_plot_1b_left_link, 
  #   aes(x = `Longitude (E)`, y = `Latitude (N)`, xend = bump_x_start, yend = bump_y_start, group=Hospital), 
  #   alpha = 0.9, smooth = 12, size = 1.2, direction = "x", color="white") +
  # geom_sigmoid(data = df_plot_1b_left_link, 
  #   aes(x = `Longitude (E)`, y = `Latitude (N)`, xend = bump_x_start, yend = bump_y_start, group=Hospital), 
  #   alpha = 0.8, smooth = 12, size = 0.8, direction = "x") +
  geom_label_repel(data = df_pos, 
                   aes(x = `Longitude (E)`, y = `Latitude (N)`, label = `Full Name`), 
                   size = 3.5, 
                   alpha = 0.8,
                   box.padding = 0.01, 
                   point.padding = 0.3,
                   force = 1,
                   min.segment.length = 0,
                   max.overlaps = 20) +
  theme_void() +
  scale_fill_manual(values = met.brewer("Degas", n=7, type="discrete")) +
  # coord_sf(xlim = c(-150, 180), ylim = c(-60, 80), expand = TRUE)+
  theme(
    legend.position = "bottom"
    ) + 
  ggtitle("d")+
  NULL
ggsave("results/fig_1d.png", p_fig_1b_right, width = 10, height = 16*(2/3.5), dpi = 300)
ggsave("results/fig_1d.pdf", p_fig_1b_right, width = 10, height = 16*(2/3.5))

# Figure 1c
## In figure 1c, we plot the mobility and hospital-related policy data for Hong Kong in 2022. Mobility data is from three sources: Google Community Mobility Report (https://www.google.com/covid19/mobility/index.html?hl=en), governmental cross-border and local mobility data (https://www.td.gov.hk/en/transport_in_hong_kong/transport_figures/monthly_traffic_and_transport_digest/2022/202212/index.html, https://data.gov.hk/en-data/dataset/hk-immd-set5-statistics-daily-passenger-traffic).

data_google <- read_csv("data/mobility_data/2022_HK_Region_Mobility_Report.csv")
names(data_google)
data_google <- data_google %>% transmute(date, `Retail and recreation` = `retail_and_recreation_percent_change_from_baseline`, `Grocery and pharmacy` = `grocery_and_pharmacy_percent_change_from_baseline`, `Parks` = `parks_percent_change_from_baseline`, `Transit stations` = `transit_stations_percent_change_from_baseline`, `Workplaces` = `workplaces_percent_change_from_baseline`, `Residential` = `residential_percent_change_from_baseline`)

data_google_long <- data_google %>% pivot_longer(cols = -date, names_to = "Location", values_to = "Mobility") %>% filter(!is.na(Location))
data_google_long$Type <- NA
data_google_long$Type[data_google_long$Location %in% c("Retail and recreation", "Grocery and pharmacy", "Parks")] <- "Public"
data_google_long$Type[data_google_long$Location %in% c("Transit stations", "Workplaces", "Residential")] <- "Private"
dir.create("data/mobility_data/parsed", showWarnings = FALSE)
write_csv(data_google_long, "data/mobility_data/parsed/google_mobility_long.csv")

fig_1c_1 <- ggplot(data_google_long) +
  geom_path(aes(x = date, y = Mobility, color = Location), position = "identity", alpha=0.8, size = 1) +
  geom_vline(xintercept = as.Date("2022-05-28"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  geom_vline(xintercept = as.Date("2022-08-18"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  scale_x_date(date_labels = "%b-%d", date_breaks = "30 days", expand = c(0.01, 0), limits = c(as.Date("2022-02-21"), as.Date("2022-08-18"))) +
  # highlight the period of interest
  labs(
    title = "Google Community Mobility",
    x = "Date",
    y = "Mobility change (%)",
    fill = "Location"
  ) +
  scale_color_manual(values = met.brewer("Austria", n=6, type="discrete")) +
  facet_wrap(~Type, ncol = 2, scales = "free_y", strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = unit(0.2, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey", size = 0.5),
    panel.grid.minor = element_line(color = "lightgrey", size = 0.25),
    # axis.text.x = element_text(angle = 30, hjust = 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(nrow = 1))
fig_1c_1

data_cross_border <- read_csv("data/mobility_data/statistics_on_daily_passenger_traffic.csv")
data_cross_border_2022 <- data_cross_border %>% mutate(Date=dmy(Date)) %>% filter(Date >= as.Date("2022-01-01") & Date <= as.Date("2022-12-31"))
unique(data_cross_border_2022$`Control Point`)
data_cross_border_2022$Border <- data_cross_border_2022$`Control Point` 
data_cross_border_2022$Border[data_cross_border_2022$Border %in% c("Express Rail Link West Kowloon", "Hung Hom", "Lo Wu", "Lok Ma Chau Spur Line", "Heung Yuen Wai", "Lok Ma Chau", "Man Kam To", "Sha Tau Kok", "Shenzhen Bay", "Hong Kong-Zhuhai-Macao Bridge", "China Ferry Terminal")] <- "Mainland border"
data_cross_border_2022$Border[data_cross_border_2022$Border %in% c("Harbour Control", "Kai Tak Cruise Terminal", "Macau Ferry Terminal", "Tuen Mun Ferry Terminal")] <- "Harbour"
data_cross_border_2022 <- data_cross_border_2022 %>% group_by(Date, Border, `Arrival / Departure`) %>% summarise(Total = sum(Total)) %>% filter(Date>="2022-02-11", Date<= "2022-09-01")
write_csv(data_cross_border_2022, "data/mobility_data/parsed/cross_border_2022.csv")

fig_1c_2 <- ggplot(data_cross_border_2022)+
  geom_path(aes(x = Date, y = Total, color = Border), position = "identity", alpha=0.8, size = 1) +
  geom_vline(xintercept = as.Date("2022-05-28"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  geom_vline(xintercept = as.Date("2022-08-18"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  scale_x_date(date_labels = "%b-%d", date_breaks = "30 days", expand = c(0.01, 0), limits = c(as.Date("2022-02-21"), as.Date("2022-08-18"))) +
  # highlight the period of interest
  labs(
    title = "Cross-border passenger traffic",
    x = "Date",
    y = "No. of passengers",
    fill = "Location"
  ) +
  scale_color_manual(values = met.brewer("Kandinsky", n=3, type="discrete")) +
  scale_linetype_manual(values = c("Arrival" = "solid", "Departure" = "dashed")) +
  facet_wrap(~`Arrival / Departure`, ncol = 2, scales = "free_y", strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey", size = 0.5),
    panel.grid.minor = element_line(color = "lightgrey", size = 0.25),
    # axis.text.x = element_text(angle = 30, hjust = 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "lines")
  )
fig_1c_2

data_local_transport <- readxl::read_excel("data/mobility_data/public_transport.xlsx")
data_local_transport <- data_local_transport %>% mutate(`Taxi and Other Services` = Taxi + `Residential Services` + `MTR buses`, Buses = `Franchised Buses` + `Public Light Buses`) %>% select(-Taxi, -`Residential Services`, -`MTR buses`, -`Franchised Buses`, -`Public Light Buses`)
data_local_transport <- left_join(
  tibble(Date = seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by = "day"),
  Year_Month = paste0(year(Date), "-", month(Date))),
  # Create a reference date (first day of each month) for joining
  data_local_transport %>% mutate(Year_Month = paste0(Year, "-", Month)) %>% select(-Year, -Month) 
)
data_local_transport_plot <- data_local_transport %>% select(-Year_Month) %>% pivot_longer(cols = -Date, names_to = "Transport", values_to = "No. of passengers") %>% filter(Date>="2022-02-11", Date<= "2022-09-01")
write_csv(data_local_transport_plot, "data/mobility_data/parsed/local_transport_long.csv")

fig_1c_3 <- ggplot(data_local_transport_plot %>% filter(Date>="2022-02-11", Date<= "2022-09-01")) +
  geom_path(aes(x = Date, y = `No. of passengers`, color = Transport), position = "identity", alpha=0.8, size = 1) +
  geom_vline(xintercept = as.Date("2022-05-28"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  geom_vline(xintercept = as.Date("2022-08-18"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  geom_vline(xintercept = as.Date("2022-05-06"), linetype = "dashed", color = "#053C29", size = 0.8, alpha = 0.8) + # Relaxed visit scheduel first round
  geom_label_repel(aes(x = as.Date("2022-05-06"), y = 2000, label = "Relaxed visit schedule, 2022-05-06\n(1st round, 26 public hospitals)"), size = 2.5, alpha=0.8, color = "#053C29", nudge_x = -0.5, show.legend = FALSE, data = tibble()) +
  geom_vline(xintercept = as.Date("2022-05-31"), linetype = "dashed", color = "#053C29", size = 0.8, alpha = 0.8) + # Relaxed visit scheduel second round
  geom_label_repel(aes(x = as.Date("2022-05-31"), y = 2000, label = "Relaxed visit schedule, 2022-05-31\n(2nd round, remaining major hospitals)"), size = 2.5, alpha=0.8, color = "#053C29", nudge_x = 0.5, show.legend = FALSE, data = tibble()) +
  scale_x_date(date_labels = "%b-%d", date_breaks = "15 days", expand = c(0.01, 0), limits = c(as.Date("2022-02-21"), as.Date("2022-08-18"))) +
  # highlight the period of interest
  labs(
    title = "Local transport",
    x = "Date",
    y = "No. of passengers (thousands)",
    fill = "Location"
  ) +
  scale_color_manual(values = met.brewer("Lakota", n=4, type="discrete")) +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey", size = 0.5),
    panel.grid.minor = element_line(color = "lightgrey", size = 0.25),
    # axis.text.x = element_text(angle = 30, hjust = 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "lines")
  )
fig_1c_3

fig_1c <- gridExtra::grid.arrange(fig_1c_1, fig_1c_2, fig_1c_3, ncol = 1)
ggsave("results/fig_1c.png", fig_1c, width = 10, height = 8, dpi = 300)
ggsave("results/fig_1c.pdf", fig_1c, width = 10, height = 8)

# Other statistics

HA_reported_data <- readxl::read_xlsx("data/hospital_data/Hospital Authority announces clusters of nosocomial infections of COVID-19.xlsx")

names(HA_reported_data)
HA_reported_data_nosocomial_infection <- HA_reported_data %>% filter(Type=="Nosocomial infection") %>% mutate(sequenced = !is.na(`Ward ID in our study`)) %>% mutate(Total_cases= as.numeric(`Nosocomial infection patient`) + as.numeric(`Nosocomial infection staff`))
HA_reported_data_nosocomial_infection %>% group_by(sequenced) %>% summarise(n = n()) # the number of nosocomial infections clusters we have sequenced
HA_reported_data_nosocomial_infection %>% group_by(sequenced) %>% summarise(Total_cases = sum(Total_cases, na.rm=TRUE)) %>% ungroup() %>% mutate(Prop_cases = Total_cases/sum(Total_cases, na.rm=TRUE)) # the number and prop of nosocomial infections clusters we have sequenced

HAI_own_data <- readxl::read_xlsx("data/hospital_data/HAI-28MAY2022-16AUG2022.xlsx")
HAI_own_data %>% group_by(`Confirmed cluster`) %>% summarise(n = n(), Total_cases=sum(`TOTAL_CASES#`)) # the number of confirmed nosocomial infections clusters we have sequenced
HAI_own_data %>% filter(`Confirmed cluster`) %>% group_by(HOSPITAL) %>% summarise(n = n(), Total_cases=sum(`TOTAL_CASES#`)) # the number of confirmed nosocomial infections clusters we have sequenced
table(HAI_own_data$`Reason for not identified as nosocomial infection / Bootstrap support for confirmed clusters`)

metadata_hai %>% filter(`Confirmed cluster`) %>% .$HA_CLUSTERS %>% table() # the number of confirmed nosocomial infections clusters we have sequenced

table(metadata_cases$LINEAGE)
table(metadata_cases$LINEAGE)/nrow(metadata_cases)*100 # the proportion of BA.2.2 in our data

metadata_cases %>% group_by(HOSPITAL, Ward_ID, LINEAGE) %>% summarise(n = n()) %>% ungroup() %>% group_by(HOSPITAL, Ward_ID) %>% summarise(n = n(), Lineages = paste(unique(LINEAGE), collapse = ",")) %>% arrange(desc(n), desc(Lineages)) # the number of lineages in each ward

