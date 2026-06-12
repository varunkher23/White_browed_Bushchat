library(janitor)
library(dplyr)

form <- read.csv("habitat_23-24/form-1__gib-habitat-ecology.csv") %>%
  clean_names() %>%
  select(-utm_northing_9_start_location:-x10_start_location_ma,-utm_northing_11_end_location:-x13_temperature,
         -created_at,-uploaded_at, -x24_notes:-x52_vegetation_shrubs)

trees <- read.csv("habitat_23-24/branch-1__vegetation-trees.csv") %>%
  clean_names() %>%
  mutate(foliar_volume = x47_number * 1.57 * (x44_height_in_m - x46_canopy_start_heig) * (x45_diameter_in_m / 2) * (x45_diameter_in_m / 2)) %>% ### Half-Oval
  mutate(x45_diameter_in_m_total = x47_number * x45_diameter_in_m) %>%
  group_by(ec5_branch_owner_uuid, x43_species) %>%
  summarise(dbh_1ha = sum(x45_diameter_in_m_total,na.rm = T),
            foliar_volume_1h = sum(foliar_volume, na.rm = T),
            n_trees = sum(x47_number, na.rm=T))

shrubs <- read.csv("habitat_23-24/branch-2__vegetation-shrubs.csv") %>%
  clean_names() %>%
  mutate(area_covered = x56_diameter_in_cm / 100 * x57_number) %>%
  mutate(volume = (x56_diameter_in_cm / 100) * x57_number * (x55_height_in_cm / 100)) %>%
  group_by(ec5_branch_owner_uuid, x54_species) %>%
  summarise(area_covered_1ha = sum(area_covered,na.rm=T) * 2.5,
            volume_1ha = sum(volume, na.rm =T) * 2.5)

understory <- read.csv("habitat_23-24/branch-5__understory.csv") %>%
  clean_names()%>%
  mutate(
    midpoint_height_class = case_when(
      x31_hight_class == "1-20" ~ 10,
      x31_hight_class == "21-40" ~ 30,
      x31_hight_class == "41-60" ~ 50,
      x31_hight_class == "61-80" ~ 70,
      x31_hight_class == "81-100" ~ 90,
      x31_hight_class == "101-120" ~ 110,
      x31_hight_class == "121-140" ~ 130,
      x31_hight_class == "141-160" ~ 150,
      x31_hight_class == "161-180" ~ 170,
      x31_hight_class == "181-200" ~ 190
    ),
    x28_species = gsub(" ", "_", gsub("\"", "", sub("_.*", "", x28_species)))
  )%>%
  group_by(ec5_branch_owner_uuid, x27_quadrat, x28_species, midpoint_height_class) %>%
  summarise(
    mean_cover = mean(as.numeric(sub("-.*", "", x32_percent_cover_cla)), na.rm = TRUE),
    n_quadrats = n(),
    .groups = "drop"
  )
