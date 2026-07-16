library(move2)
library(sf)
library(dplyr)
library(tigris)
library(motus)
library(purrr)
library(here)
library(openxlsx)


# MOVEBANK

# get MN boundary
mn <- states(cb = TRUE) |>
  filter(NAME == "Minnesota") |>
  st_transform(4326)


username <- Sys.getenv("MOVEBANK_USERNAME")
password <- Sys.getenv("MOVEBANK_PASSWORD")

# Set up Movebank credentials and get data
move2::movebank_store_credentials(
  username = username,
  password = password
)


all_studies_sf <- st_as_sf(all_studies, sf_column_name = "main_location")

mn <- states(cb = TRUE) |>
  filter(NAME == "Minnesota") |>
  st_transform(st_crs(all_studies_sf))

mn_studies <- all_studies_sf |>
  filter(!st_is_empty(main_location)) |>
  st_filter(mn)


studies_in_archive <- move2::movebank_download_study_info(i_am_collaborator=TRUE)


total_movebank_studies <- full_join(mn_studies, studies_in_archive, by = "id") |>
  mutate(
    across(
      ends_with(".x"),
      ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column()))),
      .names = "{sub('.x$', '', .col)}"
    )
  ) |>
  select(-ends_with(".x"), -ends_with(".y")) |>
  mutate(
    in_mama_archive = ifelse(id %in% studies_in_archive$id, TRUE, FALSE),
    movebank = TRUE,
    movebank_id = id,
    movebank_url = paste0("https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study", movebank_id),
    lon = st_coordinates(main_location.x)[, 1],
    lat = st_coordinates(main_location.x)[, 2]
  ) |>
  select(-main_location.x)


glimpse(total_movebank_studies)

# Motus


# let this one actually download/build the schema
sql_motus <- tagme(176, new = TRUE)   # update defaults to TRUE

recv_deps <- tbl(sql_motus, "recvDeps") |> collect()
projs <- tbl(sql_motus, "projs") |> collect()

# spatial filter to MN
recv_sf <- recv_deps |>
  filter(!is.na(latitude), !is.na(longitude), latitude != 0, longitude != 0) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


mn_motus_projects <- recv_sf |>
  st_filter(mn) |>
  left_join(projs, by = c("projectID" = "id")) |>
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry() |>
  group_by(projectID, name.y) |>
  summarise(
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    n_receivers = n(),
    .groups = "drop",
  ) |>
  mutate(
    motus = TRUE,
    motus_id = projectID,
    name = name.y,
    motus_url = paste0("https://motus.org/data/project?id=", motus_id),
    in_mama_archive = FALSE
  ) |>
  select(-c(projectID, name.y))


glimpse(mn_motus_projects)  


all_data <- bind_rows(total_movebank_studies, mn_motus_projects) |>
  st_drop_geometry() |>
  dplyr::select(name,
                principal_investigator_name,
                timestamp_first_deployed_location,
                timestamp_last_deployed_location,
                number_of_deployed_locations,
                number_of_individuals,
                sensor_type_ids,
                taxon_ids,
                in_mama_archive,
                movebank_id,
                motus_id,
                movebank_url,
                motus_url,
                lon,
                lat,
                n_receivers)



glimpse(all_data)
write.csv(all_data, here("data/studies.csv"))

write.xlsx(all_data, here("data/studies.xlsx"))

# LCCMR Studies

# load from airtable
airtable_mm_studies <- read.csv(here("data/movebank_motus_airtable.csv"))
glimpse(airtable_mm_studies)

lccmr_studies <- read.csv(here("data/lccmr_studies_airtable.csv"))
glimpse(lccmr_studies)



