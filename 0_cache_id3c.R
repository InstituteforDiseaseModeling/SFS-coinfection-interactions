### ------------------------------------------------------------------------------------------------
### PURPOSE: pull from the id3c DB and save locally for analysis.
### NOTE: - This script only relevant for those with access for SFS ID3C database, if that is not you
###         skip directly to (2_manuscript_analysis.R) script which uses a pre-loaded dataset
###       - For this to work your IP must be white-listed by the ID3C team (contact team)
###         Additionally, you will need to create a directory in your repo called ./credentials. 
###         Inside this directory you will need a populated .pgpass and a .pg_service.conf file
###       - Explore DB schema here: https://backoffice.seattleflu.org/metabase/browse/5/schema/shipping
### ------------------------------------------------------------------------------------------------



# Set  local repo location. You'll have to change this here for your code to run. 
repo_root <- <<WORKING DIRECTORY PATH HERE>> 
repo      <- sprintf('%s/SFS-coinfection-interactions/',     repo_root)
setwd(repo)

# Load in the libraries 
libs <- c('RPostgres', 'DBI', 'data.table') 
for(lib in libs) library(lib, character.only = TRUE)


# Set up permissions, pulled from your credentials folder
credentials_path <- file.path(repo,'credentials')
Sys.setenv(PGSERVICEFILE = file.path(credentials_path, ".pg_service.conf"),
           PGPASSFILE    = file.path(credentials_path, ".pgpass"))

# create a connection to the database
conn <- DBI::dbConnect(RPostgres::Postgres(), service = "seattleflu-production", timezone = NULL)

## Grab data from a few different database views and save them as .RDS files within your repo
## (NOTE that .gitignore ignores these files to avoid accidentally uploading them)

# all COV lab results
d <- dbGetQuery(conn, "select distinct * from shipping.hcov19_observation_v1;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_db.RDS', repo))

# results from scan, including those that just arrived at the lab (not really used anymore)
d <- dbGetQuery(conn, "select distinct * from shipping.incidence_model_observation_v2;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_db_scan_realtime.RDS', repo))

# questionnaire data
d <- dbGetQuery(conn, "select distinct * from shipping.fhir_encounter_details_v2;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_db_questionnaire.RDS', repo))

# follow-up survey
d <- dbGetQuery(conn, "select distinct * from shipping.scan_follow_up_encounters_v1;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_db_fu.RDS', repo))

# multipathogen stuff 
d <- dbGetQuery(conn, "select distinct * from shipping.presence_absence_result_v2;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_all_pathogens.RDS', repo))

# all obs
d <- dbGetQuery(conn, "select distinct * from shipping.incidence_model_observation_v3;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_all_obs.RDS', repo))

# all data for multipathogen 
d <- dbGetQuery(conn, "select distinct * from shipping.observation_with_presence_absence_result_v1;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_all_mp.RDS', repo))

# join on this! starting point all sample
d <- dbGetQuery(conn, "select distinct * from shipping.sample_with_best_available_encounter_data_v1;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_all_samples.RDS', repo))

## new full view (added later):
## https://seattle-flu-study.slack.com/archives/CCAA9RBFS/p1608752968265700
d <- dbGetQuery(conn, "select distinct * from shipping.observation_with_presence_absence_result_v3;")
d <- as.data.table(d)
saveRDS(d, sprintf('%s/id3c_cache/cache_observation_with_presence_absence_result_v3.RDS', repo))
