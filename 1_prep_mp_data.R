### ----------------------------------------------------------------------------------------
### Purpose: Prepare data sourced from 0_cache_id3c.R script preceeding 2_manuscript_analysis.R
### NOTE:    - This script only relevant for those with access for SFS ID3C database, if that is not you
###            skip directly to (2_manuscript_analysis.R) script which uses a pre-loaded dataset
### ------------------------------------------------------------------------------------------------


# Load in the libraries 
library(here)
library(lubridate)
library(jsonlite)
library(data.table)


# set directory
repo_root <- <<WORKING DIRECTORY PATH HERE>> 
repo      <- sprintf('%s/SFS-coinfection-interactions/',     repo_root)
setwd(repo)



## ---------------------------------------
## LOAD DATA. Note: bringing together several DB views here. 

# all samples. High level metadata on all samples
allsamp <- readRDS('id3c_cache/cache_all_samples.RDS')

# all covid
covid <- readRDS('id3c_cache/cache_db.RDS')

# all (non-covid) presense/absence
aps <- readRDS('id3c_cache/cache_all_pathogens.RDS') 

# non-covid metadata
ncmd <- readRDS('id3c_cache/cache_all_obs.RDS')   

# this view was added later, its got everything, grabbing swab from it
tt <- readRDS('id3c_cache/cache_observation_with_presence_absence_result_v3.RDS')   


## ---------------------------------------
## RE-FORMAT  DATA 

#1 shipping.observation_with_presence_absence_result_v* performs an inner join on encounters and samples. 
#2 shipping.presence_absence_result_v* has no hcov19 result (purposefully)
#3 shipping.hcov19_observation_v1 has hcov19 result and some metadata.
#4 need to join non-COVID metadata on age/sex/geography

## 1. allsamp identifies all samples we should have info of sort for so this is our 'master list'
d <- copy(allsamp)
d <- d[!is.na(sample)] # sample_id is in there too, but cant have sample being NA

# trim some variables I do not need
d[, c('best_available_site_id','best_available_site_category'):= list(NULL,NULL)]

# TEST: confirm that all sample are unique
stopifnot(length(unique(d$sample))==nrow(d))

## 2. combine presence.absence for all samples and covid
labvars <- c('sample','target','organism','present','ct')

# covid requires a bit of cleaning first
labc <- copy(covid)
labc[, present := hcov19_present_bbi]
labc[is.na(present), present := hcov19_present_uw] 
labc[, target   := 'SARS_CoV_2']
labc[, organism := 'Coronavirus_2019']

# keep CT from covid (take the mean of the values)
# Note: first two are Orf1b, second to are S. Using orf1b align with JAMA Peds Paper (Chung 2021)
labc[, ct := as.character(crt_values)]
labc[, ct := gsub('\\{','',gsub('\\}','',ct))]
x <- do.call('rbind',strsplit(labc$ct,','))
labc$ct <- round(rowMeans(matrix(as.numeric(x),ncol=4)[,1:2],na.rm=T))
labc[present!=TRUE, ct := NA]

# which covid samples were on the open array?
labc[, device_oa := grepl('OpenArray', hcov19_device_bbi)]
nrow(labc[device_oa==TRUE&present==TRUE])
labc[device_oa==FALSE & sample_id %in% labc[device_oa==TRUE&present==TRUE]$sample_id]


# keep CT from details in the MP data.. parse JSON.. slow, takes a few minutes
aps[,tmpid:=1:.N]
tmp <- aps[present==TRUE]
tmp[, c('device','ct1','ct2','b','c') := tstrsplit(details,'"crt\": \"')]
tmp[, ct11 := (as.numeric(substr(ct1,1,5)))]
tmp[, ct22 := (as.numeric(substr(ct2,1,5)))]
tmp[, ct := (ct11+ct22)/2]
tmp <- tmp[,c('tmpid','ct')]
aps[, device_oa := grepl('OpenArray',details)]
aps <- merge(aps,tmp,by='tmpid',all.x=TRUE)

# bring the two together
aps[, organism := as.character(organism)]
aps <- aps[organism!='Human_coronavirus.2019'] # avoid duplication with labc after ncov added to this view
lab <- rbind(copy(aps)[, labvars, with=FALSE], labc[, c(labvars), with=FALSE])
lab <- lab[order(sample)]

# clean up some stuff
lab <- lab[!is.na(sample)]  # remove rows with no sample number
lab <- lab[!is.na(present)] # remove rows without a lab result


## 3. individual-level (non encounter-level) metadata if available (non-covid)

# clean up age data. Use various age binnings
ncmd[, age := age_range_fine_lower + ((age_range_fine_upper-age_range_fine_lower)/2)]

ncmd[, tract := residence_census_tract]
ncmd <- ncmd[, c('sample','sex','age', 'tract')] #, 'flu_shot', 'race', 'symptoms')]
ncmd <- unique(ncmd)
ncmd[, source := 'nc']

## 4. covid meta-data
# names
covid[, tract := census_tract]

# clean up age data. Use various age binnings
covid[age_range_fine_lower == '00:00:00', arl := 0]
covid[grep('year',age_range_fine_lower),  arl := as.numeric(substr(age_range_fine_lower,1,2))]
covid[grep('mon', age_range_fine_lower),  arl := as.numeric(substr(age_range_fine_lower,1,1))/12]
covid[grep('year',age_range_fine_upper),  aru := as.numeric(substr(age_range_fine_upper,1,2))] # TODO add check for 100+
covid[grep('mon', age_range_fine_upper),  aru := as.numeric(substr(age_range_fine_upper,1,1))/12]
covid[, age := arl + ((aru-arl)/2)]

#add on questionnaire
cmd <- covid[, c('sample','sex','age', 'tract')]
cmd[, source := 'c']

# bring together covid and non-covid metadata and drop any duplicates
# we get actual age in the non-covid source, so go with that one
md <- rbind(cmd, ncmd)
md <- md[order(sample,-source)]
md <- md[duplicated(sample)==FALSE]
md <- md[!is.na(sample)]

## SWAB, collection (dry or utm), symptoms 
tmp <- unique(tt[!is.na(sample),c('sample','swab_type','collection_matrix','symptoms','target',
                                  'individual','manifest_origin','device')])
tmp[target %in% c('nCov','COVID-19','Human_coronavirus.2019'), target := 'SARS_CoV_2']
lab[target %in% c('nCov','COVID-19','Human_coronavirus.2019'), target := 'SARS_CoV_2']

## 5. bring together the master dataset now, lab and metadata info, finalize and clean data
d <- merge(d, md,  by = 'sample', all.x = TRUE)
d <- merge(d, lab, by = 'sample', all.x = TRUE)
d <- merge(d, tmp, by = c('sample','target'), all.x = TRUE)

## Finalize some cleaning things
# rename some things
setnames(d, c('best_available_encounter_date','best_available_site','best_available_site_type'),
            c('date','site','site_type'))

# remove where present or the organism is missing
d <- d[!is.na(present)]
d <- d[!is.na(target)]

# final variables
finalvars <- c('sample','sample_id','target','organism','present','date',
               'site',  'site_type', 'tract','age','sex','ct','device',
               'swab_type','collection_matrix','symptoms','individual','manifest_origin')
d <- d[, finalvars, with = FALSE]
d <- d[order(date,sample)]
d[, date := ymd(date)]

# spot fix, incorrect date
d[sample=='8c3d94e2-2ff0-491c-b5df-e66e5945beab', date := ymd('2019-08-20')]

# name missing sites so they dont drop yet. They are missing other metadata so get excluded in a future script
d[is.na(site), site := 'NA']

# filter out all site=='self_test'  as those are RDT results not PCR so no CT and thus irrelevant for analysis
d <- d[(site != 'self-test')]

# drop duplicates sars2 rows. privilege taqman cts if dual, some duplicated rows with NA device should be removed as well
d[,iid := 1:.N] #rowid
tmp<-d[target =='SARS_CoV_2']
tmp[,n:=.N,by=sample]
tmp[, both := all(c("OpenArray","TaqmanQPCR") %in% device), by = sample]
dropiids <- tmp[both==TRUE & device=='OpenArray']$iid
tmp <- tmp[!iid %in% dropiids]
tmp[,n:=.N,by=sample]
tmp[n==2]
dropiids <- c(dropiids,tmp[is.na(device)&n==2]$iid)
d <- d[!iid %in% dropiids]

## ##################################################
# clean up target and organisms that are messy in the original data
# see how Mike did in in the past here:  
#   https://github.com/seattleflu/incidence-mapper/blob/ef4540fcfa40f9f96784d2ffcb6be937c0248f19/dbViewR/R/selectFromDB.R#L122
# Also, informed by some notes from Lea:
#   https://docs.google.com/spreadsheets/d/1MolKJawXrPsr86TLXMaVnGqmPFNpOexR_NRXmWxl-6A/edit

# Lea says, always combine AdV1 and AdV2
d[tolower(target) %in% tolower(c('Adenovirus_pan_2','Adenovirus_pan_1','AdV_1of2','AdV_2of2')), organism := 'AdV']
d[tolower(target) %in% tolower(c('http://snomed.info/id/440930009')), organism := 'AdV'] 

# Lea says, always combine RV1 and RV2
d[tolower(target) %in% tolower(c('12 Rhinovirus_pan_2','11 Rhinovirus_pan_1','RV_1of2','RV_2of2')), organism := 'RV']
d[tolower(target) %in% tolower(c('http://snomed.info/id/440925005')), organism := 'RV']

# Streptococcus pneumonaie
d[tolower(target) %in% tolower(c('S. pneumoniae_APZTD4A','APZTD4A')), organism := 'Streptococcus_pneumoniae'] 
d[tolower(target) %in% tolower(c('APZTD4A','S. pneumoniae_APZTD4A','S. pneumoniae')), organism := 'Streptococcus_pneumoniae']

# C and M Pneumo
d[tolower(target) %in% tolower(c('C. pneumoniae_AI1RW2H')), organism := 'Chlamydophila_pneumoniae']
d[tolower(target) %in% tolower(c('AI5IRK5','M. pneumoniae_AI5IRK5','M. pneumoniae')), organism := 'Mycoplasma_pneumoniae']
d[tolower(target) %in% tolower(c('M. pneumoniae_AI5IRK5','AI5IRK5')), organism := 'Mycoplasma_pneumoniae'] 

# Bortedella Pertussis
d[tolower(target) %in% tolower(c('AI20U8U')), organism :='Bordetella_pertussis']

# Flu B and C
d[tolower(target) %in% tolower(c('flu_A_pan','Flu_a_pan',"http://snomed.info/id/181000124108")), organism := 'IAV']
d[tolower(target) %in% tolower(c('Influenza_B','Flu_b_pan','http://snomed.info/id/441345003')), organism := 'IBV']
d[tolower(target) %in% tolower(c('AP324NU')), organism := 'ICV']

# Combining the 4 coronaviruses. Targets have changed a bit over time, so just combine for now
#d[tolower(target) %in% tolower(c('CoV_229E_CoV_OC43','CoV_HKU1_CoV_NL63','CoV_OC43','CoV_229E','CoV_HKU1','CoV_NL63')), organism := 'CoV']
d[tolower(target) %in% tolower(c('CoV_HKU1_CoV_NL63','CoV_HKU1','CoV_NL63')), organism := 'CoV..HKU1_NL63']
d[tolower(target) %in% tolower(c('CoV_229E_CoV_OC43','CoV_OC43','CoV_229E')), organism := 'CoV..229E_OC43']

# Similar thing with 4 parainfluenzas
#d[tolower(target) %in% tolower(c('hPIV1','hPIV2','hPIV1_hPIV2','hPIV3','hPIV4','hPIV3_hPIV4')), organism := 'hPIV']
d[tolower(target) %in% tolower(c('hPIV1','hPIV2','hPIV1_hPIV2')), organism := 'hPIV..1_2']
d[tolower(target) %in% tolower(c('hPIV3','hPIV4','hPIV3_hPIV4')), organism := 'hPIV..3_4']

# Mumps
d[tolower(target) %in% tolower(c('APKA3DE')), organism := 'Mumps']

# Enterovirus is kind of a weird one,
d[tolower(target) %in% tolower(c('AP7DPVF','EnterovirusA_B 1_AP7DPVF','Ev_pan')), organism :='EV_pan']
d[tolower(target) %in% tolower(c('Enterovirus-D_APFVK4U','enterovirus-D_APFVK4U')), organism := 'Enterovirus.D.68']
d[tolower(target) %in% tolower(c('Enterovirus.D.68')), organism := 'Enterovirus.D.68'] 

# Notes on entero + RV from lea:
#  Entero_pan + RV1 or RV2 == Rhino
d[organism=='EV_pan'&present==TRUE, tmp_ERV := 'EPAN']
d[organism=='RV'&present==TRUE, tmp_ERV := 'RV']
d[, tmp_ERV2 := all(c("EPAN","RV") %in% tmp_ERV), by = sample]
d[tmp_ERV2==TRUE & organism=='EV_pan', present := FALSE]
d[tmp_ERV2==TRUE & organism=='EV_pan', ct      := NA]

# Entero_Pan alone or with EVD68 (high Ct) == enterovirus.
#  Threshold delta of 10 (lea)
entero_thresh <- 10
d[organism=='EV_pan'&present==TRUE,       tmp_E := 'EPAN']
d[organism=='Enterovirus.D.68'&present==TRUE, tmp_E := 'Enterovirus.D.68']
d[, tmp_E2 := all(c("EPAN","Enterovirus.D.68") %in% tmp_E), by = sample]
tmp <- d[tmp_E2==TRUE&organism%in%c('EV_pan','Enterovirus.D.68') & !is.na(ct)][order(sample,organism)]
tmp[,i:=1:.N,by=sample]
tmp[,ctdiff:=as.numeric(ct)-as.numeric(data.table::shift(ct,1,type='lag')),by=sample]
Esamp <- tmp[ctdiff>entero_thresh]$sample # NOTE THERE ARE NONE AT THIS TIME!
d[sample %in% Esamp & organism=='Enterovirus.D.68', present := FALSE]
d[sample %in% Esamp & organism=='Enterovirus.D.68', ct      := NA]

# entero + EVD68 (lower Ct (anything remaining after the above exclusions)) == EV D68
d[organism=='EV_pan'&present==TRUE,        tmp_EE68 := 'EPAN']
d[organism=='Enterovirus.D.68'&present==TRUE, tmp_EE68 := 'Enterovirus.D.68']
d[, tmp_EE682 := all(c("EPAN","Enterovirus.D.68") %in% tmp_EE68), by = sample]
d[tmp_EE682==TRUE & organism=='EV_pan', present := FALSE]
d[tmp_EE682==TRUE & organism=='EV_pan', ct      := NA]

# Eliminate IAV as we have nested targets
d <- d[organism != 'IAV']
d[organism=='Influenza.A.H1N1', organism := 'IAV..H1N1']
d[organism=='Influenza.A.H3N2', organism := 'IAV..H3N2']

# These seem to be a few hundred duplicated RSV samples, drop those as they are covered by A or B nested targets
d <- d[organism != 'RSV']

# RSV A and B can have some cross-reactivity. 
#   The lower Ct value for RSV should be used because probe sets overlap somewhat 
#   and will cross react at high cycles. So if you see a screaming low Ct RSVB, with a high Ct RSVA, 
#   assume cross reactivity of the probe sets not co-infection. Use a delta threshold of 10 (Lea)
rsv_thresh <- 10
tmp <- d[organism %in% c('RSV.A', 'RSV.B') & !is.na(ct)]
tmp[, bothRSVs := sum(present), by=sample]
tmp <- tmp[bothRSVs==2][order(sample,organism)] # under 150 cases
tmp[,i:=1:.N,by=sample]
tmp[,ctdiff:=as.numeric(ct)-as.numeric(data.table::shift(ct,1,type='lag')),by=sample]
falseRSVA <- tmp[ctdiff< -1*rsv_thresh]$sample
falseRSVB <- tmp[ctdiff>    rsv_thresh]$sample
d[sample%in%falseRSVA & organism == 'RSV.A', present := FALSE]
d[sample%in%falseRSVA & organism == 'RSV.A', ct := NA]
d[sample%in%falseRSVB & organism == 'RSV.B', present := FALSE]
d[sample%in%falseRSVB & organism == 'RSV.B', ct := NA]

## Finally, clean up some organism names
# Dropping really rare (<50 positives) or ones that had poor temporal coverage
d <- d[!organism %in% c('Measles','Mumps','Chlamydophila_pneumoniae','Human_bocavirus','Human_parechovirus',
                        'Bordetella_pertussis','Mycoplasma_pneumoniae')]

# rename
d[organism == 'Coronavirus_2019',         organism := 'SARS_CoV_2']
d[organism == 'Human_coronavirus.2019',   organism := 'SARS_CoV_2']
d[organism == 'Enterovirus.D.68',         organism := 'EV..D68']
d[organism == 'EV_pan',                   organism := 'EV']
d[organism == 'Human_metapneumovirus',    organism := 'hMPV']
d[organism == 'Human_bocavirus',          organism := 'HBoV']
d[organism == 'Human_parechovirus',       organism := 'HPeV']
d[organism == 'Streptococcus_pneumoniae', organism := 'SPn']
d[organism == 'RSV.A',                    organism := 'RSV..A']
d[organism == 'RSV.B',                    organism := 'RSV..B']


## COLLAPSE TO THE ORGANISM LEVEL
# Note: collapse to True if any present, and use min Ct
table(d$organism)
d <- d[,.(present = any(present==TRUE,na.rm=TRUE),ct=min(ct,na.rm=TRUE)),
           by = c(finalvars[!finalvars%in%c('target','ct','present')])]
d[present==FALSE, ct := NA]
d[ct == Inf,      ct := NA] # A few weird rows with present==TRUE and no Ct Value
table(d$organism)
table(d[present==TRUE]$organism)


d[, organism := factor(organism, levels = sort(unique(d$organism)))]


## label week for easy aggregation later
d[, week := ymd(paste0(year(date),"-01-01" )) + weeks(week(date) - 1 )]

# add one hanger day back onto the previous week so its not on its own
d[week == '2020-12-31', week := ymd('2020-12-24')]
d[week == '2019-12-31', week := ymd('2019-12-24')]

# Some 'weirdos' come flu positive are likely flumist, set these to false:
# https://docs.google.com/spreadsheets/d/1oEyALYuowLHynnDH8CIq1Ov1t_z-UwvO80PBGTtkqko/edit#gid=0
# I think these match under investigator ID (grep sample for it)
weirdo_ids <- c('bc6bdb11','fb5c850f','75c562de','c0ca3772','d4db8dcf',
                '9a94c6d7','cf814209','b5e87ec2',
                '2c14d08d','32b34f99','697e2074','821f0ee8','86ee8dfc','cdf8f495')
d[grepl(paste(weirdo_ids,collapse='|'),sample) & grepl('IAV|IBV|ICV',organism), present := FALSE]


# add a few helpful variables
d[, site2 := site_type]
d[!site2 %in% c('retrospective','clinic','SCAN','NA'), site2 := 'community']

# aggregate age bins
d[age<1,             agebin := '<1']
d[age>=1 & age<=5,   agebin := '1-5']
d[age>5  & age<=16,  agebin := '6-16']
d[age>16 & age<=45,  agebin := '17-45']
d[age>45 & age<=65,  agebin := '46-65']
d[age>65,            agebin := '66+']
d[, agebin := factor(agebin, levels=c('<1','1-5','6-16','17-45','46-65','66+'))]

# secondary age grouping similar to CDC
d[age<1,             agebin2 := '<1']
d[age>=1 & age<=4,   agebin2 := '1-4']
d[age>4  & age<=17,  agebin2 := '5-17']
d[age>17 & age<=49,  agebin2 := '18-49']
d[age>49 & age<=64,  agebin2 := '50-64']
d[age>64,            agebin2 := '65+']
d[, agebin2 := factor(agebin2, levels=c('<1','1-4','5-17','18-49','50-64','65+'))]

# set site back to missing if its so
d[site2=='NA', site2 := NA]

# We will track SCAN as community
d[ site2 == 'SCAN']$site2 <- 'community'

# Kaiser is a separate type (outpatient clinic, but retrospective)
d[ site == 'KaiserPermanente']$site2 <- 'clinic retrospective'
d[ site2 == 'clinic']$site2 <- 'clinic kiosk'

# year.month
d[, ym := paste0(year(date),'-',month(date,label=TRUE))]

# sort columns 
d <- d[order(date,sample,organism)]

# save the final data set to be used for analysis
saveRDS(d, 'model_data.RDS')

