## ----------------------------------------------------------------------------------------
## ROY BURSTEIN
## 2021
## Full Replication code for the paper
## Note, this script follows cache_id3c.R and prep_mp_data.R scripts, which pull 
##  and prep the data, respectively
## Working Draft: https://docs.google.com/document/d/1Xo3OpqmdSFkqNtSP0s0W9Bia9WeIlDVQ9q4EcillXtw/edit
## ----------------------------------------------------------------------------------------

# RUN FAKE DATA?
# SET TO TRUE IF YOU DO NOT HAVE PERMISSION TO ACCESS ID3C SFS DATABASE
RUNFAKE <- TRUE


## ----------------------------------------------------------------------------
## SCRIPT SET UP

## Libraries
library(magrittr)
library(data.table)
library(ggplot2)
library(extrafont)
library(lubridate)
library(jtools)
library(gridExtra)
library(knitr)

# Several Color palettes to use 
# Carto
palprism <- c('#5F4690','#1D6996','#38A6A5','#0F8554','#73AF48','#EDAD08',
              '#E17C05','#CC503E','#94346E','#6F4070','#994E95','#666666')
palsafe  <- c('#88CCEE','#CC6677','#DDCC77','#117733','#332288','#AA4499',
              '#44AA99','#999933','#882255','#661100','#6699CC','#888888')   
palantiq <- c('#855C75','#D9AF6B','#AF6458','#736F4C','#526A83','#625377',
              '#68855C','#9C9C5E','#A06177','#8C785D','#467378','#7C7C7C')
palvivid <- c('#E58606','#5D69B1','#52BCA3','#99C945','#CC61B0','#24796C',
              '#DAA51B','#2F8AC4','#764E9F','#ED645A','#CC3A8E','#A5AA99')
# other
newcols <- rev(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                 '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
                 '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462',
                 '#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f'))
clz <- c('#58EFECff', '#75D2DAff', '#92B4C7ff', '#AE97B5ff', '#CB79A2ff', '#E85C90ff')


# read in data (DL'd in cache_id3c.R and prepped in the prep_mp_data.R script)
setwd(<<WORKING DIRECTORY PATH HERE>> )
dOrig <- readRDS(paste0(ifelse(RUNFAKE==TRUE,'fake_',''),'model_data_final.RDS'))

# should see consistency among included organisms
table(dOrig$organism)

# Track exclusions as they are made through data cleaning steps. 
# Track number of samples remaining each step
exclusions <- c('Original'=length(unique(dOrig$sample)))
d     <- copy(dOrig)

## ----------------------------------------------------------------------------
## BASIC DATA PREP

# keep only samples with present pathogen
d <- d[present==TRUE]
exclusions <- append(exclusions, c('Present Pathogens Only' = length(unique(d$sample))))

# drop cohorts and repeats
d <- d[!site %in% c('StMartins','DESC','HutchKids','UWDaycare','WestCampusChildcareCenter')]
exclusions <- append(exclusions, c('Drop Cohorts' = length(unique(d$sample))))


# drop anyone with repeat measures (keep only first positive)
tmp <- unique(d[present==TRUE,c('sample','individual','date')])[order(individual,date)]
tmp[,order:=1:.N,by=individual]
firstsamples <- unique(c(tmp[order==1]$sample,tmp[is.na(individual)]$sample)) # assume NAs are individuals.
d <- d[sample %in% firstsamples] 
exclusions <- append(exclusions, c('Keep only first positives' = length(unique(d$sample))))


# Many samples were tested for covid only so not relevant to the current study, 
#  make sure we only keep samples tested for multiple pathogens
dOrig[,ntestsinsample:=.N,by=sample]
table(dOrig$ntestsinsample) 
d <- d[sample %in% dOrig[ntestsinsample>=16]$sample] # only those tested for full panel (some 16, pre-COVID)
exclusions <- append(exclusions, c('Dropping those not tested for full array' = length(unique(d$sample))))


## CREATE SOME VARIABLES USED IN THE ANALYSIS


# make a few useful coinfection variables we can use later
#  coinf indicator, number of coinfections, as nice factor, 
#  number of viral coinfections (not used in analysis)
d[, coinf   := sum(present,na.rm=T)>1, by = sample]
d[, coinfN  := sum(present,na.rm=T),   by = sample]
d[coinfN>4, coinfN := 4]
d[, coinfNf  := factor(coinfN,labels = c('Mono','2','3','4+'))]

d[, presentV := present==TRUE & organism != 'SPn']
d[, coinfV   := sum(presentV,na.rm=T)>1, by=sample]
d[, coinfNV  := sum(presentV,na.rm=T),   by=sample]
d[present==TRUE & organism == 'SPn', coinfNV := coinfNV+1]
d[coinfNV>4, coinfNV := 4]    
d[, coinfNVf := factor(coinfNV,
                       levels = c(1,2,3,4),
                       labels = c('Mono','2','3','4+'))]
                                   
# year-month as ordered factor
d[, ym := factor(ym, levels = unique(d[order(date)]$ym))]

# keep only needed variables 
d <- d[,c( 'sample', 'organism', 'agebin2', 'site2', 'ym', 'site','device', 
           'present', 'coinf', 'coinfN', 'coinfNf', 'coinfNV', 'ct', 'coinfV')]
d <- na.omit(d, cols = c('ct','agebin2', 'site2', 'ym'))
exclusions <- append(exclusions, c('Dropped those missing metadata' = length(unique(d$sample))))

# dall is the complete dataset (including negatives) for all included samples
dall <- copy(dOrig)[sample %in% unique(d$sample)]

table(d$device,d$organism,useNA='always') # about 10% of covid positives are on OpenArray
table(d[device=='OpenArray'&organism=='SARS_CoV_2']$site) # note research only samples from retros


## ----------------------------------------------------------------------------
## MAKE TABLE 1: SAMPLE CHARACTERISTICS BY PATHOGEN

# sumt() function formats a nice table1
sumt <- function(d, byvar=NULL, reshape=TRUE){
  if(!is.null(byvar))
    byvars <- c('organism',byvar)
  else 
    byvars <- 'organism'
  ds <- d[,.(N=.N,crt_mean=round(mean(ct),1),crt_sd=round(sd(ct),1)),by=byvars]
  ds <- ds[order(organism)]
  ds[, val := sprintf('%s (%s, %s)',crt_mean,crt_sd,N)]
  ds <- ds[,c(byvars,'val'),with=F]

  if(!is.null(byvar)){
    if(reshape==TRUE){
      ds <- reshape(ds, direction = 'wide', idvar = byvar, timevar = c('organism'))
    }
    setnames(ds,byvar,'var')
    ds <- ds[order(var)]
    ds[, var := paste(byvar,'-',var)]
  } else {
    ds$var <- 'Overall'
    if(reshape==TRUE){
      ds <- reshape(ds, direction = 'wide', idvar = 'var', timevar = c('organism'))
    }
  }
  return(ds)
}

# pull out summaries
dsumA <- sumt(d,'agebin2')  
dsumS <- sumt(d,'site2')
dsumC <- sumt(d,'coinf')
dsumO <- sumt(d, NULL)

# bind them as one large table
table1 <- rbindlist(list(dsumO,dsumC,dsumA,dsumS))
setnames(table1,names(table1),gsub('val.','',names(table1)))

# Save output
write.csv(table1, 'table1.csv')

## save a copy of the analysis dataset
write.csv(d[,c('sample','organism','agebin2','site2','ym','device','coinf','ct')],
          'interaction_analysis_clean_dataset.csv')



## --------------------------------------------------------------------------
## MAIN ANALYSIS 

## Set up analysis dataset
# wide by pathogen co-presence
tmp <- copy(d)
for(v in unique(tmp$organism)){
  tmp[, vy := present==TRUE & organism==v]
  tmp[, paste0('inf_',v) := as.numeric(sum(vy,na.rm=T)>0), by = sample]
}

# dataset has monoinfections plus each coinf pair. 
combl <- data.table()
for(V in unique(tmp$organism)){
  tmpl <- tmp[organism==V] # & coinfN <= 2]  ## <- note, uncomment this to do sensitivity analysis for SI Fig 1
  tmpl <- melt(tmpl[,c('ct','sample','coinfN','agebin2','site2','ym','device',paste0('inf_',unique(tmp$organism))),with=FALSE],
               id.vars = c('sample','ct','coinfN','agebin2','site2','device','ym'),
               variable.name = 'coinfector')[order(sample)][value>0]
  tmpl[, coinfector := gsub('inf_','',coinfector)]
  tmpl[, value := NULL]
  tmpl[, single_infection := coinfector == V & coinfN  == 1]
  tmpl <- tmpl[!(single_infection==FALSE & coinfector==V)] # remove coinf marked as V of interest
  tmpl[single_infection==TRUE, coinfector := 'x MONO-INFECTION'] 
  tmpl[, organism := V]
  combl <- rbind(combl,tmpl)
}
combl[, coinfected := single_infection == FALSE]


## Make a matrix of above or below expected CRT for pair when co-infected. 
#   The comparison is CT of Y when coinf with X to CT of Y when singly infected
#   Every pair (v, v2), is looped through here and a regression is run. 
#   Also do a 'marginal' version.

# make objects to store results and a function to do it in a standarized way
storeres  <- function(Y,X,RES,M,DAT)
  rbind(RES,
     data.table(
       Y            = Y,
       X            = X,
       intercept    = coef(summary(M))[1,1],
       diff         = coef(summary(M))[2,1], 
       diff_se      = coef(summary(M))[2,2],
       diff_pval    = coef(summary(M))[2,4],
       N_singleinf  = sum(DAT$single_infection==TRUE),
       N_coinf      = sum(DAT$single_infection==FALSE)
     ))

# adjusted formula for regressions to run (SARS-CoV-2 Needs an additional adjustment for PCR device)
getformula <- function(v){
  if(v=='SARS_CoV_2'){
    formla <- as.formula('ct ~ factor(coinfected)+factor(site2)+ym+factor(agebin2)+factor(device)')
  } else if(v=='(All Viruses)-->'){
    formla <- as.formula('ct ~ factor(coinfected)+factor(site2)+ym+factor(agebin2)+factor(organism)')
  } else if(v=='-->(All Viruses)'){
    formla <- as.formula('ct ~ factor(coinfected)+factor(site2)+ym+factor(agebin2)+factor(organism)')
  } else {
    formla <- as.formula('ct ~ factor(coinfected)+factor(site2)+ym+factor(agebin2)')
  }
} 

res         <- data.table() # Pathogen-Pathogen interaction Results
resmar      <- data.table() # Marginal interaction results store
modlist     <- list()       # store for model objects
modlistmar  <- list()


# loop over each pathogen
for(v in unique(combl$organism)){
  message(v)
  
  # keep mono and co-inf data for each
  tmp1 <- combl[organism == v & single_infection==TRUE]
  tmp2 <- combl[organism == v & single_infection==FALSE] 
  

  ## Marginal Interaction Analyses 
  # (ANY VIRUS --> V) 
  tmp     <- rbind(tmp1, tmp2[coinfector!='SPn'])
  # should remove duplicates in cases of >1 coinfector, collapse and add
  tmp <- unique(tmp[,names(tmp)[names(tmp)!='coinfector'],with=F])
  m_noadj <- lm(ct ~ factor(coinfected), data = tmp) 
  m       <- lm(getformula(v), data = tmp) 
  
  modlistmar[[paste0('(All Viruses)','-->',v)]] <- list('Unadjusted' = m_noadj, 'Adjusted' = m, 'data' = tmp)
  resmar <- storeres(v,'(All Viruses)',resmar,m,tmp)
  
  # (V --> ANY VIRUS) 
  # When this virus is coinfector, how does it impact other viruses broadly?
  tmp     <- rbind(combl[!organism%in%c(v,'SPn')&single_infection==TRUE], combl[coinfector==v & !organism%in%c(v,'SPn') &single_infection==FALSE])
  m_noadj <- lm(ct ~ factor(coinfected), data = tmp) 
  m       <- lm(getformula('-->(All Viruses)'), data = tmp) 
  
  modlistmar[[paste0(v,'-->','(All Viruses)')]] <- list('Unadjusted' = m_noadj, 'Adjusted' = m, 'data' = tmp)
  resmar <- storeres('(All Viruses)',v,resmar,m,tmp)


  
  ## Pathogen - Pathogen specific interaction
  for(v2 in unique(combl$organism)){
    
    message(paste0(v,'  -----> ',v2))
    tmp <- rbind(tmp1, tmp2[coinfector==v2])
    
    if(v2 %in% unique(tmp2$coinfector) & nrow(tmp)>20){ # under 20 total obs itll usually break, dont bother running
      
      m_noadj <- lm(ct ~ factor(coinfected), data = tmp) 
      m       <- lm(getformula(v), data = tmp) 
      
      modlist[[paste0(v2,'-->',v)]] <- list('Unadjusted' = m_noadj, 'Adjusted' = m, 'data' = tmp)
      res <- storeres(v,v2,res,m,tmp)
      
    } else { # if(v!=v2) {
      res <- rbind(res,
                   data.table(
                     Y=v,
                     X=v2,
                     intercept=NA,diff=NA,diff_se=NA,diff_pval=NA,
                     N_singleinf=sum(tmp$single_infection==TRUE),
                     N_coinf=sum(tmp$single_infection==FALSE)
                   ))
    }
  }
}



## --------------------------------------------------------------------------
## FIGURE 3

# add res-marginal to res so it plots the margins.
res <- rbind(res, resmar)


## graphing parameters
res[, ss_toosmall := N_singleinf < 10 | N_coinf < 10]
mx <- ceiling(max(abs(range(res[ss_toosmall==FALSE]$diff))))
PVALTHRESHOLD <- 0.01 # setting significance level here to 0.01

# order pathogens for plotting
lvls <- c("(All Viruses)"     ,'', "EV..D68"       , "SARS_CoV_2"    ,  "EV"       ,    
          "hPIV..1_2",
          "ICV"           , "hMPV"          ,  "IAV..H3N2"  ,    "IAV..H1N1",
          "IBV"           , "hPIV..3_4"     ,  "RSV..B"     ,    "RSV..A"   ,    
          "CoV..HKU1_NL63", "CoV..229E_OC43",  "AdV"        ,    "RV"       , 
          "SPn")
res[, X:= factor(X,levels=lvls)]
res[, Y:= factor(Y,levels=lvls)]

res <- rbind(res,expand.grid(X=lvls,Y='',ss_toosmall=F),fill=T)
res <- rbind(res,expand.grid(Y=lvls,X='',ss_toosmall=F),fill=T)

eq <- data.table(X=factor(unique(res$X),levels=lvls),
                 Y=factor(unique(res$X),levels=lvls))

# set sig threshold
res[,pthresh := diff_pval<PVALTHRESHOLD]


## SAVE THE FIGURE
png('./figs/fig 3 with Marginals.png',height=800,width=1100)

ggplot(res[ss_toosmall==FALSE],aes(X,Y)) +
  theme_classic() + 
  geom_point(pch=4, size = 3, color ='grey', alpha = .2) +
  geom_tile(aes(fill=diff),color='black',lwd=.01) +
  #geom_point(data=eq, pch='o') +
  geom_text(aes(label=N_coinf), 
            position = position_nudge(x=0.25,y=-0.22),
            pch = 1, size = 3.5, alpha = 0.25) + 
  geom_text(aes(label=round(diff)), 
            family = "Garamond",
            position = position_nudge(x=-0.2,y=0),
            pch = 1, size = 6.25, alpha = 0.25) + 
  geom_point(aes(shape=pthresh), size = 13, alpha = .65,
             position = position_nudge(x=0.27,y=0.27)) + 
  scale_shape_manual(values=c('','*')) +
  scale_fill_gradientn(
    colours = c("#AF4034","#F16B6F", "white", "#00b9f1","#005f6b"),  
    limits  = c(-mx,mx), breaks = c(-mx,-0.5*mx,0,0.5*mx,mx),
    name    = "\u0394 Ct") +
  xlab('When coinfected with:') +
  ylab('Average difference in Ct of:') +
  guides(shape = 'none') +
  theme(text = element_text(size=25 , family = "Garamond"), # , face = 'bold'
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=16),
        axis.text.x = element_text(angle=315,hjust=0), 
        legend.position = 'right')
  
dev.off()                                                                                                                                                                                                                                                                                                             


## --------------------------------------------------------------------------
## FIGURE -- MARGINALS

resmar[,pthresh := diff_pval<PVALTHRESHOLD]

ggplot(resmar,aes(X,Y)) +
  theme_classic() + 
  geom_point(pch=4, size = 3, color ='grey', alpha = .2) +
  geom_tile(aes(fill=diff),color='black',lwd=.01) +
  #geom_point(data=eq, pch='o') +
  geom_text(aes(label=N_coinf), 
            position = position_nudge(x=0.25,y=-0.22),
            pch = 1, size = 3.5, alpha = 0.25) + 
  geom_text(aes(label=round(diff)), 
            family = "Garamond",
            position = position_nudge(x=-0.2,y=0),
            pch = 1, size = 6.25, alpha = 0.25) + 
  geom_point(aes(shape=pthresh), size = 13, alpha = .65,
             position = position_nudge(x=0.27,y=0.27)) + 
  scale_shape_manual(values=c('','*')) +
  scale_fill_gradientn(
    colours = c("#AF4034","#F16B6F", "white", "#00b9f1","#005f6b"),  
    limits  = c(-mx,mx), breaks = c(-mx,-0.5*mx,0,0.5*mx,mx),
    name    = "\u0394 Crt") +
  xlab('When coinfected with:') +
  ylab('Average difference in Crt of:') +
  guides(shape = 'none') +
  theme(text = element_text(size=25 , family = "Garamond", face = 'bold'),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=16),
        axis.text.x = element_text(angle=315,hjust=0), 
        legend.position = 'right')


## --------------------------------------------------------------------------
## NUMBER PULL -- ALL THE NUMBERS CALLED OUT IN THE MANUSCRIPT

## Total pathogens in the analysis
length(unique(d$organism))

## Total pairs in the analysis
nrow(res)
sum(!is.na(res[ss_toosmall==FALSE]$diff))

## Date range
min(d$date)
max(d$date)

# original samples with at least one positive
dOrig[, N := .N, by = sample]
tmp <- dOrig[N>=16,.(anypos=sum(present,na.rm=T)>0),by=.(N,sample)]
nrow(tmp); sum(tmp$anypos); mean(tmp$anypos)

# total samples
length(unique(d$sample))
nrow(d)

# by pathogen
dall[,.(P=sum(present),pct=round(mean(present)*100,1)),by=organism][order(-pct)]

# co-infections
tmp <- unique(d[,c('sample','coinfN')])[,.(P=.N),by=coinfN]
tmp[,pct:=P/sum(P)][order(coinfN)]


## Number of samples
N=length(unique(dall$sample))
exclusions

# samples from sites, specific queries
table(unique(d[,c('sample','site2')])$site2)
table(unique(d[,c('sample','site2')])$site2)/sum(table(unique(d[,c('sample','site2')])$site2))
table(unique(d[,c('sample','site2')])$site2)

table(unique(d[,c('sample','site','site2')])[site2=='retrospective']$site)
table(unique(d[,c('sample','site','site2')])[site2=='community']$site)
table(unique(d[,c('sample','site','site2')])[grepl('clinic',site2)]$site)
table(unique(d[,c('sample','site','site2')])$site)

table((d[site=='SCAN']$organism))

## results for main analysis called out in results:
head(res[order(-N_coinf)],10) # most coinfections
nrow(res); length(unique(res$Y))*(length(unique(res$X))-1) # possible 
nrow(res[ss_toosmall==FALSE]) # how many with > threshold (10) observations
head(res[ss_toosmall==FALSE&diff_pval < PVALTHRESHOLD][order(-diff)],10) # largest effect
head(res[ss_toosmall==FALSE&diff_pval < PVALTHRESHOLD][order(diff)], 10) # largest effect
nrow(res[ss_toosmall==FALSE&diff_pval < PVALTHRESHOLD]) # signif pairs
res[ss_toosmall==FALSE&diff_pval <  PVALTHRESHOLD&diff<0][order(diff)] # all  synergy SPn
res[ss_toosmall==FALSE&diff_pval >= PVALTHRESHOLD&Y=='SPn'][order(diff)] # SPn viruses not sig
res[ss_toosmall==FALSE&diff_pval <  PVALTHRESHOLD&X=='SPn'][order(diff)] # Spn assoc with sig inhibition in these
res[ss_toosmall==FALSE&diff_pval <  PVALTHRESHOLD&diff<0&Y!='SPn'&X!='SPn'] # no V-V synergistics
res[ss_toosmall==FALSE&diff_pval <  PVALTHRESHOLD&diff>0&Y!='SPn'&X!='SPn'] # V-V Interferers
res[ss_toosmall==FALSE&Y!='SPn'&X!='SPn'] # total assessed V-V interactions
res[ss_toosmall==FALSE&diff_pval<PVALTHRESHOLD&diff>0& X %in% c('IAV..H1N1','IAV..H3N2','IBV') & Y!='SPn'] # sig IAV, IBV inx
res[ss_toosmall==FALSE& X %in% c('IAV..H1N1','IAV..H3N2','IBV') & Y!='SPn'] # total IAV IBV inx
res[ss_toosmall==FALSE&diff_pval<PVALTHRESHOLD&diff>0& Y %in% c('IAV..H1N1','IAV..H3N2','IBV') & X!='SPn'] # sig IAV, IBV inx
res[ss_toosmall==FALSE& Y %in% c('IAV..H1N1','IAV..H3N2','IBV')& X!='SPn'] # total IAV IBV inx
#summary(modlist[['hMPV-->IAV..H3N2']]$Adjusted) # notable difference

res[ss_toosmall==FALSE & (X == 'SARS_CoV_2' | Y == 'SARS_CoV_2')] 
res[ss_toosmall==TRUE &  (X == 'SARS_CoV_2' | Y == 'SARS_CoV_2')] 

# of XX positive samples from children under 5, XX (XX%) were positive for multiple pathogens, 
# while of XX samples from adults 17-45, only XX (XX%) were coinfections. 
(x <- unique(d[agebin2 %in% c('<1','1-4'),c('sample','coinfN'),with=F])[,.(p=.N),by=coinfN])
sum(x$p[2:4])/sum(x$p)
sum(x$p)
sum(x$p[2:4])

(x <- unique(d[agebin2 %in% c('18-49'),c('sample','coinfN'),with=F])[,.(p=.N),by=coinfN])
sum(x$p[2:4])/sum(x$p)
sum(x$p)
sum(x$p[2:4])


## total cases
nrow(d)
length(unique(d$sample))






## OTHER FIGURES
## -----------------------------------------------------------------------
## FIGURE 1: age distribution by pathogen

# breack apart ct and crt 

fig1d <- d[!(organism == 'SARS_CoV_2' & device == 'OpenArray')]

nrow(d[(organism == 'SARS_CoV_2' & device == 'OpenArray')])
nrow(d[(organism == 'SARS_CoV_2' & device != 'OpenArray')])
nrow(d[(organism == 'SARS_CoV_2' & device == 'OpenArray')])/
  nrow(d[(organism == 'SARS_CoV_2')])



png('./figs/fig1.png',height=800,width=1200)
ggplot(fig1d) +
  geom_jitter(aes(ct,agebin2,color=agebin2), alpha = .24,size = 1.2,
              position = position_jitter(height = 0.5, width = 0.5))+
  facet_wrap(.~organism, scales = 'free_x', nrow=4) +
  geom_point(data=fig1d[,.(m=mean(ct,na.rm=T)),by=.(agebin2,organism)],
             aes(m,agebin2,color=agebin2),size=8)+
  geom_point(data=fig1d[,.(m=mean(ct,na.rm=T)),by=.(agebin2,organism)],
             aes(m,agebin2,color=agebin2),size=8,
             shape = 1,colour = "black") +
  scale_fill_manual(values = clz) + 
  scale_color_manual(values = clz) + 
  theme_classic() + 
  xlab('Ct') + ylab('Age') +
  scale_x_continuous(breaks = c(5,10,15,20,25,30,35), limits = c(0,NA)) +
  theme(text = element_text(size=33 , family = "Garamond"),
        legend.position = 'none',
        panel.grid.major.x = element_line( size=.4, color="black" ) )
dev.off()

# distributions over co-variables
vars <- 'coinfNVf' # c('coinfNf','coinf','site2','sex')
for(var in vars){
  png(sprintf('./figs/%s_ct_distribution.png',var),height=800,width=1200)
  plot(
    ggplot(d) +
      geom_jitter(aes(ct,factor(get(var)),color=factor(get(var))), alpha = .07,size = 1.2,
                  position = position_jitter(height = 0.3, width = 0.5))+
      facet_wrap(.~organism) +
      geom_point(data=d[,.(m=mean(ct,na.rm=T)),by=c(var,'organism')],
                 aes(m,factor(get(var)),color=factor(get(var))),size=8)+
      geom_point(data=d[,.(m=mean(ct,na.rm=T)),by=c(var,'organism')],
                 aes(m,factor(get(var)),color=factor(get(var))),size=8,
                 shape = 1,colour = "black") +
      scale_fill_manual(values  = clz[2:6]) + 
      scale_color_manual(values = clz[2:6]) + 
      theme_classic() + 
      xlab('Ct') + ylab('Number of Viral coinfections') + #ylab(var) +
      scale_x_continuous(breaks = c(5,10,15,20,25,30)) +
      theme(text = element_text(size=30 , family = "Garamond"),
            legend.position = 'none',
            panel.grid.major.x = element_line( size=.4, color="black" ) )
  )
  dev.off()
}


## --------------------------------------------------------------------------
## FIGURE 2: SAMPLES OVER TIME


# get total samples and proprotion positive per organism
dd <- unique(copy(dall)[present==TRUE,c('sample','organism','week')])[,.(N=.N),by=.(organism,week)] %>%
  na.omit()
dd <- dd[!organism%in% dd[,.(N=sum(N)),by=organism]] 
dd<-merge(dd,
          unique(expand.grid(unique(dd[,c('week','organism')]))),
          by = c('week','organism'), all = TRUE)
dd[is.na(N), N:= 0]
dd[,pct := N/sum(N),by=week]

total <- copy(dOrig)[present==TRUE&duplicated(sample)==FALSE,.(N = .N), by = week]
ratio <- max(dd[,.(N=sum(N)),by=week]$N,na.rm=T)/max(total$N,na.rm=T)


lab <- data.table(
  x=c(ymd('2021-01-30'),ymd('2019-01-30'),
      ymd('2020-01-30')), y=c(750,750,750),lab=c(2021,2019,2020))

lab2 <- data.table(
  x=c(ymd('2019-04-15'), ymd('2019-06-15'), ymd('2020-11-15'),
      ymd('2019-03-01'), ymd('2019-12-01'), ymd('2019-01-30'),
      ymd('2019-05-15'), ymd('2019-05-15'), ymd('2020-01-30'),
      ymd('2019-02-20'), ymd('2019-01-01'), ymd('2019-12-15')),
  y=c(.125,              .4,                .4,
      .5,                .55,               .62,
      .98,               .65,                .92,
      .9,               .33,               .4),
  lab=c('SPn',           'RV',             'SARS-CoV-2',
        'IAV H3N2',      'IBV',            'IAV H1N1',
        'AdV',           'hPIV 3-4',       'CoV HKu1, NL63',
        'CoV 229E, OC43','RSV B',          'RSV A') 
)


g1 <- ggplot(dd[week<'2021-09-01']) + 
  theme_classic() +
  ylab('Positive Weekly Samples') + 
  xlab('') +
  ylim(0,800)+
  theme(legend.position = 'none', axis.text.x = element_blank()) +
  scale_fill_manual(values = rev(c(palprism,palvivid))) +
  geom_bar(aes(week,N,fill=organism), stat='identity', width = 7.2) +
  theme(text = element_text(size=45 , family = "Garamond")) +
  geom_vline(xintercept = c(ymd('2020-12-28'),ymd('2018-12-28'),
                            ymd('2019-12-28')), color='grey', lwd=1.2) +
  geom_text(data=lab,aes(x,y,label=lab), family = 'Garamond', color = 'grey', size = 10)

g2 <- ggplot(dd[week<'2021-09-01']) + 
  theme_classic() +
  scale_y_continuous(labels=scales::percent, name = '% of Positive Weekly Samples') + 
  scale_fill_manual(values = rev(c(palprism,palvivid))) +
  geom_area(aes(week,pct,fill=organism), stat='identity', width = 7.2) + # color = 'white', lwd=.01,
  scale_x_date(breaks='1 month', name = '',labels=scales::date_format("%b")) +
  labs(fill='') +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle=45,  hjust=1),
        text = element_text(size=45 , family = "Garamond"))  +
  geom_vline(xintercept = c(ymd('2021-01-01'),ymd('2019-01-01'),
                            ymd('2020-01-01')), color='grey', lwd=1.2) +
  geom_text(data=lab2,aes(x,y,label=lab), family = 'Garamond', color = 'black', size = 8)

png('./figs/fig2 fin.png',height=1500,width=1700)
egg::ggarrange(g1,g2)
dev.off()



## -----------------------------------------------------------------------
## FIGURE 4: PREVALENCE VERSUS CT (SI FIGURE)

tmp <- dOrig[N>=16,.(pp = mean(present), N=sum(present), ct = mean(ct,na.rm=T)), by = .(agebin2,ym,organism)]
tmp <- na.omit(tmp)

png('./figs/fig 4.png',height=700,width=1000)
ggplot(tmp[!is.na(agebin2) & N >= 10]) +
  geom_smooth(aes(pp,ct),se=T,method='lm',color='grey',lwd=1.2, alpha=0.3) +
  geom_point(aes(pp,ct,color=agebin2,size=N), alpha = .75) + 
  facet_wrap(.~organism, scales= 'free') +
  scale_color_manual(values=clz) +
  theme_classic() + 
  labs(color='Age', size='Positive\nSamples') +
  ylab('Average Crt/Ct') +
  scale_x_continuous(labels=scales::percent, name = 'Proportion Positive') +
  theme(text = element_text(size=32 , family = "Garamond"),
        legend.key.size = unit(1, 'cm'),
        axis.text.x  = element_text(size=20),
        legend.text = element_text(size=16),
        legend.position = 'bottom')
dev.off()

# is prevalence related outside of age?
tmp <- dOrig[N>=16,.(pp = mean(present), N=sum(present), ct = mean(ct,na.rm=T)), by = .(agebin2,ym,organism,site2)]
tmp <- na.omit(tmp)

rs <- data.table()
for(o in unique(d$organism)){
  rs <- rbind(rs,data.table(
    org = o,
    B   = summary(lm(pp*100~ct+agebin2+site2,tmp[organism==o]))$coefficients[2,c(1)],
    z   = summary(lm(pp*100~ct+agebin2+site2,tmp[organism==o]))$coefficients[2,c(3)]
  ))
}
rs <- rs[order(B,z)] 
rs[abs(z)>2]
rs


## -----------------------------------------------------------------------
##  SI Figure: distribution of coinfections over ages. 

tmp <- copy(dOrig)

tmp[, coinfN := sum(present,na.rm=T), by = sample]
tmp <- unique(tmp[,c('site2','site','N','coinfN','sample','agebin2')])
tmp[coinfN>3, coinfN := 3]
tmp[,coinfN := factor(coinfN, levels = 3:0, labels = rev(c('0',"1","2",'3+')))]

tmp[ site2 == 'SCAN']$site2 <- 'community'
tmp[ site == 'KaiserPermanente']$site2 <- 'clinic retrospective'
tmp[ site2 == 'clinic']$site2 <- 'clinic kiosk'

tmp  <- tmp[N>=16, .(n=.N), by = .(agebin2, site2, coinfN)] %>% na.omit()
tmp2 <- tmp[N>=16, .(n=sum(n)), by = .(agebin2, coinfN)] %>% na.omit()
tmp[, pct := n/sum(n), by = .(agebin2,site2)]
tmp2[, pct := n/sum(n), by = .(agebin2)]
tmp2[, site2 := 'Overall']
tmp <- rbind(tmp,tmp2)
tmp[, site2 := factor(site2, 
                      levels = c("clinic retrospective", "clinic kiosk", 
                                 "community", "retrospective", "Overall" ),
                      labels = c('Clinic Resid.', 'Clinic Kiosk', 'Community',
                                 'Hosp. Resid.', 'TOTAL'))]
tmp[,ss := paste0('(N=',prettyNum(sum(n),big.mark=','),')'),by=.(site2,agebin2)]

png('./figs/si fig 2.png',height=1000,width=1300)
ggplot(tmp) +
  geom_bar(aes(x=0,y=pct,fill=coinfN), stat = 'identity') +
  scale_y_continuous(labels = scales::percent, name = 'Percent of Samples') +
  scale_fill_manual(values=c('#FFBC42','#D81159','#8F2D56','#218380'),
                    name='Organisms\nDetected\nPer\nSample') +
  theme_classic() + xlab('') +
  theme(text = element_text(size=30 , family = "Garamond"),
        axis.text.x = element_blank()) +
  facet_grid(site2~agebin2) +
  geom_text(aes(x=0,y=.08,label=ss), size=6, color = 'white')
dev.off()


## -----------------------------------------------------------------------
## SI FIGURE: SPn Crt over time versus viral coinfections
d[site2=='retrospective', site3 := 'Hospital Residual']
d[site2=='community',     site3 := 'Community']

tmp <- d[!is.na(site3) & organism=='SPn',
         .(nv = mean(coinfNV), ct=mean(ct), N = .N), by = .(ym,site3)]

dOrig[, preCOVID := ymd(paste0(ym,'-01')) < ymd('2020-03-01')]

dOrig[site2=='retrospective', site3 := 'Hospital Residual']
dOrig[site2=='community', site3 := 'Community']
tmp2 <- dOrig[!is.na(site3) & organism=='SPn' & N>=16,
              .(pp = mean(present)), by = .(ym,site3,preCOVID)]

png('./figs/si fig spn and viral coinf.png',height=800,width=1200)
ggplot(tmp[N>10]) +
  geom_point(aes(nv,ct,color=ym), size = 4, alpha = .5) + 
  facet_wrap(.~site3) +
  theme_classic() +
  geom_smooth(aes(nv,ct), se = F, color = 'black', method = lm) +
  labs(color='') +
  #scale_color_manual(values=clz) +
  ylab('Average SPn Crt during the month') +
  scale_x_continuous(name = 'Average Viral Coinfections During the Month') +
  theme(text = element_text(size=32 , family = "Garamond"),
        legend.key.size = unit(1, 'cm'),
        axis.text.x  = element_text(size=20),
        legend.text = element_text(size=16),
        legend.position = 'bottom')
dev.off()



## -----------------------------------------------------------------------------
## SI TABLE SHOWING RECRUITMENT SITES

sitab <- d[,.(N=.N), by = .(site,site2)][order(site2,-N)]
write.csv(sitab, 'SI Table Sites.csv')

