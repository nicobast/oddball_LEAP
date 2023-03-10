###

# INFO BY LUKE MASON ON DATA STRUCTURE
# These are the event codes used for marking up the EEG in the MMN task: 
#   
# STD	201
# DEV_FREQ	202
# DEV_DUR	203
# DEV_COMB	204
# 
# In the eye tracking data these become {'EEG_EVENT', XXX}, e.g. {'EEG_EVENT', 201} for the standards. 
# 
# You can get these from the session_events.csv or from the eventBuffer.mat files -- the data is the same, so you can pick whatever format works best for you. 
# 
# Each event has a timestamp. The "RemoteTime" column in the session_events.csv matches with the column of the same name in session_gaze_data.csv. Likewise in the mat files, you'll find the timestamp in eventBuffer.mat links to that in timeBuffer.mat. 

## SETUP ####

#required packages
require(data.table) #fread uses parallelization and thus mucd faster than read.csv
require(readxl) # read demographics data in Excel format
require(pbapply) # progress bar for apply functions
require(zoo) #na.approx


#detect system
ifelse(Sys.info()['sysname']=='Linux',
       home_path<-'~',
       home_path<-'C:/Users/Nico')

##search for the data
data_path<-'Z:/nico/backup_data/data_AIMS/mismatch_negativity/ET_data'
# individual_data_folders<-list.dirs(data_path,recursive=T) #list data folders
# individual_data_folders<-individual_data_folders[which(nchar(individual_data_folders)>(nchar(data_path)+6))] # remove topdirectories without data

#identify relevant data files
individual_data_files<-list.files(data_path,recursive=T,full.names=T)
individual_event_files<-individual_data_files[grepl('session_events.csv',individual_data_files)]
individual_gaze_data_files<-individual_data_files[grepl('session_gaze_data.csv',individual_data_files)]

#only participants with gaze and eventdata
participants_with_gaze_and_events<-substr(individual_gaze_data_files,nchar(data_path)+2,nchar(data_path)+19) %in% substr(individual_event_files,nchar(data_path)+2,nchar(data_path)+19)
individual_gaze_data_files<-individual_gaze_data_files[participants_with_gaze_and_events]


#user-defined functions
fun_rename<-function(x,variable_position,new_name){
  names(x)[variable_position]<-new_name
  return(x)}


### READ DATA + MATCH event to gaze data ####    

#empty list to fill with data
df_list<-list()
#random_participants<-sample(length(individual_data_folders),10)
#random_participants

for(participant in 1:length(individual_event_files)){

#read_data
events<-fread(individual_event_files[participant],sep=',',header=T,fill=T,integer64='numeric') #large data is converted to numeric
et_data<-fread(individual_gaze_data_files[participant],sep=',',header=T,fill=T)
print(paste('read file nr:',participant))

#remove implausible timestamps
events<-events[events$RemoteTime>=0,] 

#reduce events to relevant task events
events<-events[grepl('SCREENFLASH',events$Event) |
                 grepl('BREAK',events$Event) |
                   (grepl('EEG_EVENT',events$Event) & events$EventData %in% c('201','202','203','204'))]

#create empty array of event data to fill with loop
na_array_length_data<-as.numeric(rep(NA,nrow(et_data)))
eventlog<-data.frame(Event=na_array_length_data,
                     EventData=na_array_length_data)

#matching loop
for(row_number in 1:nrow(events)){
which_event<-which(et_data[,'RemoteTime']>=as.numeric(events[row_number,'RemoteTime']) & et_data[,'RemoteTime']<as.numeric(events[row_number+1,'RemoteTime']))
eventlog[which_event,'Event']<-events[row_number,'Event']
eventlog[which_event,'EventData']<-events[row_number,'EventData']

  # #check for et data after end of relevant event data - label as END
  # if(row_number==nrow(events)){
  # which_after_event<-which(et_data[,'RemoteTime']>=as.numeric(max(events$RemoteTime)+median(diff(events$RemoteTime)))) #after 600ms
  # eventlog[which_after_event,'Event']<-'END'
  # eventlog[which_event,'EventData']<-'000'
  # }
  
#print(paste('read:',row_number))
}

#concatenate et data and events
et_data<-data.frame(et_data,eventlog)

#retrieve id
length_pic<-12
id<-substr(individual_data_folders[participant],nchar(individual_data_folders[participant])-length_pic+1,nchar(individual_data_folders[participant]))
wave<-substr(individual_data_folders[participant],nchar(individual_data_folders[participant])-length_pic-5,nchar(individual_data_folders[participant])-length_pic-1)
individual_id<-paste(id,wave,sep='_')

#write to list
df_list[[participant]]<-et_data
names(df_list)[participant]<-individual_id
print(paste('processed:',individual_id))

}

###save temp to file####
save(df_list, file='Z:/nico/backup_data/data_AIMS/mismatch_negativity/mmn_leap_merged') #save temp on NAS

# ##-->plausible
# et_data<-df_list[[sample(1:length(df_list),1)]]
# require(ggplot2)
# ggplot(et_data[et_data$EventData %in% c('201','202','203','204') & et_data$Left.Validity==1,],aes(x=EventData,y=Left.Diameter,group=EventData))+geom_boxplot()
# ###--> oddballs (202-204) usually larger than standard (201)

## ADD  data: demographics + sensorysubgroups + data quality ####
demfile<-paste(home_path,"/PowerFolders/data_LEAP/corelclinical_final050919/LEAP_t1_Core clinical variables_03-09-19-withlabels.xlsx",sep='')
df_dem<-read_excel(demfile, 1, col_names = T, na = c('999','777'))

table(unique(substr(names(df_list),1,12)) %in% unique(df_dem$subjects))
#--> n = 253: Participants with MMN eye-tracking and demographic data

table(df_dem$subjects %in% unique(substr(names(df_list),1,12)),df_dem$t1_group)
# NTC: n = 101; ASD: n = 152 - participants with MMN eye-tracking by group

### process demographics file
# selected_vars<-c('subjects','t1_group','t1_diagnosis','t1_asd_thresh','t1_site',
#                  't1_schedule_adj','t1_sex','t1_ageyrs',
#                  't1_viq','t1_piq','t1_fsiq','t1_ssp_total','t1_rbs_total',
#                  "t1_srs_rawscore_combined","t1_css_total_all","t1_sa_css_all","t1_rrb_css_all",
#                  "t1_adi_social_total","t1_adi_communication_total","t1_adi_rrb_total")
# 
# df_dem_select<-df_dem[,names(df_dem) %in% selected_vars]
# 
# ####--> mental health comorbidities ##
# adhd_inatt<-with(df_dem,ifelse(!is.na(t1_adhd_inattentiv_parent),t1_adhd_inattentiv_parent,t1_adhd_inattentiv_self)) #get ADHD rating from parent and self ratings
# adhd_hyper<-with(df_dem,ifelse(!is.na(t1_adhd_hyperimpul_parent),t1_adhd_hyperimpul_parent,t1_adhd_hyperimpul_self)) #get ADHD rating from parent and self ratings
# 
# 
# anx_beck<-with(df_dem,ifelse(!is.na(t1_beck_anx_adulta_self),t1_beck_anx_adulta_self,
#                              ifelse(!is.na(t1_beck_anx_youthb_self),t1_beck_anx_youthb_self,t1_beck_anx_youthcd_parent
#                              ))) #get ADHD rating from parent and self ratings
# 
# dep_beck<-with(df_dem,ifelse(!is.na(t1_beck_dep_adulta_self),t1_beck_dep_adulta_self,
#                              ifelse(!is.na(t1_beck_dep_youthb),t1_beck_dep_youthb,
#                                     ifelse(!is.na(t1_beck_dep_youthcd),t1_beck_dep_youthcd,t1_beck_dep_adultd_parent)
#                              ))) #get ADHD rating from parent and self ratings
# 
# 
# #MICE imputation of mental health covariates based on sex, age, iq, group, and other covariates
# data_imp<-mice(data.frame(df_dem_select,adhd_inatt,adhd_hyper,anx_beck,dep_beck)[,c('adhd_inatt','adhd_hyper','anx_beck','dep_beck','t1_ageyrs','t1_fsiq','t1_sex','t1_diagnosis')],m=5,maxit=50,meth='pmm',seed=500, printFlag = F)
# df_imputed<-complete(data_imp,5)[,c('adhd_inatt','adhd_hyper','anx_beck','dep_beck')]
# df_dem_select<-data.frame(df_dem_select,df_imputed)
# 
# ###MATCH sensory subgroups 
# df_ssp<-read_xlsx(paste(home_path,"/PowerFolders/data_LEAP/LEAP_t1_sensorysubgroupsTILLMANN.xlsx",sep=''))
# df_sac<-merge(df_sac,df_ssp,by.x='id',by.y='subjects')
# df_fix<-merge(df_fix,df_ssp,by.x='id',by.y='subjects')
# 
# ###MATCH data quality
# df_quality<-read_xlsx(paste(home_path,'/PowerFolders/Paper_AIMS-LEAP_ETcore/LEAP 672+60 Cluster and quality scores.xlsx',sep=''))
# df_quality<-df_quality[,c('ParticipantID','Cluster','SR','Accuracy','Precision','Flicker')]
# df_sac<-merge(df_sac,df_quality,by.x='id',by.y='ParticipantID')
# df_fix<-merge(df_fix,df_quality,by.x='id',by.y='ParticipantID')

## CREATE TIMESTAMP VARIABLE (in s) ####

#create long format variable: ts (in seconds format)
ts<-lapply(df_list,function(x){x<-x$RemoteTime})
ts<-lapply(ts,function(x){abs(head(x,n=1)-x)/1000000})

#add ts and change name of the variable
df_list<-mapply(cbind,df_list,ts,SIMPLIFY = FALSE) #simplify needs to be in capital letters
df_list<-lapply(df_list,fun_rename,variable_position=31,new_name='ts')

#### CREATE TRIAL_NUMBER and ts_trial variable ####

fun_define_trials<-function(single_participant){
  
  #create index variable - indicates when event changes #
  event.change<-which(diff(as.numeric(as.factor(interaction(single_participant$Event,single_participant$EventData))))!=0) #index event change per participant
  length.dataset<-nrow(single_participant) #length of each data set (per participant)
  event.change<-c(0,event.change,length.dataset)
  event.change<-diff(event.change) #create a difference value
  
  #create index for these new trials
  index_trial<-rep(seq_along(event.change),times=event.change)
  
  #sequence over each new trial
  ts_event<-as.numeric(do.call(c,by(index_trial,index_trial,seq_along))) #tts over each trial
  
  #add index and sequence data to list
  single_participant<-data.frame(single_participant,index_trial,ts_event)
  return(single_participant)
}

df_list<-pblapply(df_list,fun_define_trials)

rm(index_trial,ts,ts_event)

### GAZE AND PUPIL PREPROCESSING (TO DO)####

#-- estimate gaze position + calculate center deviation
fun_gaze_position<-function(single_participant){
  
  attach(single_participant)
  #exclude implausible values
  xl <- ifelse((Left.3D.REL.x<0|Left.3D.REL.x>1), NA, Left.3D.REL.x)
  xr <- ifelse((Right.3D.REL.x<0|Right.3D.REL.x>1), NA, Right.3D.REL.x)
  yl <- ifelse((Left.3D.REL.y<0|Left.3D.REL.y>1), NA, Left.3D.REL.y)
  yr <- ifelse((Right.3D.REL.y<0|Right.3D.REL.y>1), NA, Right.3D.REL.y)
  #take offset between left and right into account
  x.offset<-xl-xr
  x.offset<-na.approx(x.offset,rule=2)
  y.offset<-yl-yr
  y.offset<-na.approx(y.offset,rule=2)
  #mean gaze across both eyes
  xl <- ifelse(is.na(xl)==FALSE, xl, xr+x.offset)
  xr <- ifelse(is.na(xr)==FALSE, xr, xl-x.offset)
  yl <- ifelse(is.na(yl)==FALSE, yl, yr+y.offset)
  yr <- ifelse(is.na(yr)==FALSE, yr, yl-y.offset)
  gazepos.x<-(xl+xr)/2
  gazepos.y<-(yl+yr)/2
  
  #remove outside screen
  gazepos.x<-ifelse(gazepos.x>1 | gazepos.x<0,NA,gazepos.x)
  gazepos.y<-ifelse(gazepos.y>1 | gazepos.y<0,NA,gazepos.y)
  
  #estimate center deviation
  center_deviation<-sqrt((gazepos.x-0.5)^2 + (gazepos.y-0.5)^2)
  
  single_participant[,'gazepos.x']<-gazepos.x
  single_participant[,'gazepos.y']<-gazepos.y
  single_participant[,'center_dev']<-center_deviation
  
  detach(single_participant)
  return(single_participant)
}

df_list<-pblapply(df_list,fun_gaze_position)


## -- retrieve screen distance ####
fun_screen_dist<-function(single_participant){
  
  dist_L<-single_participant$Left.3D.UCS.z #in mm from tracker
  dist_R<-single_participant$Right.3D.UCS.z #in mm from tracker
  
  #exclude implausible values (smaller 500mm and larger 800mm is outside track box)
  dist_L<-ifelse(dist_L > 800 | dist_L < 500, NA, dist_L)
  dist_R<-ifelse(dist_R > 800 | dist_R < 500, NA, dist_R)
  
  #take offset between left and right into account
  offset<-dist_L-dist_R
  offset<-na.approx(offset,rule=2)
  
  #mean gaze across both eyes
  dist_L <- ifelse(is.na(dist_L)==FALSE, dist_L, dist_R+offset)
  dist_R <- ifelse(is.na(dist_R)==FALSE, dist_R, dist_R-offset)
  
  screen_dist<-(dist_L+dist_R)/2
  single_participant[,'screen_dist']<-screen_dist
  
  return(single_participant)
  
}

df_list<-pblapply(df_list,fun_screen_dist)

## -- drop unecessary data + particiapnts without data --> subsequent preprocessing requires far less RAM####
fun_required_necessary_data<-function(x){
  
  #drop raw eye tracking data
  x<-x[,!(grepl('.3D.',names(x)))]
  return(x)
  
}
df_list<-pblapply(df_list,fun_required_necessary_data)


#mean of 500k data points per participant
hist(sapply(df_list,nrow))
#most participants have ~450 trials
hist(sapply(df_list,function(x){length(unique(x$index_trial))}))
#number of measurements
paste('date:',paste(Sys.Date()),',individual measurements: n=',length(df_list))







### - calculate baseline pupil size - corrected pupil size ####
fun_baseline<-function(x){
  
  split_by_trial<-split(x,as.factor(x$index_trial))
  baseline_data<-sapply(split_by_trial,function(x){mean(x$pd[x$ts_event<=150],na.rm=T)}) #select between event
  #rpd<-unlist(mapply(function(x,y){x$pd-y},x=split_by_trial,y=baseline_data)) #correct for baseline by subtraction
  rpd<-unsplit(mapply(function(x,y){x$pd-y},x=split_by_trial,y=baseline_data),f=as.factor(x$index_trial)) #correct for baseline by subtraction
  baseline_pd<-rep(baseline_data,sapply(split_by_trial,nrow)) #calculate baseline pd
  
  x[,'rpd']<-rpd
  x[,'baseline_pd']<-baseline_pd
  return(x)
  
}

list_jointatt<-lapply(list_jointatt,fun_baseline)



