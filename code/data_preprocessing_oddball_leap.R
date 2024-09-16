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
require(data.table) #rbindlist - fast row binding of lists
require(mice) #imputation of clincial data
require(RColorBrewer) #color palettes

require(MatchIt) #matching

require(lme4) #mixed models
require(lmerTest)
require(emmeans)

require(ggplot2) #visualization
require(gridExtra)

require(DescTools) #AUC function
require(performance) #r²_nakagawa - marginalized and conditionalk R²

# require(parallel)#
# require(doSNOW)
# #setup parallel processing - specific to windows
# numCores <- detectCores() #number of cores for parallization
# cl<-makeCluster(numCores-2,type="SOCK") #set cluster -under windows
# registerDoSNOW(cl) #
# #stopCluster(cl)

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
start<-Sys.time()

for(participant in 1:length(individual_event_files)){

  #testing
  #participant<-sample(1:length(individual_event_files),1)

#read_data
events<-fread(individual_event_files[participant],sep=',',header=T,fill=T,integer64='numeric') #large data is converted to numeric
et_data<-fread(individual_gaze_data_files[participant],sep=',',header=T,fill=T)
print(paste('read file nr:',participant))

#remove implausible timestamps
events<-events[events$RemoteTime>=0,]

#reduce events to relevant task events - only retain oddball data
  # events<-events[grepl('SCREENFLASH',events$Event) |
  #                  grepl('BREAK',events$Event) |
  #                    (grepl('EEG_EVENT',events$Event) & events$EventData %in% c('201','202','203','204'))]
events<-events[grepl('EEG_EVENT',events$Event) & events$EventData %in% c('201','202','203','204')]

#create empty array of event data to fill with loop
na_array_length_data<-as.numeric(rep(NA,nrow(et_data)))
eventlog<-data.frame(Event=na_array_length_data,
                     EventData=na_array_length_data)

#create an event counter
events$EventCounter<-seq_along(events$Event)

  #matching loop
  for(row_number in 1:nrow(events)){
  which_event<-which(et_data[,'RemoteTime']>=as.numeric(events[row_number,'RemoteTime']) & et_data[,'RemoteTime']<as.numeric(events[row_number+1,'RemoteTime']))
  eventlog[which_event,'Event']<-events[row_number,'Event']
  eventlog[which_event,'EventData']<-events[row_number,'EventData']
  eventlog[which_event,'EventCounter']<-events[row_number,'EventCounter']

  #print(paste('read:',row_number))
  }

#concatenate et data and events
et_data<-data.frame(et_data,eventlog)

#retrieve id
length_pic<-12
id<-substr(individual_event_files[participant],nchar(individual_event_files[participant])-length_pic-nchar("session_events.csv"),nchar(individual_event_files[participant])-nchar("session_events.csv")-nchar('/'))
wave<-substr(individual_event_files[participant],nchar(individual_event_files[participant])-length_pic-nchar("session_events.csv")-nchar('waveX')-nchar('/'),nchar(individual_event_files[participant])-length_pic-nchar("session_events.csv")-nchar('/')*2)
individual_id<-paste(id,wave,sep='_')

#write to list
df_list[[participant]]<-et_data
names(df_list)[participant]<-individual_id
print(paste('processed:',individual_id))

}
Sys.time()-start ###--> 2.5 hours


###save temp to file####
save(df_list, file='C:/Users/nico/Desktop/mmn_leap_merged_24012023') #save temp on NAS


    # ##-->plausible
    et_data<-df_list[[sample(1:length(df_list),1)]]
    require(ggplot2)
    ggplot(et_data[et_data$EventData %in% c('201','202','203','204') & et_data$Left.Validity==1,],aes(x=EventData,y=Left.Diameter,group=EventData))+geom_boxplot()
    ###--> oddballs (202-204) usually larger than standard (201)

## CREATE TIMESTAMP VARIABLE (in s) ####

#create long format variable: ts (in seconds format)
ts<-lapply(df_list,function(x){x<-x$RemoteTime})
ts<-lapply(ts,function(x){abs(head(x,n=1)-x)/1000000})

#add ts and change name of the variable
df_list<-mapply(cbind,df_list,ts,SIMPLIFY = FALSE) #simplify needs to be in capital letters
df_list<-lapply(df_list,fun_rename,variable_position=ncol(df_list[[1]]),new_name='ts') #rename LAST variabel to ts

## ADD ID variable ####
id_list<-names(df_list)
id<-pbmapply(function(x,y){rep(x,nrow(y))},id_list,df_list,SIMPLIFY=FALSE)
df_list<-pbmapply(cbind,df_list,id,SIMPLIFY = FALSE) #simplify needs to be in capital letters
df_list<-lapply(df_list,fun_rename,variable_position=ncol(df_list[[1]]),new_name='id') #rename LAST variabel to ts

#### CREATE TRIAL_NUMBER and ts_trial variable ####

fun_define_trials<-function(single_participant){

  # #create index variable - indicates when event changes #
  # event.change<-which(diff(as.numeric(as.factor(interaction(single_participant$Event,single_participant$EventData))))!=0) #index event change per participant
  # length.dataset<-nrow(single_participant) #length of each data set (per participant)
  # event.change<-c(0,event.change,length.dataset)
  # event.change<-diff(event.change) #create a difference value
  #
  # #create index for these new trials
  # index_trial<-rep(seq_along(event.change),times=event.change)
  ##--> event change not necessary as EventCounter already includes this information

  # #create index for these new trials
  index_trial<-single_participant[,'EventCounter']

  #sequence over each new trial
  split_list<-split(single_participant,f=list(index_trial))
  ts_event<-lapply(split_list,function(x){seq_along(x$EventCounter)})

  #unlist
  single_participant<-rbindlist(split_list)
  ts_event<-unlist(ts_event)

  #add index and sequence data to list
  single_participant<-data.frame(single_participant,ts_event)


  # #sequence over each new trial
  # ts_event<-as.numeric(do.call(c,by(index_trial,index_trial,seq_along))) #tts over each trial
  #
  # #add index and sequence data to list
  # single_participant<-data.frame(single_participant,index_trial,ts_event)
  #
  # #drop previous NA data
  # single_participant<-single_participant[single_participant$index_trial!=0,]


  return(single_participant)
}

df_list<-pblapply(df_list,fun_define_trials)

rm(index_trial,ts,ts_event)

#### remove participatns without MMN data ####
which(sapply(df_list,nrow)==0)
df_list<-df_list[sapply(df_list,nrow)!=0]

### calculate sampling rate and define time_event based on it ####

frequency_rate<-pbsapply(df_list,function(x){median(diff(x$ts),na.rm=T)})
frequency_rate<-ifelse(round(frequency_rate,3)==0.003,0.003333,
                      ifelse(round(frequency_rate,3)==0.008,0.008333,NA))

fun_define_time_event<-function(x,y){

  frequency_rate_i<-rep(y,nrow(x))
  x[,'frequency_rate']<-frequency_rate_i
  x[,'time_event']<-(x$ts_event-1)*x$frequency_rate
  return(x)

}

df_list<-pbmapply(fun_define_time_event,x=df_list,y=frequency_rate, SIMPLIFY = F)
rm(frequency_rate,id)

### REMOVE event=NA data ####

## --> not required - is already dropped during create of TRIAL_NUMBER
#df_list<-pblapply(df_list,function(x){x<-x[!is.na(x$Event),]})

### GAZE AND PUPIL PREPROCESSING ####

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

## -- drop unecessary data --> subsequent preprocessing requires far less RAM####
fun_required_necessary_data<-function(x){

  #drop raw eye tracking data
  x<-x[,!(grepl('.3D.',names(x)))]
  return(x)

}

df_list<-pblapply(df_list,fun_required_necessary_data)


#mean of 30k data points per participant
hist(sapply(df_list,nrow))
#most participants have ~1399 trials
summary(sapply(df_list,function(x){length(unique(x$EventCounter))}))
#number of measurements - 343 (24.01.2023)
paste('date:',paste(Sys.Date()),',individual measurements: n=',length(df_list))

## ----------> PUPIL Preprocessing (Kret, 2018) -------------- ####

  #required function
  fun_blink_cor <- function(signal,lower_threshold,upper_threshold,samples_before,samples_after) {
    #change NA to 999 for rle()-function
    findna <- ifelse(is.na(signal),999,signal)
    #find blinks:
    #output of rle(): how many times values (NA) are repeated
    repets <- rle(findna)
    #stretch to length of PD vector for indexing
    repets <- rep(repets[["lengths"]], times=repets[["lengths"]])
    #difference between two timestamps~3.33ms -> 75/3.333=22.5 -> wenn 23 Reihen PD=NA, dann blink gap
    #if more than 150ms (45 rows) of NA, missing data due to blink unlikely
    #dummy coding of variables (1=at least 23 consecutive repetitions, 0=less than 23 repetitions)
    repets <- ifelse(repets>=lower_threshold & repets<=upper_threshold, 1, 0)
    #exclude cases where other values than NA (999) are repeated >=23 times by changing dummy value to 0:
    repets[findna!=999 & repets==1] <- 0
    #gives out where changes from 0 to 1 (no NA at least 23xNA) or 1 to 0 (23x NA to no NA) appear
    changes <- c(diff(repets),0)
    #define start (interval before blink/missing data)
    changes.start<-which(changes==1) #where NA-sequence starts
    #gives out row numbers of NA (blink) and previous 8 frames
    start.seq<-unlist(lapply(changes.start, function(x) {seq(max(x-(samples_before-1),1), x)}))
    repets[start.seq]<-1
    #define end (interval after blink/missing data)
    changes.end<-which(changes==-1)+1 #where NA.sequence ends
    #gives out row numbers of NA (blink) and subsequent 8 frames
    end.seq<-unlist(lapply(changes.end, function(x) {seq(x, min(x+(samples_before-1),length(repets)))}))
    repets[end.seq]<-1
    #replace PD data in blink interval (start to end) with NA
    signal[repets==1]<-NA
    return(signal)
  }

func_pd_preprocess<-function(x){

  #define variables
  Left_Diameter<-x$Left.Diameter
  Right_Diameter<-x$Right.Diameter
  RemoteTime<-x$RemoteTime

  #constant for MAD caluclation
  constant<-3 ##--> if change speed is higher than constant * median change --> values are excluded
  #constant<-3 #default value

  # STEP 1 - exclude invalid data ####
  pl <- ifelse((Left_Diameter<2|Left_Diameter>8), NA, Left_Diameter)
  pr <- ifelse((Right_Diameter<2|Right_Diameter>8), NA, Right_Diameter)
  #table(is.na(pl))
  #table(is.na(pr))

  # STEP 2 - filtering ####
  ## A) normalized dilation speed, take into account time jumps with Remotetimestamps: ####
  #maximum change in pd compared to last and next pd measurement
  #Left
  pl.speed1<-diff(pl)/diff(RemoteTime) #compared to last
  pl.speed2<-diff(rev(pl))/diff(rev(RemoteTime)) #compared to next
  pl.speed1<-c(NA,pl.speed1)
  pl.speed2<-c(rev(pl.speed2),NA)
  pl.speed<-pmax(pl.speed1,pl.speed2,na.rm=T)
  rm(pl.speed1,pl.speed2)
  #Right
  pr.speed1<-diff(pr)/diff(RemoteTime)
  pr.speed2<-diff(rev(pr))/diff(rev(RemoteTime))
  pr.speed1<-c(NA,pr.speed1)
  pr.speed2<-c(rev(pr.speed2),NA)
  pr.speed<-pmax(pr.speed1,pr.speed2,na.rm=T)
  rm(pr.speed1,pr.speed2)
  #median absolute deviation -SPEED
  #constant<-3
  pl.speed.med<-median(pl.speed,na.rm=T)
  pl.mad<-median(abs(pl.speed-pl.speed.med),na.rm = T)
  pl.treshold.speed<-pl.speed.med+constant*pl.mad #treshold.speed units are mm/microsecond
  #plot(abs(pl.speed))+abline(h=pl.treshold.speed)
  pr.speed.med<-median(pr.speed,na.rm=T)
  pr.mad<-median(abs(pr.speed-pr.speed.med),na.rm = T)
  pr.treshold.speed<-pr.speed.med+constant*pr.mad #treshold.speed units are mm/microsecond
  #plot(abs(pr.speed))+abline(h=pr.treshold.speed)
  #correct pupil dilation for speed outliers
  pl<-ifelse(abs(pl.speed)>pl.treshold.speed,NA,pl)
  pr<-ifelse(abs(pr.speed)>pr.treshold.speed,NA,pr)

  ## B) delete data around blinks ####
  #gaps=missing data sections > 75ms; Leonie: also <=250ms, otherwise not likely to be a blink
  #to be excluded: samples within 50 ms of gaps -> +-25 (8 data points) oder 50?
    #take sampling rate into account (300 vs. 120) - define individual values:
    deletion_span<-25
    shortest_blink<-75 #measured in ms
    longest_blink<-250 #measured in ms
    individial_lower_threshold<-round(shortest_blink/median(diff(RemoteTime),na.rm=T))
    individial_upper_threshold<-round(longest_blink/median(diff(RemoteTime),na.rm=T))
    individial_deletion_span<-round(deletion_span/median(diff(RemoteTime),na.rm=T))
   pl<-fun_blink_cor(pl,lower_threshold=individial_lower_threshold,upper_threshold=individial_upper_threshold,samples_before=individial_deletion_span,samples_after=individial_deletion_span)
   pr<-fun_blink_cor(pr,lower_threshold=individial_lower_threshold,upper_threshold=individial_upper_threshold,samples_before=individial_deletion_span,samples_after=individial_deletion_span)

  ## C) normalized dilation size - median absolute deviation -SIZE ####
  #applies a two pass approach
  #first pass: exclude deviation from trend line derived from all samples
  #second pass: exclude deviation from trend line derived from samples passing first pass
  #-_> reintroduction of sample that might have been falsely excluded due to outliers
  #estimate smooth size based on sampling rate
  smooth.length<-150 #measured in ms
  #take sampling rate into account (300 vs. 120):
  smooth.size<-round(smooth.length/median(diff(RemoteTime),na.rm=T)) #timestamp resolution in milliseconds
  is.even<-function(x){x%%2==0}
  smooth.size<-ifelse(is.even(smooth.size)==T,smooth.size+1,smooth.size) #make sure to be odd value (see runmed)
  #Left
  pl.smooth<-na.approx(pl,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - FIRST pass
  pl.pass1<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
  pr.pass1<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
  #Left
  pl.smooth<-na.approx(pl.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - SECOND pass
  pl.pass2<-ifelse((pl>pl.smooth+constant*pl.mad)|(pl<pl.smooth-constant*pl.mad),NA,pl)
  pr.pass2<-ifelse((pr>pr.smooth+constant*pr.mad)|(pr<pr.smooth-constant*pr.mad),NA,pr)
  pl<-pl.pass2
  pr<-pr.pass2

  ## D) sparsity filter - not applied ####
  # STEP 3 - processing valid samples  ####
  #take offset between left and right into account
  pd.offset<-pl-pr
  pd.offset<-na.approx(pd.offset,rule=2)
  #mean pupil dilation across both eyes
  pl <- ifelse(is.na(pl)==FALSE, pl, pr+pd.offset)
  pr <- ifelse(is.na(pr)==FALSE, pr, pl-pd.offset)

  #interpolation of NA (for <=300ms)
  interpolation_span<-100 #in ms
  individual_interpolation_samples<-round(interpolation_span/median(diff(RemoteTime),na.rm=T))
  pl<-na.approx(pl, na.rm=F, maxgap=individual_interpolation_samples, rule=2)
  pr<-na.approx(pr, na.rm=F, maxgap=individual_interpolation_samples, rule=2)

  pd <- (pl+pr)/2

  #calculate baseline pd - mean of first 250ms
  baseline_pd<-mean(pd[x$time_event<=0.250],na.rm=T)
  baseline_pd<-rep(baseline_pd,nrow(x))

  #return pupillometry values
  x[,'pd_baseline']<-baseline_pd
  x[,'pd']<-pd
  x[,'rpd']<-pd-baseline_pd
  return(x)
}

#split each participant to trials
df_list_split<-pblapply(df_list,function(x){split(x,f=as.factor(x$EventCounter))})

df_list_split<-pblapply(df_list_split,function(x){lapply(x,func_pd_preprocess)})
#--> rerun with revised, currently takes 60 minutes

df_list<-pblapply(df_list_split,rbindlist)
#rm(df_list_split)

##--> LOAD preprocessed data ####

#save(df_list,file='C:/Users/nico/Desktop/mmn_leap_pd_preproccesed_24012023')
load(file='C:/Users/nico/Desktop/mmn_leap_pd_preproccesed_24012023')


  # # ##-->plausible responses per individual
  # et_data<-df_list[[sample(1:length(df_list),1)]]
  # ggplot(et_data,aes(x=EventData,y=pd,group=EventData))+geom_boxplot()
  # ###--> oddballs (202-204) usually larger than standard (201)
  #   ggplot(et_data,aes(x=time_event,y=pd,group=as.factor(EventData),color=as.factor(EventData)))+geom_smooth()
  #   ggplot(et_data,aes(x=EventCounter,y=pd))+geom_smooth()

## bind to data frame
df<-rbindlist(df_list) #very fast
length(unique(df$id))

    # #testing
    # selected_sample<-sample(1:length(df_list),120)
    # test<-df[sample(1:nrow(df),nrow(df)/10),]
    #
    # ggplot(test[test$time_event<0.6,],aes(x=time_event,y=rpd,group=EventData,color=EventData))+geom_smooth()+theme_bw()
    # ggplot(test[test$EventCounter<1400,],aes(x=EventCounter,y=pd_baseline))+geom_smooth()
    #
    # ggplot(test,aes(x=ts,y=pd))+geom_smooth(method='lm')


## CREATE TRIAL specific data frame ####

    #split each participant to trials ~ 2min
    df_list_split<-pblapply(df_list,function(x){split(x,f=as.factor(x$EventCounter))})

    #function to retrieve pupillary response
    fun_trial_rpd_response<-function(x){

      rpd_end<-mean(x$rpd[x$time_event>=0.5],na.rm=T)
      rpd_start<-mean(x$rpd[x$time_event<0.2],na.rm=T)
      #rpd_start<-mean(x$rpd[x$time_event<=0.2 & x$time_event>=0.1],na.rm=T)

      rpd_response<-rpd_end-rpd_start
      return(rpd_response)

    }

    #apply function to list (TRIAL split list)
    list_rpd_response<-pblapply(df_list_split,function(x,y){sapply(x,fun_trial_rpd_response)})

    # start<-Sys.time()
    # list_rpd_response2<-pblapply(df_list_split,function(x,y){parSapply(cl,x,fun_trial_rpd_response)})
    # Sys.time()-start

    ##retrieve area under curve of response
    fun_trial_rpd_auc<-function(x){

      # rel_data<-x[x$time_event>=0.2,] #only calculate AUC after 200ms
      # rpd_min_cor<-rel_data$rpd-min(rel_data$rpd,na.rm=T) #correct for negative value that are normally substracted from the AUC
      rpd_auc<-AUC(x$time_event,x$rpd,na.rm=T)
      return(rpd_auc)

    }

    list_rpd_auc<-pblapply(df_list_split,function(x,y){sapply(x,fun_trial_rpd_auc)})

    ##return missing data per trial
    fun_return_missing_trial<-function(x){

      sum(is.na(x$rpd)==T)/length(x$rpd)

    }

    list_missing_data<-pblapply(df_list_split,function(x,y){sapply(x,fun_return_missing_trial)})

    ##return sequence position counter
    fun_sequence_pos<-function(x){

    events<-as.numeric(sapply(x,function(y){head(y$EventData,n=1)}))

    # last_event<-c(0,diff(events))
    # last_event<-ifelse(last_event==-3,204,ifelse(last_event==-2,203,ifelse(last_event==-1,202,201)))
    last_event<-c(201,events)[-length(events)] #assuming first event is always a standard = 201

      sequence_position<-as.numeric(0)
      current_sequence<-0
      for(i in 1:length(events)){

              if(events[i]==201 & last_event[i]==201){current_sequence<-current_sequence+1}
              if(events[i]==201 & last_event[i]!=201){current_sequence<-1}
              if(events[i]!=201 & last_event[i]!=201){current_sequence<-1}
              if(events[i]!=201 & last_event[i]==201){current_sequence<-current_sequence+1}
              # if(events[i]!=201){current_sequence<-current_sequence+1}
              sequence_position[i]<-current_sequence

      }

    return(sequence_position)
    }

    list_sequence_position<-pblapply(df_list_split,fun_sequence_pos)


    #fun_stimulus_occurance
    fun_stimulus_occurance<-function(x){

      events<-as.numeric(sapply(x,function(y){head(y$EventData,n=1)}))
      specific_event_counter<-rep(NA,length(events))
      specific_event_counter[events==201]<-1:length(events[events==201])
      specific_event_counter[events==202]<-1:length(events[events==202])
      specific_event_counter[events==203]<-1:length(events[events==203])
      specific_event_counter[events==204]<-1:length(events[events==204])
      return(specific_event_counter)

    }

    list_specificevents_counter<-pblapply(df_list_split,fun_stimulus_occurance)

    #function to aggregate to trials
    fun_aggregate_trial<-function(x){

      #aggregation functions
      subject_by_trial_factors<-aggregate(x[,c('id','EventData')],by=list(x$EventCounter),FUN=head,n=1)
      subject_by_trial_numeric<-aggregate(x[,c('gazepos.x','gazepos.y','center_dev','screen_dist','pd','pd_baseline','frequency_rate')],by=list(x$EventCounter),FUN=mean,na.rm=T)
      subject_by_trial_numeric<-subject_by_trial_numeric[,!names(subject_by_trial_numeric) %in% c('Group.1')]
      subject_by_trial_samples<-aggregate(x[,c('ts_event')],by=list(x$EventCounter),max,na.rm=T)
      subject_by_trial_duration<-aggregate(x[,c('time_event')],by=list(x$EventCounter),max,na.rm=T)

      subject_by_trial<-data.frame(subject_by_trial_factors,
                                   subject_by_trial_numeric,
                                   number_of_samples=subject_by_trial_samples[,'ts_event'],
                                   trial_duration=subject_by_trial_duration[,'time_event'])
      names(subject_by_trial)[1]<-'EventCounter'
      return(subject_by_trial)

    }

    #apply function to list (PARTICIPANT split list)
    df_trial<-pblapply(df_list,fun_aggregate_trial) ##takes three minutes

    df_trial<-pbmapply(data.frame,df_trial,list_rpd_response,SIMPLIFY=F)
    df_trial<-lapply(df_trial,fun_rename,variable_position=ncol(df_trial[[1]]),new_name='rpd_response') #rename LAST variable

    df_trial<-pbmapply(data.frame,df_trial,list_rpd_auc,SIMPLIFY=F)
    df_trial<-lapply(df_trial,fun_rename,variable_position=ncol(df_trial[[1]]),new_name='rpd_auc') #rename LAST variable

    df_trial<-pbmapply(data.frame,df_trial,list_missing_data,SIMPLIFY=F)
    df_trial<-lapply(df_trial,fun_rename,variable_position=ncol(df_trial[[1]]),new_name='missing_data_trial') #rename LAST variabele

    df_trial<-pbmapply(data.frame,df_trial,list_sequence_position,SIMPLIFY=F)
    df_trial<-lapply(df_trial,fun_rename,variable_position=ncol(df_trial[[1]]),new_name='sequence_position') #rename LAST variable

    df_trial<-pbmapply(data.frame,df_trial,list_specificevents_counter,SIMPLIFY=F)
    df_trial<-lapply(df_trial,fun_rename,variable_position=ncol(df_trial[[1]]),new_name='specific_events_counter') #rename LAST variable

    df_trial<-rbindlist(df_trial)

    #add subjects variable
    df_trial$subjects<-substr(df_trial$id,1,12)

      #remove participants with low or high number of trials
      subjects_with_different_trials<-names(sapply(df_list_split,length)[sapply(df_list_split,length)<1100 | sapply(df_list_split,length)>1700])
      #they have already been dropped, see:
      df_trial<-df_trial[!df_trial$subjects %in% subjects_with_different_trials,]

    # #diagnosistic
    # with(df_trial,by(rpd_response,EventData,summary))
    # #missings in trials
    # with(df_trial,table(is.na(rpd_response),EventData))

    ##testing
      tmp<-df_trial

      df_trial<-tmp

## REMOVE trials with low data quality ####
with(df_trial,hist(missing_data_trial))
with(df_trial,hist(center_dev))
with(df_trial,hist(rpd_response))
with(df_trial,hist(rpd_auc[rpd_auc<0.5 & rpd_auc>-0.5]))

table(df_trial$missing_data_trial<=0.2)
table(df_trial$center_dev<=0.3)
table(df_trial$rpd_response<1 & df_trial$rpd_response>-1)
table(df_trial$rpd_auc<0.5 & df_trial$rpd_auc>-0.5)

df_trial<-df_trial[df_trial$missing_data_trial<=0.2,]
df_trial<-df_trial[df_trial$center_dev<=0.3,]
df_trial<-df_trial[df_trial$rpd_response<1 & df_trial$rpd_response>-1,]
df_trial<-df_trial[df_trial$rpd_auc<0.5 & df_trial$rpd_auc>-0.5,]

with(df_trial,hist(EventCounter[EventCounter<1400]))
df_trial<-df_trial[df_trial$EventCounter<=1400,]

## create wave + sampling rate variable ####
df_trial$wave<-substr(df_trial$id,14,18)
df_trial$sampling_rate<-with(df_trial,ifelse(frequency_rate==0.003333,'300Hz',ifelse(frequency_rate==0.008333,'120Hz','Other')))


## LOAD + ADD DEMOGRAPHICS DATA  ####

#demfile<-paste0(home_path,"/PowerFolders/data_LEAP/corelclinical_final050919/LEAP_t1_Core clinical variables_03-09-19-withlabels.xlsx")
#t1 and t2 data
demfile<-paste0(home_path,"/PowerFolders/data_LEAP/corelclinical_final050919/LEAP_t1_t2_Core clinical variables_03-09-19-withlabels.xlsx")

df_dem<-read_excel(demfile, 1, col_names = T, na = c('999','777'))

ids_et_data<-unique(substr(df_trial$id,1,12))

table(ids_et_data %in% unique(df_dem$subjects))
#--> n = 251: unique participants with MMN eye-tracking and demographic data
# NTC: n = 100; ASD: n = 151 - participants with MMN eye-tracking by group

#change to characters
df_dem$subjects<-as.character(df_dem$subjects)


#select only participants with ET data
df_dem<-df_dem[unique(df_dem$subjects) %in% ids_et_data,]

###--> CREATE DF TIMEPOINT - per timepoint data ####
fun_df_timepoints<-function(x){

  #split into wave 1 and wave 2 data
  subjects<-x[,names(x)=='subjects']
  df_dem_wave1<-x[,grepl('t1',names(x))]
  df_dem_wave2<-x[,grepl('t2',names(x))]

  #separate variables that are unique and non-unique in wide format
  unique_variables_wave1<-!(substr(names(df_dem_wave1),4,nchar(names(df_dem_wave1))) %in% substr(names(df_dem_wave2),4,nchar(names(df_dem_wave2))))
  unique_variables_wave2<-!(substr(names(df_dem_wave2),4,nchar(names(df_dem_wave2))) %in% substr(names(df_dem_wave1),4,nchar(names(df_dem_wave1))))

  df_dem_unique_wave1<-df_dem_wave1[,unique_variables_wave1]
  df_dem_notunique_wave1<-df_dem_wave1[,!unique_variables_wave1]
  df_dem_unique_wave2<-df_dem_wave2[,unique_variables_wave2]
  df_dem_notunique_wave2<-df_dem_wave2[,!unique_variables_wave2]

  #allign names for non unique variables
  names(df_dem_notunique_wave1)<-substr(names(df_dem_notunique_wave1),4,nchar(names(df_dem_notunique_wave1)))
  names(df_dem_notunique_wave2)<-substr(names(df_dem_notunique_wave2),4,nchar(names(df_dem_notunique_wave2)))

  df_dem_wave1<-data.frame(subjects,df_dem_unique_wave1)
  df_dem_wave2<-data.frame(subjects,df_dem_unique_wave2)

  #add wave variable
  df_dem_wave1$wave<-'wave1'
  df_dem_wave2$wave<-'wave2'

  df_dem_notunique<-rbind(df_dem_notunique_wave1,df_dem_notunique_wave2)
  df_timepoint<-rbind(df_dem_wave1[,c('subjects','wave')],df_dem_wave2[,c('subjects','wave')])
  df_timepoint<-cbind(df_timepoint,df_dem_notunique)

  df_dem_unique_wave1<-rbind(df_dem_unique_wave1,df_dem_unique_wave1)
  df_dem_unique_wave2<-rbind(df_dem_unique_wave2,df_dem_unique_wave2)
  df_timepoint<-cbind(df_timepoint,df_dem_unique_wave1,df_dem_unique_wave2)


  #create merging variable
  df_timepoint$id<-paste(df_timepoint$subjects,df_timepoint$wave,sep='_')

  # df_dem_wave1$id<-paste(df_dem_wave1$subjects,df_dem_wave1$wave,sep='_')
  # df_dem_wave2$id<-paste(df_dem_wave2$subjects,df_dem_wave2$wave,sep='_')
  #
  # #exclude variable that are also in merged df
  # df_dem_wave1<-df_dem_wave1[,!(names(df_dem_wave1) %in% c('subjects','wave'))]
  # df_dem_wave2<-df_dem_wave2[,!(names(df_dem_wave2) %in% c('subjects','wave'))]
  #
  # #merge
  # df_timepoint<-merge(df_timepoint,df_dem_wave1,by='id',all.x=T)
  # df_timepoint<-merge(df_timepoint,df_dem_wave2,by='id',all.x=T)

  return(df_timepoint)

}

df_timepoint<-fun_df_timepoints(df_dem)

#drop those that have no et data
df_timepoint<-df_timepoint[df_timepoint$id %in% unique(df_trial$id),]
###--> n = 339

with(df_timepoint,table(wave,t1_diagnosis))
with(df_timepoint,length(unique(subjects)))
with(df_timepoint,by(ageyrs,wave,summary))
with(df_timepoint,table(wave))

names(df_timepoint)

###--> define timepoint variable - identify which participants have which measurement timepoints completed ####
subjects<-rownames(with(df_timepoint,table(subjects,wave)))
timepoints<-ifelse(with(df_timepoint,table(subjects,wave)[,1]==1 & table(subjects,wave)[,2]==0),'only_wave1',
                   ifelse(with(df_timepoint,table(subjects,wave)[,1]==0 & table(subjects,wave)[,2]==1),'only_wave2','wave1+2'))

df_subjects<-data.frame(subjects,timepoints)
df_dem<-merge(df_dem,df_subjects,by='subjects')
table(df_dem$timepoints)
###--> n=75 only wave 1
###--> n=88 only wave 2
###--> n=88 wave 1 & wave 2

###- ADD sensory subgroups (see Tillmann, 2021) ####
df_ssp<-read_excel(paste0(home_path,"/PowerFolders/data_LEAP/LEAP_t1_sensorysubgroupsTILLMANN.xlsx"))
df_dem<-merge(df_dem,df_ssp,by='subjects')

###- ADD data quality ####
df_quality<-read_excel(paste0(home_path,'/PowerFolders/Paper_AIMS-LEAP_ETcore/data/LEAP 672+60 Cluster and quality scores.xlsx'))
df_quality<-df_quality[,c('ParticipantID','Cluster','SR','Accuracy','Precision','Flicker')]

table(df_dem$subjects %in% df_quality$ParticipantID)
df_dem<-merge(df_dem,df_quality,by.x='subjects',by.y='ParticipantID',all.x=T)
###-->would lose 31 subject that have no data quality scores - ergo also no information on sampling rate

## - ADD data quality (from df_trial) ####
df_missingdata<-data.frame(subjects=unique(df_trial$subjects),
                           missing_data_trial=with(df_trial,as.numeric(by(missing_data_trial,subjects,mean,na.rm=T))),
                           sampling_rate=with(df_trial,as.factor(by(sampling_rate,subjects,head,n=1))))
df_dem<-merge(df_dem,df_missingdata,by='subjects')

df_missingdata<-data.frame(id=unique(df_trial$id),
                           missing_data_trial=with(df_trial,as.numeric(by(missing_data_trial,id,mean,na.rm=T))),
                           sampling_rate=with(df_trial,as.factor(by(sampling_rate,id,head,n=1))))
df_timepoint<-merge(df_timepoint,df_missingdata,by='id')

## - ADD eye tracking variables ####
df_et_vars<-aggregate(df_trial[,c('center_dev','screen_dist','pd','pd_baseline','rpd_auc','rpd_response')],by=list(df_trial$subjects),FUN=mean,na.rm=T)
df_et_vars_rpd<-aggregate(df_trial[,c('rpd_response')],by=list(df_trial$subjects,df_trial$EventData),FUN=mean,na.rm=T)
df_et_vars_rpd<-reshape(df_et_vars_rpd, idvar = "Group.1", timevar = "Group.2", direction = "wide") #long to wide format
df_et_vars_rpd_auc<-aggregate(df_trial[,c('rpd_auc')],by=list(df_trial$subjects,df_trial$EventData),FUN=mean,na.rm=T)
df_et_vars_rpd_auc<-reshape(df_et_vars_rpd_auc, idvar = "Group.1", timevar = "Group.2", direction = "wide") #long to wide format
names(df_et_vars)[1]<-'subjects'
names(df_et_vars_rpd)[1]<-'subjects'
names(df_et_vars_rpd_auc)[1]<-'subjects'

df_dem<-merge(df_dem,df_et_vars,by='subjects')
df_dem<-merge(df_dem,df_et_vars_rpd,by='subjects')
df_dem<-merge(df_dem,df_et_vars_rpd_auc,by='subjects')


df_et_vars<-aggregate(df_trial[,c('center_dev','screen_dist','pd','pd_baseline','rpd_auc')],by=list(df_trial$id),FUN=mean,na.rm=T)
df_et_vars_rpd<-aggregate(df_trial[,c('rpd_response')],by=list(df_trial$id,df_trial$EventData),FUN=mean,na.rm=T)
df_et_vars_rpd<-reshape(df_et_vars_rpd, idvar = "Group.1", timevar = "Group.2", direction = "wide") #long to wide format
df_et_vars_rpd_auc<-aggregate(df_trial[,c('rpd_auc')],by=list(df_trial$id,df_trial$EventData),FUN=mean,na.rm=T)
df_et_vars_rpd_auc<-reshape(df_et_vars_rpd_auc, idvar = "Group.1", timevar = "Group.2", direction = "wide") #long to wide format
names(df_et_vars)[1]<-'id'
names(df_et_vars_rpd)[1]<-'id'
names(df_et_vars_rpd_auc)[1]<-'id'

df_timepoint<-merge(df_timepoint,df_et_vars,by='id')
df_timepoint<-merge(df_timepoint,df_et_vars_rpd,by='id')
df_timepoint<-merge(df_timepoint,df_et_vars_rpd_auc,by='id')

### IMPUTE missing clincial data ####
selected_vars<-c('subjects','t1_group','t1_diagnosis','t1_asd_thresh','t1_site',
                 't1_schedule_adj','t1_sex','t1_ageyrs',
                 't1_viq','t1_piq','t1_fsiq','t1_ssp_total','t1_rbs_total',
                 "t1_srs_rawscore_combined","t1_css_total_all","t1_sa_css_all","t1_rrb_css_all",
                 "t1_adi_social_total","t1_adi_communication_total","t1_adi_rrb_total")

df_dem_select<-df_dem[,names(df_dem) %in% selected_vars]

####--> mental health comorbidities ##
adhd_inatt<-with(df_dem,ifelse(!is.na(t1_adhd_inattentiv_parent),t1_adhd_inattentiv_parent,t1_adhd_inattentiv_self)) #get ADHD rating from parent and self ratings
adhd_hyper<-with(df_dem,ifelse(!is.na(t1_adhd_hyperimpul_parent),t1_adhd_hyperimpul_parent,t1_adhd_hyperimpul_self)) #get ADHD rating from parent and self ratings

anx_beck<-with(df_dem,ifelse(!is.na(t1_beck_anx_adulta_self),t1_beck_anx_adulta_self,
                             ifelse(!is.na(t1_beck_anx_youthb_self),t1_beck_anx_youthb_self,t1_beck_anx_youthcd_parent
                             )))

dep_beck<-with(df_dem,ifelse(!is.na(t1_beck_dep_adulta_self),t1_beck_dep_adulta_self,
                             ifelse(!is.na(t1_beck_dep_youthb),t1_beck_dep_youthb,
                                    ifelse(!is.na(t1_beck_dep_youthcd),t1_beck_dep_youthcd,t1_beck_dep_adultd_parent)
                             )))

#missings
table(is.na(adhd_inatt))[2]/length(adhd_inatt)
table(is.na(adhd_hyper))[2]/length(adhd_hyper)
table(is.na(anx_beck))[2]/length(anx_beck)
table(is.na(dep_beck))[2]/length(dep_beck)

#MICE imputation of mental health covariates based on sex, age, iq, group, and other covariates
require(mice)
data_imp<-mice(data.frame(df_dem_select,adhd_inatt,adhd_hyper,anx_beck,dep_beck)[,c('adhd_inatt','adhd_hyper','anx_beck','dep_beck','t1_ageyrs','t1_viq','t1_piq','t1_fsiq','t1_sex','t1_diagnosis')],m=5,maxit=50,meth='pmm',seed=500, printFlag = F)
df_imputed<-complete(data_imp,5)[,c('adhd_inatt','adhd_hyper','anx_beck','dep_beck','t1_piq')]

df_imputed<-data.frame(subjects=df_dem$subjects,df_imputed)
df_dem<-df_dem[,!names(df_dem) %in% c('t1_piq')]
df_dem<-merge(df_dem,df_imputed,by='subjects')
#df_timepoint<-merge(df_timepoint,df_imputed,by='subjects')


### REMOVE participants with missing IQ ####

table(is.na(df_dem$t1_piq))
table(is.na(df_timepoint$t1_piq))

df_dem<-df_dem[!is.na(df_dem$t1_piq),]
df_timepoint<-df_timepoint[!is.na(df_timepoint$t1_piq),]

#### MATCH SAMPLE ####

groupBoo<-with(df_dem,ifelse(t1_diagnosis=='ASD',1,0))

with(df_dem,t.test(t1_piq~t1_diagnosis))
with(df_dem,t.test(missing_data_trial~t1_diagnosis))
with(df_dem,t.test(center_dev~t1_diagnosis))

#ALL - matching
set.seed(100)
all.match<-matchit(groupBoo~t1_piq+center_dev,
                   data=df_dem,
                   method='nearest',discard='both',
                   ratio=4, #match four controls to each ASD
                   replace=T,caliper=0.05)

summary(all.match)[['nn']]

#remove unmatched cases - from eye tracking data
all.match<-match.data(all.match)

with(all.match,t.test(t1_piq~t1_diagnosis))
with(all.match,t.test(missing_data_trial~t1_diagnosis))
with(all.match,t.test(center_dev~t1_diagnosis))

#reduce data sets based on matched sample
df_dem<-df_dem[df_dem$subjects %in% all.match$subjects,]
df_trial<-df_trial[df_trial$subjects %in% all.match$subjects,]

# #also do for df --> happens down below in visualization

#### sample descriptives - PARTICIPANT specific ####

fun_return_descriptive<-function(variable,group,rounding=2,scaling=1){
  mean_values<-by(variable,group,function(x){round(mean(x,na.rm=T),rounding)})
  sd_values<-by(variable,group,function(x){round(sd(x,na.rm=T),rounding)})
  min_values<-by(variable,group,function(x){round(min(x,na.rm=T),rounding)})
  max_values<-by(variable,group,function(x){round(max(x,na.rm=T),rounding)})
  mean_values<-mean_values*scaling
  sd_values<-sd_values*scaling
  min_values<-min_values*scaling
  max_values<-max_values*scaling
  paste0(mean_values,'/',sd_values,' ','[',min_values,'-',max_values,']')
}

#standard descriptives
n_group<-with(df_dem,by(subjects,t1_diagnosis,function(x){length(unique(x))}))
gender<-with(df_dem,by(t1_sex,t1_diagnosis,function(x){paste0(table(x)[2],'/',table(x)[1])}))
timepoints<-with(df_dem,by(timepoints,t1_diagnosis,function(x){paste0(table(x)[1],'/',table(x)[2],'/',table(x)[3])}))
age<-with(df_dem,fun_return_descriptive(variable=t1_ageyrs,group=t1_diagnosis))
iq<-with(df_dem,fun_return_descriptive(variable=t1_fsiq,group=t1_diagnosis))
piq<-with(df_dem,fun_return_descriptive(variable=t1_piq,group=t1_diagnosis))
viq<-with(df_dem,fun_return_descriptive(variable=t1_viq,group=t1_diagnosis))

#clinical variables
srs<-with(df_dem,fun_return_descriptive(variable=t1_srs_rawscore,group=t1_diagnosis))
rbs_total<-with(df_dem,fun_return_descriptive(variable=t1_rbs_total,group=t1_diagnosis))
sdq_total<-with(df_dem,fun_return_descriptive(variable=t1_sdq_total_difficulties_p,group=t1_diagnosis))
adhd_inatt<-with(df_dem,fun_return_descriptive(variable=adhd_inatt,group=t1_diagnosis))
adhd_hyper<-with(df_dem,fun_return_descriptive(variable=adhd_hyper,group=t1_diagnosis))
anx_beck<-with(df_dem,fun_return_descriptive(variable=anx_beck,group=t1_diagnosis))
dep_beck<-with(df_dem,fun_return_descriptive(variable=dep_beck,group=t1_diagnosis))

##data quality
#precision<-with(df_dem,fun_return_descriptive(variable=Precision,group=t1_diagnosis))
#accuracy<-with(df_dem,fun_return_descriptive(variable=Accuracy,group=t1_diagnosis))
#sampling_rate<-with(df_dem,by(sampling_rate,t1_diagnosis,function(x){paste0(table(x)[2],'/',table(x)[1])}))
missing_data_trial<-with(df_dem,fun_return_descriptive(variable=missing_data_trial,group=t1_diagnosis,rounding=4,scaling = 100))
sampling_rate<-with(df_timepoint,by(sampling_rate,t1_diagnosis,function(x){paste0(table(x)[2],'/',table(x)[1])}))

#et vars
center_dev_id<-with(df_dem,fun_return_descriptive(variable=center_dev,group=t1_diagnosis,rounding=4,scaling=100))
screen_dist_id<-with(df_dem,fun_return_descriptive(variable=screen_dist,group=t1_diagnosis))
pd_baseline_id<-with(df_dem,fun_return_descriptive(variable=pd_baseline,group=t1_diagnosis))
rpd_response201_id<-with(df_dem,fun_return_descriptive(variable=rpd_response.201,group=t1_diagnosis,rounding=4,scaling=100))
rpd_response202_id<-with(df_dem,fun_return_descriptive(variable=rpd_response.202,group=t1_diagnosis,rounding=4,scaling=100))
rpd_response203_id<-with(df_dem,fun_return_descriptive(variable=rpd_response.203,group=t1_diagnosis,rounding=4,scaling=100))
rpd_response204_id<-with(df_dem,fun_return_descriptive(variable=rpd_response.204,group=t1_diagnosis,rounding=4,scaling=100))
rpd_auc201_id<-with(df_dem,fun_return_descriptive(variable=rpd_auc.201,group=t1_diagnosis,rounding=4,scaling=100))
rpd_auc202_id<-with(df_dem,fun_return_descriptive(variable=rpd_auc.202,group=t1_diagnosis,rounding=4,scaling=100))
rpd_auc203_id<-with(df_dem,fun_return_descriptive(variable=rpd_auc.203,group=t1_diagnosis,rounding=4,scaling=100))
rpd_auc204_id<-with(df_dem,fun_return_descriptive(variable=rpd_auc.204,group=t1_diagnosis,rounding=4,scaling=100))


rpd_auc_id<-with(df_dem,fun_return_descriptive(variable=rpd_auc,group=t1_diagnosis,rounding=4,scaling=100))


descriptives_table<-rbind(n_group,gender,timepoints,age,iq,piq,viq,
                          srs,rbs_total,sdq_total,adhd_inatt,adhd_hyper,anx_beck,dep_beck,
                          #precision,accuracy,
                          #sampling_rate,
                          missing_data_trial,sampling_rate,
                          center_dev_id,screen_dist_id,pd_baseline_id,
                          rpd_response201_id,rpd_response202_id,rpd_response203_id,rpd_response204_id,
                          rpd_auc201_id,rpd_auc202_id,rpd_auc203_id,rpd_auc204_id)

##group differences
groupdiff_p<-c(NA,
               round(as.numeric(with(df_dem,chisq.test(t1_sex,t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,chisq.test(timepoints,t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(t1_ageyrs~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(t1_fsiq~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(t1_piq~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(t1_viq~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(t1_srs_rawscore~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(t1_rbs_total~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(t1_sdq_total_difficulties_p~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(adhd_inatt~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(adhd_hyper~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(anx_beck~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(dep_beck~t1_diagnosis))['p.value']),3),
               #round(as.numeric(with(df_dem,t.test(Precision~t1_diagnosis))['p.value']),3),
               #round(as.numeric(with(df_dem,t.test(Accuracy~t1_diagnosis))['p.value']),3),
               #round(as.numeric(with(df_dem,chisq.test(sampling_rate,t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(missing_data_trial~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_timepoint,chisq.test(sampling_rate,t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(center_dev~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(screen_dist~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(pd_baseline~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_response.201~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_response.202~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_response.203~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_response.204~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_auc.201~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_auc.202~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_auc.203~t1_diagnosis))['p.value']),3),
               round(as.numeric(with(df_dem,t.test(rpd_auc.204~t1_diagnosis))['p.value']),3)
               )

groupdiff_p[which(groupdiff_p==0)]<-'<0.001'
groupdiff_p[which(is.na(groupdiff_p))]<-'-'

row_names<-c('n','gender (F/M)','timepoints (t1/t2/t1+t2)','age (yrs)','IQ','perceptual IQ','verbal IQ',
             'SRS (total)','RBS (total)','SDQ (total)','ADHD inattention','ADHD hyperactivity','anxiety symptoms','depressive symptoms',
             #'data quality: precision','data quality: accuracy',
             #'sampling rate (120Hz/300Hz)',
             'missing data (%)','sampling rate (300Hz / 120Hz)',
             'center deviation (%)','screen distance (mm)','pupil size (mm)',
             'pupillary response: amplitude - standard (%)','pupillary response: amplitude - oddball pitch (%)','pupillary response: amplitude - oddball length (%)','pupillary response: amplitude - oddball both (%)',
             'pupillary response: AUC - standard (%)','pupillary response: AUC - oddball pitch (%)','pupillary response: AUC - oddball length (%)','pupillary response: AUC - oddball both (%)')

descriptives_table<-cbind(row_names,descriptives_table,groupdiff_p)

#add group differences
#descriptives_table<-ifelse(descriptives_table=='0','<0.001',descriptives_table)

require(kableExtra)
table_sample<-descriptives_table %>%
  kbl(caption = "Sample description",
      col.names = c('','autistic individuals (ASD)','neurotypical individuals (NTC)','group difference (p)'),
      row.names = F) %>%
  kable_classic(full_width = F, html_font = "Cambria")


descriptives_table

table_sample

###--> merge TRIAL data to demographics ####

#remove variables that are in both dfs
# exclusive_variables<-names(df_dem)[(names(df_dem) %in% names(df_trial))]
# exclusive_variables<-exclusive_variables[-1]
# df_dem_merge<-df_dem[,!names(df_dem) %in% exclusive_variables]
#
# df_trial<-merge(df_trial,df_dem_merge,by='subjects')

exclusive_variables<-names(df_timepoint)[(names(df_timepoint) %in% names(df_trial))]
exclusive_variables<-exclusive_variables[-1]
df_timepoint_merge<-df_timepoint[,!names(df_timepoint) %in% exclusive_variables]

df_trial<-merge(df_trial,df_timepoint_merge,by='id')

  table(is.na(df_timepoint$t1_piq))
  table(is.na(df_timepoint$sex))
  table(is.na(df_timepoint$ageyrs))


### trials in sequence - remove oddballs after oddballs ####

  ## --> retain for Bayesian modelling
  # with(df_trial,by(sequence_position,EventData,table))
  # ###--> remove oddball trials that are preceeded by oddballs
  # df_trial<-df_trial[(df_trial$sequence_position!=1 & df_trial$EventData!=201) | df_trial$EventData==201,]

###SAVE final df #####

  save(df, df_trial, df_timepoint, df_dem,file='C:/Users/nico/Desktop/mmn_leap_pd_final_dfs_21022023')

  ### MIXED MODELS ####
  with(df_dem,table(t1_schedule_adj,t1_diagnosis))

  hist(df_trial$rpd_response,50)
  hist(df_trial$rpd_auc,50)


  trials_by_event<-table(df_trial$id,df_trial$EventData)
      hist(trials_by_event[,1],20)
      hist(trials_by_event[,2],20)
      hist(trials_by_event[,3],20)
      hist(trials_by_event[,4],20)

      #trial per event
      apply(trials_by_event,2,mean)
      apply(trials_by_event,2,sd)

  table(df_trial$sampling_rate)
  hist(df_trial$center_dev)
  table(df_trial$subjects)
  table(df_trial$id)

  with(df_trial,hist(trial_duration[trial_duration<0.7],30))
  with(df_trial,hist(number_of_samples,30))

  ### - SEPR - PUPILLARY RESPONSE ####

  ###polynomial fit of trial counter
  lmm<-lmer(scale(rpd_auc)~EventData*scale(EventCounter)+
          (1|subjects)+(1|wave),data=df_trial)

  lmm2<-lmer(scale(rpd_auc)~EventData*poly(scale(EventCounter),2)+
              (1|subjects)+(1|wave),data=df_trial)

  lmm3<-lmer(scale(rpd_auc)~EventData*poly(scale(EventCounter),3)+
               (1|subjects)+(1|wave),data=df_trial)

  anova(lmm,lmm2,lmm3)
  ###>BIC does not improve - rather linear fit

  #pupillary response - technical model
  lmm<-lmer(scale(rpd_auc)~EventData+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)
  contrast(emmeans(lmm,~EventData),'revpairwise')


        #predicted values
        predicted_rpd<-predict(lmm)
        predicted_rpd_mean<-as.numeric(with(df_trial,by(predicted_rpd,interaction(id,EventData),median)))
        predicted_rpd_group<-as.factor(with(df_trial,by(t1_diagnosis,interaction(id,EventData),head,n=1)))
        predicted_rpd_eventdata<-as.factor(with(df_trial,by(EventData,interaction(id,EventData),head,n=1)))
        # levels(predicted_rpd_eventdata)<-c('low-utility','high-utility')
        # levels(predicted_RT_group)<-c('ADHD','ASD','TD')
        df_plot_predicted_rpd<-data.frame(predicted_rpd_group,predicted_rpd_eventdata,predicted_rpd_mean)


      plot_task_effect<-plot(emmeans(lmm,~EventData))[['data']]

            ggplot(plot_task_effect,aes(x=EventData,y=the.emmean))+geom_errorbar(aes(min=asymp.LCL,max=asymp.UCL))+
              geom_jitter(data=df_plot_predicted_rpd[complete.cases(df_plot_predicted_rpd),],
                          aes(x=predicted_rpd_eventdata,
                              y=predicted_rpd_mean,shape=predicted_rpd_group),alpha=0.3)+ #overplot predicted values

              theme_bw()


            ggplot(plot_task_effect,aes(x=EventData,y=the.emmean))+
              #conventional error bar
              geom_errorbar(aes(min=asymp.LCL,max=asymp.UCL),width=0.2)+
              #overplot predicted values
              geom_jitter(data=df_plot_predicted_rpd[complete.cases(df_plot_predicted_rpd),],
                          aes(x=predicted_rpd_eventdata,
                              y=predicted_rpd_mean,color=predicted_rpd_eventdata),alpha=0.3,width=0.4,show.legend=F)+
              #boxplot
              geom_boxplot(aes(fill=EventData,
                               middle=the.emmean,
                               lower=the.emmean-1.5*SE,
                               upper=the.emmean+1.5*SE,
                               ymin=asymp.LCL,
                               ymax=asymp.UCL),stat = "identity",alpha=0.7)+
              scale_fill_discrete(name = "stimulus", labels = c('standard','oddball pitch','oddball length','oddball both'))+
              scale_x_discrete(labels=c('standard','oddball pitch','oddball length','oddball both'))+

              theme_bw()


  table(df_trial$EventData)/sum(table(df_trial$EventData))
  #~oddball propability = 15%

  ## include trial characteristics
  lmm<-lmer(scale(rpd_auc)~EventData*(sequence_position+EventCounter)+
              scale(ageyrs)+scale(t1_piq)+sex+
              scale(center_dev)+
              as.factor(sampling_rate)+scale(missing_data_trial)+
              (1|subjects),data=df_trial)
  anova(lmm)
  # task effect: higher responses in oddballs
  contrast(emmeans(lmm,~EventData),'revpairwise')
  fixef(lmm)['sequence_position'] # later position associated with increased response

  #covariate effects
  fixef(lmm)['scale(missing_data_trial)']
  fixef(lmm)['scale(center_dev)']

  #making categorical variables
  hist(df_trial$sequence_position)
  summary(df_trial$sequence_position)
  df_trial$sequence_position_cat<-ifelse(df_trial$sequence_position>2,'unpredictable','predictable')
  ###--> 15% oddballs indicates 1 in 7 trials is an oddball

  df_trial$EventCounter_cat<-ifelse(df_trial$EventCounter<=700,'first half','second half')

  #pupillary response by group - random intercept
  lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*(sequence_position+EventCounter)+
              # scale(ageyrs)+scale(t1_piq)+sex+
              # scale(center_dev)+
              #as.factor(sampling_rate)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)
              #--> group effect (interaction) is robust in reduced and full model

  #pupillary response by group - random slope
  lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*EventCounter+
              # scale(ageyrs)+scale(t1_piq)+sex+
              # scale(center_dev)+
              #as.factor(sampling_rate)+scale(missing_data_trial)+
              (EventCounter-1|subjects)+(1|wave),data=df_trial)
  #--> group effect (interaction) is robust in reduced and full model
  anova(lmm)

        #Bayes factor based on BIC - Wagenmakers 2007 - http://www.ejwagenmakers.com/2007/pValueProblems.pdf
          ###approximation of Bayesian without definition of priors (i.e. flat priors)
        full_lmm<-lme4::lmer(scale(rpd_auc)~EventData*t1_diagnosis*EventCounter+
                               # scale(ageyrs)+scale(t1_piq)+sex+
                               # scale(center_dev)+
                               #as.factor(sampling_rate)+scale(missing_data_trial)+
                               (1|subjects)+(1|wave),data=df_trial,REML=F)
        anova(full_lmm)
        null_lmm <- update(full_lmm, formula = ~ . -t1_diagnosis:EventData:EventCounter)
        BF_BIC <- exp((BIC(full_lmm) - BIC(null_lmm))/2)  # BICs to Bayes factor
        BF_BIC
        #### high evidence --> however inclusion of covariates makes a big difference

        # #takes forever on windows - TRUE calculation of BayesFactor
        # require(BayesFactor)
        # df_trial$rpd_auc_z<-scale(df_trial$rpd_auc)
        # df_trial$subjects_factor<-as.factor(df_trial$subjects)
        # #simplified model
        # bfMain<-lmBF(rpd_auc_z ~ EventData*t1_diagnosis+subjects_factor, data = df_trial, whichRandom = "subjects_factor")

        # # Savage-Dickey desnity ratio
        # require(bayestestR)
        # require(logspline)
        # posterior<-predict(full_lmm)
        # prior<-predict(null_lmm)
        # bayesfactor(posterior, prior = prior)
        # #### --> provides very different BF, but may not be applicabel to LMM: https://pubmed.ncbi.nlm.nih.gov/30451277/#:~:text=The%20Savage%2DDickey%20density%20ratio%20is%20a%20simple%20method%20for,effect%20on%20the%20dependent%20variable.
        #

  qqnorm(residuals(lmm))
  hist(residuals(lmm),50)
  anova(lmm)

    #task effect
    plot(emmeans(lmm,~EventData)) ## --> higher response for all oddballs
    #plot(contrast(emmeans(lmm,~EventData|t1_diagnosis),'pairwise')) ## --> higher response for all oddballs
    #covariate effects
    round(fixef(lmm)['sequence_position'],5) #positive --> in Bayesian modelling relates to PREDICTABILITY/CERTAINTY
    round(fixef(lmm)['EventCounter'],5) #negative --> in Bayesian modelling - direction of precision weighting

      #cbind(round(fixef(lmm)['scale(t1_ageyrs)'],2),round(confint(lmm,parm = 'scale(t1_ageyrs)'),2))
      ### --> higher age with lower response: beta: -0.14 [-0.19; -0.10]

    #predictability effects
    plot(emmeans(lmm,~sequence_position_cat|EventData))
    confint(emmeans(lmm,~sequence_position_cat|EventData))
    contrast(emmeans(lmm,~EventData|sequence_position_cat),'pairwise')
    ###--> better differentiation of prediction errors under unpredictability

    #GROUP DIFFERENCES
      ##predictability effect
      contrast(emmeans(lmm,~sequence_position_cat|t1_diagnosis),'pairwise') ##--> in TD stronger response to pitch compared to length
      plot(emmeans(lmm,~sequence_position_cat|t1_diagnosis)) ##--> in TD stronger response to pitch compared to length
      ##arousal effect
      plot(emmeans(lmm,~EventData|t1_diagnosis)) ##--> in TD stronger response to pitch compared to length
      contrast(emtrends(lmm,~t1_diagnosis|EventData,var='EventCounter'),'pairwise')
      plot(emtrends(lmm,~t1_diagnosis|EventData,var='EventCounter'))
      ### --> attenuated habituation in ASD for pitch oddballs (only contrast significant in pairwise comparison)
      ##--> for rpd_auc higher pitch response in controls, higher oddball length response in ASD

    ###PLOT EFFECTS
    model_plot_data2<-as.data.frame(emtrends(lmm,~t1_diagnosis|EventData,var='EventCounter'))
    labels <- c("201" = "standard", "202" = "oddball pitch", "203" = "oddball length", "204" = "oddball both")
    g1<-ggplot(model_plot_data2,aes(x=t1_diagnosis,y=EventCounter.trend,group=t1_diagnosis,color=t1_diagnosis))+
      #geom_smooth(aes(ymin=asymp.LCL,ymax=asymp.UCL))+
      geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL),width=0.4)+
      geom_point()+
      geom_hline(yintercept=0,linetype=2)+
      facet_grid(~EventData,labeller=labeller(EventData = labels))+
      scale_color_manual(name = "group", labels = c("autistic", "neurotypical"),values = brewer.pal(3, "Dark2")[1:2])+
      labs(x='',y='effect of trial number on pupillary response (beta)')+
      theme_bw()

    model_plot_data<-as.data.frame(emmeans(lmm,~t1_diagnosis+EventData+EventCounter,at=list(EventCounter = seq(1,1400,70))))
    labels <- c("201" = "standard", "202" = "oddball pitch", "203" = "oddball length", "204" = "oddball both")
    g2<-ggplot(model_plot_data,aes(x=EventCounter,y=emmean,
                               group=interaction(EventData,t1_diagnosis),color=t1_diagnosis))+
      #geom_smooth(aes(ymin=asymp.LCL,ymax=asymp.UCL))+
      geom_errorbar(aes(ymin=emmean-1*SE,ymax=emmean+1*SE))+
      geom_point()+
      geom_line()+
      geom_hline(yintercept=0,linetype=2)+
      facet_wrap(~EventData,labeller=labeller(EventData = labels))+
      scale_color_manual(name = "group", labels = c("autistic", "neurotypical"),values = brewer.pal(3, "Dark2")[1:2])+
      labs(x='trial number',y='pupillary response (z)')+
      theme_bw()+
      theme(legend.position = "none")

    grid.arrange(g1,g2)


    ###model comparison
    lmm2<-lmer(scale(rpd_auc)~EventData*(sequence_position_cat+EventCounter)+
                # scale(ageyrs)+scale(t1_piq)+sex+
                # scale(center_dev)+
                as.factor(sampling_rate)+scale(missing_data_trial)+
                (1|subjects)+(1|wave),data=df_trial)

      anova(lmm,lmm2)
      ###-->not a better model


          #drop all with CSS < 4
          df_trial_excl<-df_trial[df_trial$t1_diagnosis == 'Control' | (!is.na(df_trial$t1_css_total_all) & df_trial$t1_css_total_all>3),]
          lmm<-lmer(scale(rpd_response)~EventData*t1_diagnosis*(sequence_position+EventCounter)+
                      # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                      # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                      (1|subjects),data=df_trial_excl)
          anova(lmm)

          #justification  for baseline pd as covariate - but has a group difference
          lmm<-lmer(scale(rpd_response)~EventData*(scale(pd_baseline)+sequence_position+EventCounter)+
                      # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                      # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                      (1|subjects)+(1|wave),data=df_trial)


          anova(lmm)
          fixef(lmm)
          ###--> baseline pd
          cbind(round(fixef(lmm)['scale(pd_baseline)'],2),
                round(confint(lmm,parm = 'scale(pd_baseline)'),2))
          # beta = - 0.43 [-0.44; -0.42]


      ###effects in age subgroups

        #CHILDREN
      df_trial_children<-df_trial[df_trial$t1_schedule_adj=='Children',]
      lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects)+(1|wave),data=df_trial_children)


      anova(lmm)
      #task effect
      plot(emmeans(lmm,~EventData)) ###--> response attenuated for length & combined oddballs
      plot(emmeans(lmm,~t1_diagnosis|EventData)) ###-->


      #ADOLESCENTS
      df_trial_adol<-df_trial[df_trial$t1_schedule_adj=='Adolescents',]
      lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects)+(1|wave),data=df_trial_adol)


      anova(lmm)
      #task effect
      plot(emmeans(lmm,~EventData)) ###--> response attenuated for length oddballs
      #group difference
      plot(emtrends(lmm,~t1_diagnosis|EventData,var='EventCounter'))


      #ADULTS
      df_trial_adul<-df_trial[df_trial$t1_schedule_adj=='Adults',]
      lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects)+(1|wave),data=df_trial_adul)

      anova(lmm)
      #task effect
      plot(emmeans(lmm,~EventData)) ###--> response FOR all Oddballs


  ### stability of respones #####

      table(df_trial$wave)
      table(df_dem$timepoints)

      #select only individuals that have both waves
      subjects_with_both_waves<-df_dem$subjects[df_dem$timepoints=='wave1+2']
      df_trial_wave12<-df_trial[df_trial$subjects %in% subjects_with_both_waves,]

      #response to standards
      lmm<-lmer(scale(rpd_auc.201)~wave+EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects),data=df_trial_wave12)
      anova(lmm)
      contrast(emmeans(lmm,~wave),'pairwise')
      ###--> lower response to standards in wave 2

      #response to pitch oddballs
      lmm<-lmer(scale(rpd_auc.202)~wave+EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects),data=df_trial_wave12)
      anova(lmm)
      contrast(emmeans(lmm,~wave),'pairwise')
      ###--> higher response to pitches in wave 2

      #response to length oddballs
      lmm<-lmer(scale(rpd_auc.203)~wave+EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects),data=df_trial_wave12)
      anova(lmm)
      contrast(emmeans(lmm,~wave),'pairwise')
      ###--> lower response to length oddballs in wave 2

      #response to length oddballs
      lmm<-lmer(scale(rpd_auc.204)~wave+EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects),data=df_trial_wave12)
      anova(lmm)
      contrast(emmeans(lmm,~wave),'pairwise')
      ###--> lower response to length+pitch oddballs in wave 2




      lmm<-lmer(scale(pd_baseline)~wave*EventData*EventCounter+
                  # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
                  # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
                  (1|subjects),data=df_trial_wave12)
      anova(lmm)

      cbind(round(fixef(lmm)['wavewave2'],2),
            round(confint(lmm,parm = 'wavewave2'),2))
      ###reduction in baseline PD across waves


    #### --> ICC across waves ####

      table(df_trial_wave12$wave)

      #BPS
      baseline_pd_wave1<-with(df_trial_wave12[df_trial_wave12$wave=='wave1',],aggregate(pd_baseline,by=list(subjects),FUN=mean,na.rm=T))
      baseline_pd_wave2<-with(df_trial_wave12[df_trial_wave12$wave=='wave2',],aggregate(pd_baseline,by=list(subjects),FUN=mean,na.rm=T))
      names(baseline_pd_wave1)<-c('subjects','wave1')
      names(baseline_pd_wave2)<-c('subjects','wave2')

      baseline_pd_waves<-merge(baseline_pd_wave1,baseline_pd_wave2,by='subjects')

      cor.test(baseline_pd_waves$wave1,baseline_pd_waves$wave2)

      psych::ICC(baseline_pd_waves[,-1])
      ###--> ICC: 0.88 [0.81, 0.93]


      #SEPR
      names(df_trial_wave12)

      rpd204_wave1<-with(df_trial_wave12[df_trial_wave12$EventData=='204' & df_trial_wave12$wave=='wave1',],aggregate(rpd_auc.204,by=list(subjects),FUN=mean,na.rm=T))
      rpd204_wave2<-with(df_trial_wave12[df_trial_wave12$EventData=='204' & df_trial_wave12$wave=='wave2',],aggregate(rpd_auc.204,by=list(subjects),FUN=mean,na.rm=T))

      names(rpd204_wave1)<-c('subjects','wave1')
      names(rpd204_wave2)<-c('subjects','wave2')

      rpd_waves<-merge(rpd204_wave1,rpd204_wave2,by='subjects')

      cor.test(rpd_waves$wave1,rpd_waves$wave2)
      ##--> r =.23

      psych::ICC(rpd_waves[,-1])
      ###--> ICC: 0.37 [-0.01, 0.60]




  ### - BPS - baseline progression ####
  #lmm_data<-df_trial[df_trial$EventData=='201' & !is.na(missing_data_trial) & !is.na(center_dev) & !is.na(pd_baseline) & !is.na(t1_piq),]

  #response as a function of baseline
  lmm<-lmer(scale(rpd_auc)~EventData*scale(pd_baseline)*(sequence_position+EventCounter)+
              # scale(ageyrs)+scale(t1_piq)+sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)

      anova(lmm)
      fixef(lmm)['scale(pd_baseline)']
      ###-> ###main effect beta=-0.38


  lmm<-lmer(scale(pd_baseline)~EventData*t1_diagnosis*(sequence_position+EventCounter)+
              # scale(ageyrs)+scale(t1_piq)+sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)

  fixef(lmm)['EventCounter']
  emtrends(lmm,~t1_diagnosis,var='EventCounter')
  emtrends(lmm,~t1_diagnosis,var='sequence_position')

        #Bayes factor based on BIC - Wagenmakers 2007 - http://www.ejwagenmakers.com/2007/pValueProblems.pdf
        ###approximation of Bayesian without definition of priors (i.e. flat priors)
        full_lmm<-lme4::lmer(scale(pd_baseline)~EventData*t1_diagnosis*(sequence_position+EventCounter)+
                               # scale(ageyrs)+scale(t1_piq)+sex+
                               # scale(center_dev)+
                               #as.factor(sampling_rate)+scale(missing_data_trial)+
                               (1|subjects)+(1|wave),data=df_trial,REML=F)
        anova(full_lmm)
        null_lmm <- update(full_lmm, formula = ~ . -t1_diagnosis:EventData:EventCounter)
        BF_BIC <- exp((BIC(full_lmm) - BIC(null_lmm))/2)  # BICs to Bayes factor
        BF_BIC
        #### high evidence --> however inclusion of covariates makes a big difference



  #main effects
  fixef(lmm)['EventCounter']*1000 #across 1000 trials a reduction of 1.7% SD
  cbind(round(fixef(lmm)['EventCounter'],5),
        round(confint(lmm,parm = 'EventCounter'),5))


  #covariate effects
  fixef(lmm)['scale(ageyrs)']
    cbind(round(fixef(lmm)['scale(ageyrs)'],3),
          round(confint(lmm,parm = 'scale(ageyrs)'),3))
    ###--> higher age with lower pd_baseline: -0.56 [-0.63; -0.48]

  #task interactions
    plot(emtrends(lmm,~EventData,var='EventCounter'))
    ###--> baseline reduction is attenuated in 204
    plot(emtrends(lmm,~EventData,var='sequence_position'))

  #group effects
  emtrends(lmm,~t1_diagnosis,var='sequence_position')
  emtrends(lmm,~t1_diagnosis,var='EventCounter')
  ###--> attenuated habituation in ASD


  #plot interaction
  model_plot_data<-as.data.frame(emmeans(lmm,~t1_diagnosis+sequence_position,
                                         at=list(sequence_position = seq(1,25,1))))

  g1<-ggplot(model_plot_data,aes(x=sequence_position,y=emmean,group=t1_diagnosis,color=t1_diagnosis))+
    #geom_smooth(aes(ymin=asymp.LCL,ymax=asymp.UCL))+
    geom_errorbar(aes(ymin=emmean-1*SE,ymax=emmean+1*SE))+
    geom_line()+
    geom_point()+
    scale_color_manual(name = "group", labels = c("autistic", "neurotypical"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='position in sequence',y='effect on baseline pupil size (z)')+
    theme_bw()


  #plot interaction
  model_plot_data<-as.data.frame(emmeans(lmm,~t1_diagnosis+EventCounter,
                                         at=list(EventCounter = seq(1,1400,70))))

  g2<-ggplot(model_plot_data,aes(x=EventCounter,y=emmean,group=t1_diagnosis,color=t1_diagnosis))+
    #geom_smooth(aes(ymin=asymp.LCL,ymax=asymp.UCL))+
    geom_errorbar(aes(ymin=emmean-1*SE,ymax=emmean+1*SE))+
    geom_line()+
    geom_point()+
    scale_color_manual(name = "group", labels = c("autistic", "neurotypical"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='total number of trials',y='effect on baseline pupil size (z)')+
    theme_bw()

  grid.arrange(g1,g2)

  ###plot observed pupil size
  require(ggplot2)
  ggplot(df_trial,aes(x=EventCounter,y=scale(pd_baseline),group=EventData,color=EventData))+geom_smooth(method='lm')+theme_bw()


  ### - clinical associations ####

  ###---> SEPR ####
  df_trial$ssp_total_scaled<-scale(df_trial$ssp_total)
  lmm<-lmer(scale(rpd_auc)~EventData*ssp_total_scaled*EventCounter+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)
  fixef(lmm)

  emtrends(lmm,~EventData|EventCounter,var='ssp_total_scaled',at=list(EventCounter = seq(1,1400,140)))
  ####--> higher SSP is associated increased reaction to oddball in the first half of the task

  names(df_trial)

  df_trial$anx_beck<-with(df_trial,ifelse(!is.na(t1_beck_anx_adulta_self),t1_beck_anx_adulta_self,
                               ifelse(!is.na(t1_beck_anx_youthb_self),t1_beck_anx_youthb_self,t1_beck_anx_youthcd_parent)))
  df_trial$anx_beck_scaled<-scale(df_trial$anx_beck)

  lmm<-lmer(scale(rpd_auc)~EventData*anx_beck_scaled*EventCounter+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)
  ###--> no effect


  df_trial$dep_beck<-with(df_trial,ifelse(!is.na(t1_beck_dep_adulta_self),t1_beck_dep_adulta_self,
                                          ifelse(!is.na(t1_beck_dep_youthb),t1_beck_dep_youthb,t1_beck_dep_adultd_parent)))
  df_trial$dep_beck_scaled<-scale(df_trial$dep_beck)

  lmm<-lmer(scale(rpd_auc)~EventData*dep_beck_scaled*EventCounter+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)
  ###--> no effect

  names(df_trial)

  ##### -- transdiagnostic factors ####

  names(df_trial)
  table(is.na(df_timepoint$sdq_internalising_p))

  df_trial$sdq_internalising_p_scaled<-scale(df_trial$sdq_internalising_p)
  df_trial$sdq_externalising_p_scaled<-scale(df_trial$sdq_externalising_p)
  df_trial$sdq_total_p_scaled<-scale(df_trial$sdq_total_difficulties_p)

  hist(df_trial$sdq_internalising_p_scaled,breaks=10)
  hist(df_trial$sdq_externalising_p_scaled,breaks=10)
  hist(df_trial$sdq_total_p_scaled,breaks=10)


  lmm<-lmer(scale(rpd_auc)~EventData*sdq_internalising_p_scaled*EventCounter+
              (1|subjects)+(1|wave),data=df_trial)
  anova(lmm)
  #no effect

        ###in ASD
        lmm<-lmer(scale(rpd_auc)~EventData*sdq_internalising_p_scaled*EventCounter+
                    (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='ASD',])
        anova(lmm)
        #no effect

        ###in TD
        lmm<-lmer(scale(rpd_auc)~EventData*sdq_internalising_p_scaled*EventCounter+
                    (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='Control',])
        anova(lmm)
        #no effect


  lmm<-lmer(scale(rpd_auc)~EventData*sdq_externalising_p_scaled*EventCounter+
              (1|subjects)+(1|wave),data=df_trial)
  anova(lmm)
  #no effect

        ##in ASD
        lmm<-lmer(scale(rpd_auc)~EventData*sdq_externalising_p_scaled*EventCounter+
                    (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='ASD',])
        anova(lmm)
        #no effectr

        ###in TD ~ n = 31
        lmm<-lmer(scale(rpd_auc)~EventData*sdq_externalising_p_scaled*EventCounter+
                    (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='Control',])
        anova(lmm)
        #no effect


  lmm<-lmer(scale(rpd_auc)~EventData*sdq_total_p_scaled*EventCounter+
              (1|subjects)+(1|wave),data=df_trial)
  anova(lmm)
  #-> no effect


        ##in ASD
        lmm<-lmer(scale(rpd_auc)~EventData*sdq_total_p_scaled*EventCounter+
                    (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='ASD',])
        anova(lmm)
        ##--> no effect

        ###in TD ~ n =31
        lmm<-lmer(scale(rpd_auc)~EventData*sdq_total_p_scaled*EventCounter+
                    (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='Control',])
        anova(lmm)
        ## --> no effect



  names(df_trial)[grepl('aq',names(df_trial))]
  df_trial$aq_total<-with(df_trial,ifelse(!is.na(t1_aq_adult_total),t1_aq_adult_total,
                                          ifelse(!is.na(t1_aq_adol_total),t1_aq_adol_total,
                                                 ifelse(!is.na(t1_aq_child_total),t1_aq_child_total,NA))))

  df_trial$aq_total_z<-scale(df_trial$aq_total)
  lmm<-lmer(scale(rpd_response)~EventData*aq_total_z*(sequence_position+EventCounter)+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects),data=df_trial)

  anova(lmm)
  ###--> no effect


  names(df_trial)[grepl('srs',names(df_trial))]
  df_trial$srs_total_z<-scale(df_trial$srs_rawscore_combined)
  lmm<-lmer(scale(rpd_response)~EventData*srs_total_z*(sequence_position+EventCounter)+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects),data=df_trial[df_trial$t1_diagnosis=='Control',])

  anova(lmm)
  ###--> no
  plot(emtrends(lmm,~EventData,var='srs_total_z'))
  ####--> in TD: the response to pitch oddball scales with SRS total score






  names(df_trial)[grepl('aq',names(df_trial))]
  df_trial$aq_total<-with(df_trial,ifelse(!is.na(t1_aq_adult_total),t1_aq_adult_total,
                                          ifelse(!is.na(t1_aq_adol_total),t1_aq_adol_total,
                                                 ifelse(!is.na(t1_aq_child_total),t1_aq_child_total,NA))))

  df_trial$aq_total_z<-scale(df_trial$aq_total)
  lmm<-lmer(scale(rpd_auc)~EventData*aq_total_z*(sequence_position+EventCounter)+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects),data=df_trial)

  anova(lmm)
  ###--> no effect


  names(df_trial)[grepl('srs',names(df_trial))]
  df_trial$srs_total_z<-scale(df_trial$srs_rawscore_combined)
  lmm<-lmer(scale(rpd_auc)~EventData*srs_total_z*(sequence_position+EventCounter)+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects),data=df_trial[df_trial$t1_diagnosis=='Control',])

  anova(lmm)
  ###--> no
  plot(emtrends(lmm,~EventData,var='srs_total_z'))
  ####--> in TD: the response to pitch oddball scales with SRS total score


  ###--> BPS ####
  df_trial$ssp_total_scaled<-scale(df_trial$ssp_total)
  lmm<-lmer(scale(pd_baseline)~EventData*ssp_total_scaled*EventCounter+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)

  cbind(round(fixef(lmm)['ssp_total_scaled'],2),
        round(confint(lmm,parm = 'ssp_total_scaled'),2))

  ### higher SSP with higher pd_baseline (beta = 0.21 [0.20, 0.24])
  emtrends(lmm,~EventCounter,var='ssp_total_scaled',at=list(EventCounter = seq(1,1400,140)))
  ####--> this association is stronger int he beginning

          #in ASD and TD
          lmm<-lmer(scale(pd_baseline)~EventData*ssp_total_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='ASD',])

          anova(lmm)
          cbind(round(fixef(lmm)['ssp_total_scaled'],2),
                round(confint(lmm,parm = 'ssp_total_scaled'),2))
          #-> ASD. beta = 0.29 [0.27, 0.31]
          emtrends(lmm,~EventCounter,var='ssp_total_scaled',at=list(EventCounter = seq(1,1400,140)))
          ####--> this association is weaker at the end of the task


          #in ASD and TD
          lmm<-lmer(scale(pd_baseline)~EventData*ssp_total_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='Control',])

          anova(lmm)
          cbind(round(fixef(lmm)['ssp_total_scaled'],2),
                round(confint(lmm,parm = 'ssp_total_scaled'),2))
          #-> TD beta = -0.38 [-0.42, -0.33]
          emtrends(lmm,~EventCounter,var='ssp_total_scaled',at=list(EventCounter = seq(1,1400,140)))
          ####--> this association is stronger at the end of the task


  names(df_trial)

  df_trial$anx_beck<-with(df_trial,ifelse(!is.na(t1_beck_anx_adulta_self),t1_beck_anx_adulta_self,
                                          ifelse(!is.na(t1_beck_anx_youthb_self),t1_beck_anx_youthb_self,t1_beck_anx_youthcd_parent)))
  df_trial$anx_beck_scaled<-scale(df_trial$anx_beck)

  lmm<-lmer(scale(pd_baseline)~EventData*anx_beck_scaled*EventCounter+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)
  ###--> no effect


  df_trial$dep_beck<-with(df_trial,ifelse(!is.na(t1_beck_dep_adulta_self),t1_beck_dep_adulta_self,
                                          ifelse(!is.na(t1_beck_dep_youthb),t1_beck_dep_youthb,t1_beck_dep_adultd_parent)))
  df_trial$dep_beck_scaled<-scale(df_trial$dep_beck)

  lmm<-lmer(scale(pd_baseline)~EventData*dep_beck_scaled*EventCounter+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)

  anova(lmm)
  ###--> no effect

  ##### -- transdiagnostic factors ####

  names(df_trial)
  table(is.na(df_timepoint$sdq_internalising_p))

  df_trial$sdq_internalising_p_scaled<-scale(df_trial$sdq_internalising_p)
  df_trial$sdq_externalising_p_scaled<-scale(df_trial$sdq_externalising_p)
  df_trial$sdq_total_p_scaled<-scale(df_trial$sdq_total_difficulties_p)

  hist(df_trial$sdq_internalising_p_scaled,breaks=10)
  hist(df_trial$sdq_externalising_p_scaled,breaks=10)
  hist(df_trial$sdq_total_p_scaled,breaks=10)


  lmm<-lmer(scale(pd_baseline)~EventData*sdq_internalising_p_scaled*EventCounter+
              # scale(t1_ageyrs)+scale(t1_piq)+t1_sex+
              # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
              (1|subjects)+(1|wave),data=df_trial)
  anova(lmm)

  cbind(round(fixef(lmm)['sdq_internalising_p_scaled'],2),
        round(confint(lmm,parm = 'sdq_internalising_p_scaled'),2))
  #-> beta = 0.07 [0.05, 0.08]
  emtrends(lmm,~EventCounter,var='sdq_internalising_p_scaled',at=list(EventCounter = seq(1,1400,140)))
  ####--> this association increases with task progression

          ###in ASD
          lmm<-lmer(scale(pd_baseline)~EventData*sdq_internalising_p_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='ASD',])
          anova(lmm)
          cbind(round(fixef(lmm)['sdq_internalising_p_scaled'],2),
                round(confint(lmm,parm = 'sdq_internalising_p_scaled'),2))
          #-> beta = 0.11 [0.09, 0.13]

          ###in TD
          lmm<-lmer(scale(pd_baseline)~EventData*sdq_internalising_p_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='Control',])
          anova(lmm)
          cbind(round(fixef(lmm)['sdq_internalising_p_scaled'],2),
                round(confint(lmm,parm = 'sdq_internalising_p_scaled'),2))
          #-> beta = -0.15 [-0.18, -0.12]


  lmm<-lmer(scale(pd_baseline)~EventData*sdq_externalising_p_scaled*EventCounter+
              (1|subjects)+(1|wave),data=df_trial)
  anova(lmm)

  cbind(round(fixef(lmm)['sdq_externalising_p_scaled'],2),
        round(confint(lmm,parm = 'sdq_externalising_p_scaled'),2))
  #-> beta = -0.17 [-0.18, -0.15]


          ##in ASD
          lmm<-lmer(scale(pd_baseline)~EventData*sdq_externalising_p_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='ASD',])
          anova(lmm)
          cbind(round(fixef(lmm)['sdq_externalising_p_scaled'],2),
                round(confint(lmm,parm = 'sdq_externalising_p_scaled'),2))
          #-> beta = -0.18 [-0.20, -0.16]

          ###in TD ~ n = 31
          lmm<-lmer(scale(pd_baseline)~EventData*sdq_externalising_p_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='Control',])
          anova(lmm)
          cbind(round(fixef(lmm)['sdq_externalising_p_scaled'],2),
                round(confint(lmm,parm = 'sdq_externalising_p_scaled'),2))
          #-> beta = -0.07 [-0.12, -0.01]


  lmm<-lmer(scale(pd_baseline)~EventData*sdq_total_p_scaled*EventCounter+
              (1|subjects)+(1|wave),data=df_trial)
  anova(lmm)
  #-> no effect


          ##in ASD
          lmm<-lmer(scale(pd_baseline)~EventData*sdq_total_p_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='ASD',])
          anova(lmm)
          ##--> no effect

          ###in TD ~ n =31
          lmm<-lmer(scale(pd_baseline)~EventData*sdq_total_p_scaled*EventCounter+
                      (1|subjects)+(1|wave),data=df_trial[df_trial$t1_diagnosis=='Control',])
          anova(lmm)
          cbind(round(fixef(lmm)['sdq_total_p_scaled'],2),
                round(confint(lmm,parm = 'sdq_total_p_scaled'),2))
          #-> beta = -0.16 [-0.19, -0.12]

    ### transdiagnosic associations on per participant level

          #BPS
          with(df_timepoint,cor.test(scale(sdq_internalising_p),pd_baseline))
          with(df_timepoint,cor.test(scale(sdq_externalising_p),pd_baseline))
          with(df_timepoint,cor.test(scale(sdq_total_difficulties_p),pd_baseline))
          ### --> all weak positive associations across groups
          with(df_timepoint[df_timepoint$t1_diagnosis=='ASD',],cor.test(scale(sdq_internalising_p),pd_baseline))
          with(df_timepoint[df_timepoint$t1_diagnosis=='ASD',],cor.test(scale(sdq_externalising_p),pd_baseline))
          with(df_timepoint[df_timepoint$t1_diagnosis=='ASD',],cor.test(scale(sdq_total_difficulties_p),pd_baseline))
          ### --> not significant in ASD --> power problem
          with(df_timepoint[df_timepoint$t1_diagnosis=='Control',],cor.test(scale(sdq_internalising_p),pd_baseline))
          with(df_timepoint[df_timepoint$t1_diagnosis=='Control',],cor.test(scale(sdq_externalising_p),pd_baseline))
          with(df_timepoint[df_timepoint$t1_diagnosis=='Control',],cor.test(scale(sdq_total_difficulties_p),pd_baseline))
          ### --> not significant in TD --> power problem

          #SEPR ---> no association
          with(df_timepoint,cor.test(scale(sdq_internalising_p),rpd_auc))
          with(df_timepoint,cor.test(scale(sdq_externalising_p),rpd_auc))
          with(df_timepoint,cor.test(scale(sdq_total_difficulties_p),rpd_auc))
          ### ---> no association
          with(df_timepoint[df_timepoint$t1_diagnosis=='ASD',],cor.test(scale(sdq_internalising_p),rpd_auc))
          with(df_timepoint[df_timepoint$t1_diagnosis=='ASD',],cor.test(scale(sdq_externalising_p),rpd_auc))
          with(df_timepoint[df_timepoint$t1_diagnosis=='ASD',],cor.test(scale(sdq_total_difficulties_p),rpd_auc))
          ### --> not significant in ASD --> power problem
          with(df_timepoint[df_timepoint$t1_diagnosis=='Control',],cor.test(scale(sdq_internalising_p),rpd_auc))
          with(df_timepoint[df_timepoint$t1_diagnosis=='Control',],cor.test(scale(sdq_externalising_p),rpd_auc))
          with(df_timepoint[df_timepoint$t1_diagnosis=='Control',],cor.test(scale(sdq_total_difficulties_p),rpd_auc))
          ### --> not significant in TD --> power problem
          with(df_timepoint[df_timepoint$t1_diagnosis=='Control',],cor.test(scale(sdq_externalising_p),rpd_auc))
          ###--> higher rpd_auc assoicated with higher externalising in TD

  ### - Changes in BPS and SEPR by sequence position (habituation / repetition suppression) ####

  names(df_trial)
  with(df_trial,table(sequence_position,EventData))

  #distribution of sequence position
  labels <- c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball", "204" = "pitch & length oddball")
  ggplot(df_trial,aes(x=sequence_position,group=EventData,fill=EventData))+geom_density(adjust=3)+facet_wrap(~EventData,labeller=labeller(EventData = labels))

  table_sequence_position<-with(df_trial,table(sequence_position,EventData))
  table_sample<-table_sequence_position %>%
    kbl(caption = "trials by sequence position",
        col.names = c('standard','pitch','length','pitch & length'),
        row.names = T) %>%
    kable_classic(full_width = F, html_font = "Cambria")
  save_kable(table_sample, file=paste0(project_path,'/output/table_trials_by_sequence_position.html'))

  table_sample


  #required transformations
  df_trial$standardtrial<-with(df_trial,ifelse(EventData==201,T,F))
  df_trial$sequence_position_z<-scale(df_trial$sequence_position)

      ggplot(df_trial,aes(x=EventCounter,y=scale(pd_baseline),group=interaction(standardtrial,t1_diagnosis),linetype=t1_diagnosis,color=standardtrial,fill=standardtrial))+geom_smooth()+
        theme_bw()+labs(y='pupillary response: SEPR (z)')



  #visualization
  require(wesanderson)
  custom_condition_colors <- wes_palette('FantasticFox1',5,type='discrete')[2:5]

  hist(scale(df_trial$rpd_auc))
  hist(scale(df_trial$pd_baseline))

  ggplot(df_trial,aes(x=EventCounter,y=scale(rpd_response),group=EventData,color=EventData,fill=EventData))+geom_smooth()+
    theme_bw()+labs(y='pupillary response: SEPR (z)')+
    scale_fill_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))+
    scale_color_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))

  ggplot(df_trial,aes(x=EventCounter,y=scale(pd_baseline),group=EventData,color=EventData,fill=EventData))+geom_smooth()+
    theme_bw()+labs(y='baseline pupil size: BPS (z)')+
    scale_fill_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))+
    scale_color_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))

  ggplot(df_trial[df_trial$sequence_position<=10,],aes(x=sequence_position,y=scale(rpd_response),group=EventData,color=EventData,fill=EventData))+geom_smooth()+
    theme_bw()+labs(y='pupillary response: SEPR (z)')+
    scale_fill_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))+
    scale_color_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))

  ggplot(df_trial[df_trial$sequence_position<=10,],aes(x=sequence_position,y=scale(pd_baseline),group=EventData,color=EventData,fill=EventData))+geom_smooth()+
    theme_bw()+labs(y='baseline pupil size: BPS (z)')+
    scale_fill_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))+
    scale_color_manual(values = custom_condition_colors, labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))



  ggplot(df_trial[df_trial$sequence_position<4,],aes(x=interaction(as.factor(sequence_position),EventData),y=scale(rpd_response),color=EventData,fill=EventData))+geom_violin()+
    theme_bw()+labs(y='pupillary response: SEPR (z)')



  ggplot(df_trial[df_trial$sequence_position<=10,],aes(x=sequence_position,y=scale(pd_baseline),group=EventData,color=EventData))+geom_smooth()+theme_bw()


  #figure
  ggplot(df_trial[df_trial$sequence_position<=10,],aes(x=sequence_position,y=scale(pd_baseline),group=interaction(standardtrial,t1_diagnosis),color=standardtrial,linetype=t1_diagnosis))+
    geom_smooth()+labs(x='sequence position',y='BPS (z)')+scale_x_continuous(breaks=1:10)+theme_bw()
  ggplot(df_trial[df_trial$sequence_position<=10,],aes(x=sequence_position,y=scale(rpd_response),group=interaction(standardtrial,t1_diagnosis),color=t1_diagnosis,linetype=standardtrial))+
    geom_smooth()+labs(x='sequence position',y='SEPR (z)')+scale_x_continuous(breaks=1:10)+theme_bw()


  #BPS - continuous
  lmm<-lmer(scale(pd_baseline)~EventData*t1_diagnosis*sequence_position_z+
              (1|subjects)+(1|wave),data=df_trial[df_trial$sequence_position<10,])

  #model fit
  anova(lmm)
  r2_nakagawa(lmm)

        ### supplementary analysis
        lmm<-lmer(scale(pd_baseline)~t1_diagnosis*sequence_position_z+
                    (1|subjects)+(1|wave),data=df_trial[df_trial$EventData==201 & df_trial$sequence_position<10,])

        #model fit
        anova(lmm)


        #BIC Bayes Factor
        lmm_ML<-lmer(scale(pd_baseline)~t1_diagnosis*sequence_position_z+
                       (1|subjects)+(1|wave),data=df_trial[df_trial$EventData==201 & df_trial$sequence_position<10,],REML=F)

        null_lmm_ML <- lmer(scale(pd_baseline)~sequence_position_z+
                              (1|subjects)+(1|wave),data=df_trial[df_trial$EventData==201 & df_trial$sequence_position<10,],REML=F)

        BF_BIC <- exp((BIC(null_lmm_ML) - BIC(lmm_ML))/2)
        BF_BIC

  #BPS categorical
  lmm<-lmer(scale(pd_baseline)~EventData*t1_diagnosis*as.factor(sequence_position)+
              (1|subjects)+(1|wave),data=df_trial[df_trial$sequence_position<4,])

  anova(lmm)
  contrast(emmeans(lmm,~as.factor(sequence_position)),'pairwise')
  ###-> 1-3 effect
  ###inclusion of covariates leads to boundary singular

  lmm<-lmer(scale(pd_baseline)~t1_diagnosis*as.factor(sequence_position)+
              (1|subjects)+(1|wave),data=df_trial[df_trial$EventData==201 & df_trial$sequence_position<4,])



  #SEPR
  lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*sequence_position_z+
              (1|subjects)+(1|wave),data=df_trial[df_trial$sequence_position<10,])

  anova(lmm)


  lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*as.factor(sequence_position)+
              (1|subjects)+(1|wave),data=df_trial[df_trial$sequence_position<4,])


      #BIC Bayes Factor
      lmm_ML<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*as.factor(sequence_position)+
                     (1|subjects)+(1|wave),data=df_trial[df_trial$sequence_position<4,],REML=F)

      null_lmm_ML <- lmer(scale(rpd_auc)~1+
                            (1|subjects)+(1|wave),data=df_trial[df_trial$sequence_position<4,],REML=F)

      BF_BIC <- exp((BIC(null_lmm_ML) - BIC(lmm_ML))/2)
      BF_BIC

      BIC(null_lmm_ML)
      BIC(lmm_ML)



  anova(lmm)
  ###--> if an oddballs follows an oddball - higehr response in ASD
  ###--> if oddball follows a standard - higher response in TD

  contrast(emmeans(lmm,~t1_diagnosis|EventData+as.factor(sequence_position)),'pairwise')
  confint(contrast(emmeans(lmm,~t1_diagnosis|EventData+as.factor(sequence_position)),'revpairwise'))






  df_trial$trial_position<-ifelse(df_trial$EventCounter<=280,'first quintile',
                              ifelse(df_trial$EventCounter>1120,'fifth quintile',
                                     ifelse(df_trial$EventCounter>560 & df_trial$EventCounter<=840,'third quintile',
                                            ifelse(df_trial$EventCounter>280 & df_trial$EventCounter<=560,'second quintile',
                                                   ifelse(df_trial$EventCounter>840 & df_trial$EventCounter<=1120,'fourth quintile',NA)))))
  df_trial$trial_position<-factor(df_trial$trial_position,levels=c('first quintile','second quintile','third quintile','fourth quintile','fifth quintile'))


  lmm<-lmer(scale(rpd_auc)~t1_diagnosis*trial_position*as.factor(sequence_position)+
              (1|subjects)+(1|wave),data=df_trial[df_trial$EventData==201 & df_trial$sequence_position<4,])

  anova(lmm)
  contrast(emmeans(lmm,~as.factor(sequence_position)),'pairwise')
  contrast(emmeans(lmm,~t1_diagnosis|trial_position),'pairwise')



###--> ADDITIONAL VISUALIZATION: ####

  ##pupil signal in per trial data
  hist(df_trial$rpd_auc)
  ggplot(df_trial[df_trial$rpd_auc>-0.1 & df_trial$rpd_auc<0.1 ],aes(x=scale(rpd_auc),group=EventData))+geom_density()


  #pd baseline by age
  ggplot(df_trial,aes(x=ageyrs,y=pd_baseline,group=t1_diagnosis,color=t1_diagnosis))+geom_smooth(method='lm',formula = y ~ x + poly(x,3))+theme_bw()+geom_point(data=df_trial[sample(1:nrow(df_trial),nrow(df_trial)/100),],alpha=0.2)
  ggplot(df_timepoint,aes(x=ageyrs,y=pd_baseline,group=t1_diagnosis,color=t1_diagnosis))+geom_smooth(method='lm',formula = y ~ x + poly(x,3))+theme_bw()+geom_point(data=df_trial[sample(1:nrow(df_trial),nrow(df_trial)/100),],alpha=0.2)

  ggplot(df_trial,aes(x=ageyrs,y=rpd_auc,group=t1_diagnosis,color=t1_diagnosis))+geom_smooth(method='lm',formula = y ~ x + poly(x,3))+theme_bw()+geom_point(data=df_trial[sample(1:nrow(df_trial),nrow(df_trial)/100),],alpha=0.2)


  ###SCALING of pupillary response AUC WITH MAIN DEMOGRAPHIC VARIABLES ####

  hist(df_trial$ageyrs) ## pick 10 - 25 years
  hist(df_trial$t1_piq) ### IQ 80-135

  labels <- c("201" = "standard", "202" = "oddball pitch", "203" = "oddball length", "204" = "oddball both")
  #age
  g1<-ggplot(df_trial[df_trial$ageyrs>=10 & df_trial$ageyrs<=25,],aes(x=ageyrs,y=scale(rpd_auc),group=interaction(EventData,t1_diagnosis),color=t1_diagnosis,fill=t1_diagnosis))+
    geom_smooth(method='lm')+facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    scale_fill_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='age (years)',y='pupillary response (z)')+
    coord_cartesian(ylim = c(-0.05, 0.15))+
    theme_bw()

  #IQ
  g2<-ggplot(df_trial[df_trial$t1_piq > 80 & df_trial$t1_piq < 135,],aes(x=t1_piq,y=scale(rpd_auc),group=interaction(EventData,t1_diagnosis),color=t1_diagnosis,fill=t1_diagnosis))+
    geom_smooth(method='lm')+facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    scale_fill_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='perceptual IQ',y='pupillary response (z)')+
    coord_cartesian(ylim = c(-0.05, 0.15))+
    theme_bw()


    #summary functions
    stat_sum_df <- function(fun, geom="crossbar", ...) {
      stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
    }

    stat_sum_single <- function(fun, geom="point", ...) {
      stat_summary(fun=fun, geom=geom, size = 3, ...)
    }

  group_labels<-c("female","female","male","male")
  g3<-ggplot(df_trial,aes(x=interaction(t1_diagnosis,sex),y=scale(rpd_auc),group=interaction(EventData,t1_diagnosis,sex),color=t1_diagnosis))+
    stat_sum_df("mean_cl_normal", geom = "errorbar")+
    stat_sum_single("mean")+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_x_discrete(labels=group_labels)+
    labs(x='biological sex',y='pupillary response (z)')+
    coord_cartesian(ylim = c(-0.05, 0.15))+
    theme_bw()

  grid.arrange(g1,g2,g3,nrow=3)


  ###SCALING of baseline pupil size WITH MAIN DEMOGRAPHIC VARIABLES ####

  hist(df_trial$t1_ageyrs) ## pick 10 - 25 years
  hist(df_trial$t1_piq) ### IQ 80-135

  labels <- c("201" = "standard", "202" = "oddball pitch", "203" = "oddball length", "204" = "oddball both")
  #age
  g1<-ggplot(df_trial[df_trial$ageyrs>=10 & df_trial$ageyrs<=25,],aes(x=ageyrs,y=pd_baseline,group=interaction(EventData,t1_diagnosis),color=t1_diagnosis,fill=t1_diagnosis))+
    geom_smooth(method='lm')+facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    scale_fill_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='age (years)',y='baseline pupil size (mm)')+
    theme_bw()

  #IQ
  g2<-ggplot(df_trial[df_trial$t1_piq > 80 & df_trial$t1_piq < 135,],aes(x=t1_piq,y=pd_baseline,group=interaction(EventData,t1_diagnosis),color=t1_diagnosis,fill=t1_diagnosis))+
    geom_smooth(method='lm')+facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    scale_fill_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='perceptual IQ',y='baseline pupil size (mm)')+
    theme_bw()


  #summary functions
  stat_sum_df <- function(fun, geom="crossbar", ...) {
    stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
  }

  stat_sum_single <- function(fun, geom="point", ...) {
    stat_summary(fun=fun, geom=geom, size = 3, ...)
  }

  group_labels<-c("female","female","male","male")
  g3<-ggplot(df_trial,aes(x=interaction(t1_diagnosis,sex),y=pd_baseline,group=interaction(EventData,t1_diagnosis,sex),color=t1_diagnosis))+
    stat_sum_df("mean_cl_normal", geom = "errorbar")+
    stat_sum_single("mean")+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_x_discrete(labels=group_labels)+
    labs(x='biological sex',y='baseline pupil size (mm)')+
    theme_bw()

  grid.arrange(g1,g2,g3,nrow=3)
  ####--> higher IQ and lower age in autistic females?

      with(df_dem,by(t1_piq,interaction(t1_sex,t1_diagnosis),summary))
      with(df_dem,by(t1_ageyrs,interaction(t1_sex,t1_diagnosis),summary))

  ###SCALING of MMN WITH MAIN DEMOGRAPHIC VARIABLES ####

  hist(df_trial$t1_ageyrs) ## pick 10 - 25 years
  hist(df_trial$t1_piq) ### IQ 80-135

  labels <- c("201" = "standard", "202" = "oddball pitch", "203" = "oddball length", "204" = "oddball both")
  #age
  g1<-ggplot(df_trial[df_trial$ageyrs>=10 & df_trial$ageyrs<=25,],aes(x=ageyrs,y=mmn,group=interaction(EventData,t1_diagnosis),color=t1_diagnosis,fill=t1_diagnosis))+
    geom_smooth(method='lm')+facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    scale_fill_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='age (years)',y='mismatch negativity (uV)')+
    theme_bw()

  #IQ
  g2<-ggplot(df_trial[df_trial$t1_piq > 80 & df_trial$t1_piq < 135,],aes(x=t1_piq,y=mmn,group=interaction(EventData,t1_diagnosis),color=t1_diagnosis,fill=t1_diagnosis))+
    geom_smooth(method='lm')+facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    scale_fill_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    labs(x='perceptual IQ',y='mismatch negativity (uV)')+
    theme_bw()


  #summary functions
  stat_sum_df <- function(fun, geom="crossbar", ...) {
    stat_summary(fun.data=fun, geom=geom, width=0.2, ...)
  }

  stat_sum_single <- function(fun, geom="point", ...) {
    stat_summary(fun=fun, geom=geom, size = 3, ...)
  }

  group_labels<-c("female","female","male","male")
  g3<-ggplot(df_trial,aes(x=interaction(t1_diagnosis,sex),y=mmn,group=interaction(EventData,t1_diagnosis,sex),color=t1_diagnosis))+
    stat_sum_df("mean_cl_normal", geom = "errorbar")+
    stat_sum_single("mean")+
    scale_color_manual(name = "group", labels = c("autistic", "non-autistic"),values = brewer.pal(3, "Dark2")[1:2])+
    facet_grid(~EventData,labeller=labeller(EventData = labels))+
    scale_x_discrete(labels=group_labels)+
    labs(x='biological sex',y='mismatch negativity (uV)')+
    theme_bw()

  grid.arrange(g1,g2,g3,nrow=3)
  ####--> higher IQ and lower age in autistic females?

  with(df_dem,by(t1_piq,interaction(t1_sex,t1_diagnosis),summary))
  with(df_dem,by(t1_ageyrs,interaction(t1_sex,t1_diagnosis),summary))



### MODERATION analysis of BPS

  #moderation of BPS on SEPR interaction with group
  lmm_REML_BPS<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*EventCounter+scale(pd_baseline)+
                   (1|subjects)+(1|wave),data=df_trial)

  anova(lmm_REML_BPS)


  #moderation of BPS on MMN interaction with group
  lmm_REML_BPS<-lmer(scale(mmn)~EventData*t1_diagnosis*EventCounter+scale(pd_baseline)+
                       (1|subjects)+(1|wave),data=df_trial)

  anova(lmm_REML_BPS)



# -- full VISUALIZATION (df - rpd change)####
df_dem_merge<-df_dem[,c('subjects','t1_group','t1_diagnosis','t1_sex','t1_ageyrs','t2_ageyrs','t1_fsiq','t1_piq','t1_viq')]
df_dem_merge$subjects<-as.character(df_dem_merge$subjects)
#df$subjects<-substr(df$id,1,12)
df<-merge(df,df_dem_merge,by='subjects') ###--> takes a while

    test<-df[sample(1:nrow(df),nrow(df)/100),]
    test<-df

    ggplot(test[test$time_event<0.60,],aes(x=time_event,y=rpd,group=interaction(t1_diagnosis,EventData),color=EventData,linetype=t1_diagnosis))+geom_smooth()+theme_bw()

    ggplot(test[test$EventCounter<1400 & test$EventData==201,],aes(x=EventCounter,y=pd_baseline,group=t1_diagnosis,color=t1_diagnosis))+geom_smooth()+geom_smooth(linetype=2,method='lm')
    ##--> characteristics of the movie?


    test$trial_position<-ifelse(test$EventCounter<=280,'first quintile',
                                ifelse(test$EventCounter>1120,'fifth quintile',
                                       ifelse(test$EventCounter>560 & test$EventCounter<=840,'third quintile',
                                              ifelse(test$EventCounter>280 & test$EventCounter<=560,'second quintile',
                                                     ifelse(test$EventCounter>840 & test$EventCounter<=1120,'fourth quintile',NA)))))
    test$trial_position<-factor(test$trial_position,levels=c('first quintile','second quintile','third quintile','fourth quintile','fifth quintile'))

    ggplot(test[test$time_event<0.70 & !is.na(test$trial_position),],aes(x=time_event,y=rpd,group=interaction(t1_diagnosis,EventData),color=EventData,linetype=t1_diagnosis))+
      facet_grid(EventData ~ trial_position)+
      geom_smooth(formula = y ~ x + poly(x,4))+theme_bw()
    ###--> takes some time


    ###TASK
    test<-df[sample(1:nrow(df),nrow(df)/100),]
    ggplot(df[df$time_event<0.6,],aes(x=time_event,y=rpd,group=EventData,color=EventData,fill=EventData))+
      geom_smooth(alpha=0.4)+
      labs(x='trial duration (s)',y='change in pupil size (mm)')+
      scale_fill_discrete(labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))+
      scale_color_discrete(labels=c("201" = "standard", "202" = "pitch oddball", "203" = "length oddball ", "204" = "pitch & length oddball"))+
      theme_bw()



    #with(df_trial,by(center_dev,sampling_rate,summary)) #does not differ between sampling rate and thus eye-trackers
    #plot gaze behavior
    require(hexbin)
    require(gridExtra)
    ggplot(df_trial[df_trial$sampling_rate=='300Hz',],aes(x=gazepos.x,y=gazepos.y))+
      stat_density_2d(aes(fill=after_stat(density)), n=100, geom = "raster", contour = FALSE)+
      scale_fill_gradientn(colours=rev(rainbow(3)))+coord_fixed(ratio = 9/16)

    ggplot(df_trial[df_trial$sampling_rate=='120Hz',],aes(x=gazepos.x,y=gazepos.y))+
      stat_density_2d(aes(fill=after_stat(density)), n=100, geom = "raster", contour = FALSE)+
      scale_fill_gradientn(colours=rev(rainbow(3)))+coord_fixed(ratio = 3/4)


hist(test$time_event[test$time_event<1])



### TESTING Bayesian models ####

    # #--> Bayesian modelling using STAN
    # #install.packages('brms')
    # require(brms)
    # help('brms')
    #
    # ###requires Rtools - however Rstan ist not compatible with R4.2 right now
    # #You will need to install the preview of rstan 2.26 using:
    # # install.packages('Rtools')
    # # install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    # # install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))


### linear mixed model - frequentist approach

#response as a function of basleine

lmm<-lmer(scale(rpd_auc)~EventData*t1_diagnosis*EventCounter+
            # scale(ageyrs)+scale(t1_piq)+sex+
            # as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)+
            (1|subjects),data=df_trial)

anova(lmm)
fixef(lmm)['scale(pd_baseline)']
emmeans(lmm,~EventData)
###-> ###main effect beta=-0.44


#distributions
hist(scale(df_trial$pd_baseline))
psych::describe(scale(df_trial$pd_baseline))


table(df_trial$EventData)

psych::describe(scale(df_trial$rpd_response))

### - use BRMS ####

require(brms)

names(df_trial)

df_trial$rpd_auc_z<-scale(df_trial$rpd_auc)
df_trial$pd_baseline_z<-scale(df_trial$pd_baseline)

def_formula <- rpd_auc_z ~ pd_baseline_z + (1|subjects)

get_prior(def_formula, data  = df_trial)

pr <- c(set_prior("normal(0,1)", class = "b", coef = "pd_baseline_z"))

# verify that the priors indeed found their way into Stan's model code
make_stancode(def_formula,  data = df_trial, prior = pr)

bayesian_mixed <- brm(
  formula = def_formula,
  data  = df_trial,
  prior = pr,
  #chains = 4,
  #cores = 24,
  #threads = 8,
  #iter = 100
)


#more chains needed
summary(bayesian_mixed)
bayes_R2(bayesian_mixed) #R² for bayesian models
bayes_factor(bayesian_mixed,bayesian_mixed2) #used to compare models


### - USE greta ####

#install.packages('greta')
#install_greta_deps() ###install tensorflow on a specific python version - thus a miniconda environment is installed - so that this version of python is only seen by Greta
library(greta)
require(bayesplot) #plotting
require(DiagrammeR) #plotting

    #--- TESTING - predict response from baseline --> works fine ####
    ###vignette
    # data - create greta arrays
    x <- as_data(scale(df_trial$pd_baseline))
    y <- as_data(scale(df_trial$rpd_auc))

    # variables and priors
    int <- normal(0, 1)
    coef <- normal(0, 3)
    sd <- student(3, 0, 1, truncation = c(0, Inf))

    # operations
    mean <- int + coef * x

    # likelihood - distribution over data
    distribution(y) <- normal(mean, sd)

    # defining the model
    m <- model(int, coef, sd)

    # plotting
    require(DiagrammeR)
    plot(m)


    # sampling
    draws <- mcmc(m, n_samples = 1000)
    summary(draws)

#--- approach our main model ####
#group <- as_data(ifelse(df_trial$t1_diagnosis=='ASD',T,F))

# A. translate data - make greta arrays
y <- as_data(scale(df_trial$rpd_auc))
group <- as_data(ifelse(df_trial$t1_diagnosis=='ASD',T,F))
events <- as_data(model.matrix(~ EventData - 1, df_trial)) #alternative to dummy coding?
trials <- as_data(scale(df_trial$EventCounter))
sequence_pos <- as_data(scale(df_trial$sequence_position))

subject_id <- as.integer(factor(df_trial$subjects)) #random intercept (in LMM)

interaction_variable<-with(df_trial,interaction(t1_diagnosis,EventData))
group_by_events <- as_data(model.matrix(~interaction_variable-1, df_trial)) #alternative to dummy coding - remove intercept (-1)
#table(group_by_events[,1],group_by_events[,9])

#define priors
event_coefs<- normal(0, 1, dim=ncol(events))
group_coef<- normal(0, 1)

group_event_coefs <- normal(0,1, dim=ncol(group_by_events))
#names(group_event_coefs)<-c('ASD201','TD201','ASD202','TD202','ASD203','TD203','ASD204','TD204')

sd <- cauchy(0, 3, truncation = c(0, Inf))
a_subject <- normal(0, 1, dim = max(subject_id))

    #relevant for random slope model
    a_subject_intercept <- normal(0, 1, dim = max(subject_id))
    a_subject_slope <- normal(0, 1, dim = max(subject_id))

#OPERATIONS - MODEL
#model: group * events
#mean <- group_by_events %*% group_event_coefs + a_subject[subject_id]
#model: trials * group * events - random intervept model
#mean <- trials * group_by_events %*% group_event_coefs + a_subject[subject_id]
mean <- trials * group_by_events %*% group_event_coefs + a_subject[subject_id]

#model: trials + group * events - random slope model
#mean <- trials * group_by_events %*% group_event_coefs + a_subject_intercept[subject_id] + trials * a_subject_slope[subject_id]

#mean <- events %*% event_coefs + group * group_coef + a_subject[subject_id]

# likelihood - distribution over data
distribution(y) <- normal(mean, sd)

# defining the model
#m <- model(int, coef_201, coef_202, coef_203, coef_204, sd)
m <- model(group_event_coefs)

# plotting
plot(m)

# sampling
#draws <- mcmc(m, warmup = 100,n_samples = 100) #testing
start <-Sys.time()
draws <- mcmc(m, n_samples = 1000, warmup = 300,) ###--> test
Sys.time() - start

# ##save model
# bm_random_intercept<-draws
# save(bm_random_intercept,file='C:/Users/nico/Desktop/mmn_leap_bayesian_model_23022023')

#analyse
summary(draws)
plot(draws)
mcmc_trace(draws) #see whether chains converge and are stationary
mcmc_dens(draws) # + stable estimates

##results
mcmc_intervals(draws, prob= .66, prob_outer = .95)+scale_y_discrete(labels=c('ASD201','TD201','ASD202','TD202','ASD203','TD203','ASD204','TD204'))+coord_flip()+theme_bw()
mcmc_areas(draws, prob= .66, prob_outer = .95)+scale_y_discrete(labels=c('ASD201','TD201','ASD202','TD202','ASD203','TD203','ASD204','TD204'))+coord_flip()+theme_bw()

param_draws_df <- data.frame(draws[[1]])
names(param_draws_df)<-c('ASD201','TD201','ASD202','TD202','ASD203','TD203','ASD204','TD204')

    most_likely_coefs<-unlist(opt(m)['par'])
    predicted_y<- model.matrix(~interaction_variable-1, df_trial) %*% most_likely_coefs
    predicted_df<-data.frame(predicted_y,model.matrix(~interaction_variable-1, df_trial))
    predicted_df<-reshape2::melt(predicted_df, id = 'predicted_y')
    ggplot(predicted_df,aes(x=predicted_y))+geom_density()+facet_wrap(~variable+value)


#retrieve posterior predictve values
test <- calculate(mean, values = draws)


#posterior predictive checking
coda::gelman.diag(draws) #Rhat <= 1.01 is cutoff for reliable estimator
coda::effectiveSize(draws) #ESS

save(draws,draws2,file='C:/Users/nico/Desktop/mmn_leap_bayesian_models_24022023')
#draws - rpd_auc ~ trials * group * events
#draws2 - rpd_auc ~ sequence_pos * group * events


### - BAYESIAN TEST PLAYGROUND ####

###testing bayesian sample models ####
data(attitude)
design <- as.matrix(attitude[, 2:7])
int <- normal(0, 10)
sd <- cauchy(0, 3, truncation = c(0, Inf))

tau <- exponential(0.5, dim = ncol(design))
coefs <- normal(0, tau)

mu <- int + design %*% coefs

distribution(attitude$rating) <- normal(mu, sd)

#defining the model
m <- model(int, coef, sd)

# plotting
plot(m)

# sampling
draws <- mcmc(m, n_samples = 1000)

summary(draws)


# linear model parameters
int <- normal(0, 10)
coef <- normal(0, 10)
sd <- cauchy(0, 3, truncation = c(0, Inf))

# hierarchical model for species effect; use the first species as the baseline
# like in lm()
species_sd <- lognormal(0, 1)
species_offset <- normal(0, species_sd, dim = 2)
species_effect <- rbind(0, species_offset)
species_id <- as.numeric(iris$Species)

# model
mu <- int + coef * iris$Sepal.Width + species_effect[species_id]
distribution(iris$Sepal.Length) <- normal(mu, sd)

# defining the model
m <- model(int, coef, sd, species_effect)

# plotting
plot(m)

# sampling
draws <- mcmc(m, n_samples = 1000)

summary(draws)


draws <- mcmc(m, n_samples = number_of_iterations, warmup = number_of_warmups) ###--> test



