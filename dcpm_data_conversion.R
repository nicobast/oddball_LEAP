#convert data for DCPM

### taken from dataAnalysis.m of DCPM package

# Inputs:
# %
# %              input struct "in" with fields containing cell array with one
# %              cell per run of data. Field names:
# %
# %              - yPos: horizontal gaze position in dva (column vector)
# %              - yPos: vertical gaze position in dva (column vector)
# %                      (the coodinates (0,0) should be at fixation)
# %              - pupilArea: pupil area (column vector)
# %              - startInds: trial start indexes in samples (n x 2 matrix
#                                                             %                with trial start and end times, n is number of trials)
# %              - sampleRate: sampling rate of eye tracker (Hz)
# %              - trialTypes: trial types (e.g., easy / hard, or corr\error)
# %                integers for each trial, e.g. 1,2,3,4,5 (1 x n vector)
# %              - predictionWindow: time window beyond trial onset to make a
# %                prediction within (e.g. 4 sec in a jittered ISI expt)

require(R.matlab)
require(ggplot2)


    test_data<-df[df$id=='628118917744_wave1',] #sampleData2.mat --> works great
    test_data<-df[df$id=="976537183752_wave1",] #test_dcpm_input2.mat
    test_data<-df[df$id=="133731063641_wave1",] #test_dcpm_input3.mat
    test_data<-df[df$id=="251317189021_wave1",]
    test_data<-df[df$id=="761647460630_wave2",] #test_dcpm_input5.mat --> works great


    unique(test_data$frequency_rate)

    ggplot(test_data[test_data$EventData=='201',],aes(x=time_event,y=pd))+geom_smooth()
    ggplot(test_data[test_data$time_event<0.65 & test_data$EventData=='201',],aes(x=time_event,y=pd))+geom_smooth()
    ggplot(test_data[test_data$time_event<0.6,],aes(x=time_event,y=pd,group=EventData,color=EventData))+geom_smooth()
    #628118917744_wave1 --> with clean response

sample_ids<-sample(unique(df$id),10)
sample_ids<-unique(df$id)

for(participant in seq(sample_ids)){

#select data of a participant
test_data<-df[df$id==sample_ids[participant],] #test_dcpm_input5.mat --> works great
##fix to standard trials
test_data<-test_data[test_data$EventData=='201',]


## convert gaze coord from relative space to degress in visual angle
###formula: visual angle = 2 x atan (0.5* gaze_coord/screen_distance*0.1)
screen_width<-345 #mm fixed width of presentation screen EU  AIMS LEAP
screen_height<-259 #mm fixed height of presentation screen EU  AIMS LEAP
degrees_by_radian<-360/(2*pi) #fixed conversion facor
screen_dist<-test_data$screen_dist
x_norm <- test_data$gazepos.x-0.5
y_norm <- test_data$gazepos.y-0.5
x_dva<-degrees_by_radian*2*atan(x_norm*screen_width/(2*screen_dist))
y_dva<-degrees_by_radian*2*atan(y_norm*screen_height/(2*screen_dist))
# a = 2 * arctan(size / (2 * distance))

###convert pupil size to Arbitary units as recorded by Eye-Link
# % convert pupil size in mm to AREA
# % see https://doi.org/10.3758/s13428-015-0588-x, Experiment 1 end
# % AU = pupil_size / (rescaling_in_radians_of_AU * distance_to_screen )
radians_per_AU<-1.7*10^(-4) #see https://doi.org/10.3758/s13428-015-0588-x, Experiment 1 end
pd_AU <- with(test_data,pd/(radians_per_AU*screen_dist))

###recalculations for correct data format
trial_start<-which(test_data$ts_event==1)
trial_end<-cumsum(with(test_data,by(ts_event,EventCounter,max)))
trial_number<-as.numeric(with(test_data,by(EventCounter,EventCounter,head,n=1)))

event_types<-as.character(with(test_data,by(EventData,EventCounter,head,n=1)))
event_types<-sapply(event_types,function(x){switch(x,
                                                   "201" = 1,
                                                   "202" = 2,
                                                   "203" = 2,
                                                   "204" = 2)})

# event_types<-sapply(event_types,function(x){switch(x,
#                                                    "201" = 1,
#                                                    "202" = 1,
#                                                    "203" = 1,
#                                                    "204" = 1)})


# #put it together
matlab_input<-list()
matlab_input[[1]]<-x_dva #0,0 should be center
matlab_input[[2]]<-y_dva  #0,0 should be center
matlab_input[[3]]<-pd_AU
matlab_input[[4]]<-as.matrix(cbind(trial_start,trial_end))
matlab_input[[5]]<-with(test_data,round(1/frequency_rate)[1])
matlab_input[[6]]<-event_types


# # ###for different runs
#
# split_point<-mean(test_data$EventCounter)
# split_runs<-with(test_data,ifelse(EventCounter<split_point,1,ifelse(EventCounter>split_point,2,NA)))
#
# matlab_input<-list()
# matlab_input[[1]]<-split(x_dva,split_runs) #0,0 should be center
# matlab_input[[2]]<-split(y_dva,split_runs)  #0,0 should be center
# matlab_input[[3]]<-split(pd_AU,split_runs)
# matlab_input[[4]]<-split(data.frame(trial_start,trial_end),trial_number>split_point)
# matlab_input[[5]]<-list(with(test_data,round(1/frequency_rate)[1]),with(test_data,round(1/frequency_rate)[1]))
# matlab_input[[6]]<-split(event_types,trial_number>split_point)
# ##reset startinds of run 2 (relative to absolute)
# matlab_input[[4]][[2]]<-matlab_input[[4]][[2]]-matlab_input[[4]][[2]][1,1]+1

names(matlab_input)<-c('xPos','yPos','pupilArea','startInds','sampleRate','trialTypes')


#writeMat(paste0('data/pcdm_input/',sample_ids[participant],'.mat'),input=matlab_input) #input - is the name how it appears in matlab workspace
writeMat(paste0('C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/input/',sample_ids[participant],'.mat'),input=matlab_input) #input - is the name how it appears in matlab workspace
#writeMat('data/pcdm_input/test_data5.mat',input=matlab_input) #input - is the name how it appears in matlab workspace

print(paste0('saved: ',sample_ids[participant]))

}


# #read matlab output form DCPM
# mat_output<-readMat("C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/sample_ouput.mat")
# ?readMat

##output is saved as cell types
mat_output<-readMat("C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/output/estimates.mat")
test_mat<-mat_output[[1]][13][[1]][[1]] #extract relevant data of single participant
names(test_mat)<-dimnames(test_mat)[[1]] #assign names of dimensions
test_save<-unlist(test_mat['gain'])

