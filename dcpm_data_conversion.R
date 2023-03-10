#convert data for DCPM

### taken from dataAnalysis.m of DCPM package

# Inputs:
#   %
# %              input struct "in" with fields containing cell array with one
# %              cell per run of data. Field names:
#   %
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

unique(df$id)
test_data<-df[df$id=='323195490394_wave2',]

matlab_input<-list()
matlab_input[[1]]<-list(test_data$gazepos.x-0.5)
matlab_input[[2]]<-test_data$gazepos.y-0.5
matlab_input[[3]]<-test_data$pd

trial_start<-which(test_data$ts_event==1)
trial_end<-cumsum(with(test_data,by(ts_event,EventCounter,max)))

matlab_input[[4]]<-as.matrix(cbind(trial_start,trial_end))
matlab_input[[5]]<-with(test_data,round(1/frequency_rate)[1])
matlab_input[[6]]<-as.character(with(test_data,by(EventData,EventCounter,head,n=1)))
matlab_input[[7]]<-0.55 #prediction window - length of average trial

names(matlab_input)<-c('xPos','yPos','pupilArea','startInds','sampleRate','trialTypes','predictionWindow')

writeMat('data/test_dcpm_input.mat',input=matlab_input) #input - is the name how it appears in matlab workspace



table(df$frequency_rate)

#read matlab output form DCPM
mat_output<-readMat("C:/Users/nico/PowerFolders/project_matlab_pcdm/PCDM/sample_ouput.mat")
?readMat
