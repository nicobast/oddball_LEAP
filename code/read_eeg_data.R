
### INFO ####

# EEG data saved with matlab package EEGLAB - but separated as .set and .fdt file
# .fdt files are bianry and need to be provide with the correct data shape:
# https://www.mattcraddock.com/blog/2017/11/17/loading-eeglab-set-files-in-r-part-2/

##object chanlocs could be applied to retrieve further info on channels


#INFO ON EEGLAB file structure form eeg_checkset.m:

# % The structure of an EEG dataset under EEGLAB (as of v5.03):
#   %
# % Basic dataset information:
#   %   EEG.setname      - descriptive name|title for the dataset
# %   EEG.filename     - filename of the dataset file on disk
# %   EEG.filepath     - filepath (directory/folder) of the dataset file(s)
# %   EEG.trials       - number of epochs (or trials) in the dataset.
# %                      If data are continuous, this number is 1.
# %   EEG.pnts         - number of time points (or data frames) per trial (epoch).
# %                      If data are continuous (trials=1), the total number
# %                      of time points (frames) in the dataset
# %   EEG.nbchan       - number of channels
# %   EEG.srate        - data sampling rate (in Hz)
# %   EEG.xmin         - epoch start latency|time (in sec. relative to the
#                                                  %                      time-locking event at time 0)
# %   EEG.xmax         - epoch end latency|time (in seconds)
# %   EEG.times        - vector of latencies|times in miliseconds (one per time point)
# %   EEG.ref          - ['common'|'averef'|integer] reference channel type or number
# %   EEG.history      - cell array of ascii pop-window commands that created
# %                      or modified the dataset
# %   EEG.comments     - comments about the nature of the dataset (edit this via
#                                                                  %                      menu selection Edit > About this dataset)
# %   EEG.etc          - miscellaneous (technical or temporary) dataset information
# %   EEG.saved        - ['yes'|'no'] 'no' flags need to save dataset changes before exit
# %
# % The data:
#   %   EEG.data         - two-dimensional continuous data array (chans, frames)
# %                      ELSE, three-dim. epoched data array (chans, frames, epochs)
# %
# % The channel locations sub-structures:
#   %   EEG.chanlocs     - structure array containing names and locations
# %                      of the channels on the scalp
# %   EEG.urchanlocs   - original (ur) dataset chanlocs structure containing
# %                      all channels originally collected with these data
# %                      (before channel rejection)
# %   EEG.chaninfo     - structure containing additional channel info
# %   EEG.ref          - type of channel reference ('common'|'averef'|+/-int]
# %   EEG.splinefile   - location of the spline file used by headplot() to plot
# %                      data scalp maps in 3-D
# %
# % The event and epoch sub-structures:
#   %   EEG.event        - event structure containing times and nature of experimental
# %                      events recorded as occurring at data time points
# %   EEG.urevent      - original (ur) event structure containing all experimental
# %                      events recorded as occurring at the original data time points
# %                      (before data rejection)
# %   EEG.epoch        - epoch event information and epoch-associated data structure array (one per epoch)
# %   EEG.eventdescription - cell array of strings describing event fields.
# %   EEG.epochdescription - cell array of strings describing epoch fields.
# %   --> See the http://sccn.ucsd.edu/eeglab/maintut/eeglabscript.html for details
# %
# % ICA (or other linear) data components:
#   %   EEG.icasphere   - sphering array returned by linear (ICA) decomposition
# %   EEG.icaweights  - unmixing weights array returned by linear (ICA) decomposition
# %   EEG.icawinv     - inverse (ICA) weight matrix. Columns gives the projected
# %                     topographies of the components to the electrodes.
# %   EEG.icaact      - ICA activations matrix (components, frames, epochs)
# %                     Note: [] here means that 'compute_ica' option has bee set
# %                     to 0 under 'File > Memory options' In this case,
# %                     component activations are computed only as needed.
# %   EEG.icasplinefile - location of the spline file used by headplot() to plot
# %                     component scalp maps in 3-D
# %   EEG.chaninfo.icachansind  - indices of channels used in the ICA decomposition
# %   EEG.dipfit      - array of structures containing component map dipole models
# %
# % Variables indicating membership of the dataset in a studyset:
#   %   EEG.subject     - studyset subject code
# %   EEG.group       - studyset group code
# %   EEG.condition   - studyset experimental condition code
# %   EEG.run         - studyset run number
# %   EEG.session     - studyset session number
# %
# % Variables used for manual and semi-automatic data rejection:
#   %   EEG.specdata           - data spectrum for every single trial
# %   EEG.specica            - data spectrum for every single trial
# %   EEG.stats              - statistics used for data rejection
# %       EEG.stats.kurtc    - component kurtosis values
# %       EEG.stats.kurtg    - global kurtosis of components
# %       EEG.stats.kurta    - kurtosis of accepted epochs
# %       EEG.stats.kurtr    - kurtosis of rejected epochs
# %       EEG.stats.kurtd    - kurtosis of spatial distribution
# %   EEG.reject            - statistics used for data rejection
# %       EEG.reject.entropy - entropy of epochs
# %       EEG.reject.entropyc  - entropy of components
# %       EEG.reject.threshold - rejection thresholds
# %       EEG.reject.icareject - epochs rejected by ICA criteria
# %       EEG.reject.gcompreject - rejected ICA components
# %       EEG.reject.sigreject  - epochs rejected by single-channel criteria
# %       EEG.reject.elecreject - epochs rejected by raw data criteria


### SETUP ####

#### packages
require(R.matlab) #read matlab files
require(dplyr) #mutate function
require(pbapply)
require(zoo) #rolling mean

#plyr is implicitly called but not loaded as it overshadows some dplyr functions


### paths ####
data_path<-"Z:/nico/backup_data/data_AIMS/mismatch_negativity/EEG/05_segmented_baselinecorrected"
file_list<-list.files(data_path,recursive=T)

##files
file_list_set<-file_list[grepl(".set",file_list)] #.set files only contain the data structure
file_list_fdt<-file_list[grepl(".fdt",file_list)] #.fdt contains actual data - but as binary format

    # ###exclude participants with missing data --> FIXED BELOW
    # file_list_set<-file_list_set[-grep('196084871273',file_list_set)] #have the trial are missing for 204 binary
    # file_list_fdt<-file_list_fdt[-grep('196084871273',file_list_fdt)]
    # file_list_set<-file_list_set[-grep('321466110510',file_list_set)] #have the trial are missing for 204 binary
    # file_list_fdt<-file_list_fdt[-grep('321466110510',file_list_fdt)]
    # file_list_set<-file_list_set[-grep('322708064836',file_list_set)] #have the trial are missing for 204 binary
    # file_list_fdt<-file_list_fdt[-grep('322708064836',file_list_fdt)]
    # file_list_set<-file_list_set[-grep('495779827165',file_list_set)] #have the trial are missing for 204 binary
    # file_list_fdt<-file_list_fdt[-grep('495779827165',file_list_fdt)]



      ###analyze data structure
      number_of_conditions<-4
      length(file_list_set)/number_of_conditions
      ###--> 462 data sets


###DEFINE FUNCTION ####

fun_read_eeglab_output<-function(files_to_read){

  print(paste0('processing: ',files_to_read))

###READ DATA ####

      #testing
      # testing_set<-188
      # one_set<-readMat(paste(data_path,file_list_set[testing_set],sep='/'))
      # "combined/196084871273_204_CombinedNB.set"
      # length(file_list_fdt)
      # length(unique(file_list_fdt))
      #one_set<-readMat(paste(data_path,"combined/196084871273_204_CombinedNB.set",sep='/'))
      #files_to_read<-file_list_set[1]
#
#       one_set<-readMat(paste(data_path,file_list_set[104],sep='/'))
#       one_data_file_path<-paste0(substr(file_list_set[104],1,nchar(file_list_set[104])-4),'.fdt')
#
#       one_set<-readMat(paste(data_path,failed_files[8],sep='/'))
#       one_data_file_path<-paste0(substr(failed_files[8],1,nchar(failed_files[8])-4),'.fdt')


#read meta data
one_set<-readMat(paste(data_path,files_to_read,sep='/'))

# #name of data file in meta data
# one_data_file_path<-unlist(one_set[[1]][43]) #may has more than one entry
# #one_set<-readMat(paste(data_path,"combined/196084871273_204_CombinedNB.set",sep='/')) # this has a weird one_data_file_path output
# #alternative approach
one_data_file_path<-paste0(substr(files_to_read,1,nchar(files_to_read)-4),'.fdt')

#meta data structure
#one_set$EEG

#### data structure required to reconstruct correct matrix shape from binary file
var_names<-dimnames(one_set$EEG)[[1]] #get variable names of meta data
n_chans <- one_set$EEG[[which(var_names == 'nbchan')]] #number of channels
n_trials <- one_set$EEG[[which(var_names == 'trials')]] #number of retained trials
times <- one_set$EEG[[which(var_names == 'times')]] #time within trial (-100-500ms)

##read that data file

#create file connection
one_data_connection<-file(paste(data_path,
                          file_list_fdt[grep(one_data_file_path,file_list_fdt)],
                          sep='/'),'rb') # rb = binary file


#read data from binary in correct matrix shape
one_data_set<-readBin(one_data_connection, #data connection
                      'double', #number type
                      n = n_chans * n_trials * length(times),
                      size=4, #byte size
                      endian='little')

#close connection
close(one_data_connection)

#format to correct matrix dimensions
dim(one_data_set)<- c(n_chans, length(times) * max(n_trials,1))
times <- rep(times, max(n_trials,1))
one_data_set<-data.frame(cbind(t(one_data_set),times))
### column = one channel, row = time

##get channel information
chanlocs <- one_set$EEG[[which(var_names == "chanlocs")]]

##--> chanlocs comes out as a 3D list with 12 * number_of_chans elements, with each list element itself a list
#rbind(chanlocs[,,1]) #required information is in 3rd dimension
chanlocs <- as.data.frame(t(rbind(chanlocs[, , ])))

#name channels
names(one_data_set)<-c(unlist(chanlocs$labels),'times')

##create epoch information
one_data_set <- one_data_set %>%
  group_by(times) %>%
  mutate(epoch = 1:n()) %>%
  ungroup

#### extract trial information from epoch meta data ####

#GET all triggers that were recorded in raw data
# --for some participants - the strucutre of urevent is different - need to be formatted different to retrieve correct trigger_sequence
# --this different structure can be detected in the meta data: one_set[[1]][43][[1]]
if(length(one_set[[1]][43][[1]])!=1){
events_all<-one_set$EEG[[which(var_names == 'urevent')]] #event trigger in raw data
events_all<-data.frame(matrix(unlist(events_all),ncol=2,byrow=T)) #bring matrix into right format - 2 instead of 4 cols
trigger_sequence<-events_all[,1] #extract trigger
}

if(length(one_set[[1]][43][[1]])==1){
  events_all<-one_set$EEG[[which(var_names == 'urevent')]] #event trigger in raw data
  events_all<-data.frame(matrix(unlist(events_all),ncol=4,byrow=T)) #bring matrix into right format
  trigger_sequence<-events_all[,1] #extract trigger
}

#only extract task relevant triggers (for correct trial counter) - some data sets have additional triggers
relevant_trigger_sequence<-sort(c(grep('201',trigger_sequence),
                                  grep('202',trigger_sequence),
                                  grep('203',trigger_sequence),
                                  grep('204',trigger_sequence)))

#create trial sequence based on relevant triggers
trial_counter<-seq_along(relevant_trigger_sequence)

    #DEBUGGING
    print(length(trial_counter)) ###always be 1400


#information on retained triggers of trials is in EPOCH meta data
epochs_retained<-one_set$EEG[[which(var_names == 'epoch')]]
epochs_retained<-data.frame(matrix(unlist(epochs_retained),ncol=6,byrow=T)) #bring matrix into right format
retained_triggers<-as.numeric(epochs_retained[,6]) #information which triggers are retained for final anaylsis


# which trial are retained from the complete sequence
retained_trigger_of_type<-relevant_trigger_sequence %in% retained_triggers

#return trial counter for these retained trials
retained_trials<-trial_counter[retained_trigger_of_type]

    # #checking
      print(table(trigger_sequence[relevant_trigger_sequence])) #only 201, 202, 203, 204
      print(table(trigger_sequence[retained_triggers])) #only current condition of .fdt file

##add to data set
length_of_epoch<-table(one_data_set$epoch)[1]
one_data_set$trial_counter<-rep(retained_trials,each=length_of_epoch)

#### ADD ID, condition, timepoint variable ####

ID<-substr(one_data_file_path,10,22)
condition<-substr(one_data_file_path,1,8)
timepoint<-as.character(one_set$EEG[[which(var_names == 'euaims')]][[5]])

one_data_set$ID<-ID
one_data_set$condition<-condition
one_data_set$timepoint<-timepoint


### DROP CHANNEL to reduce file size ####

one_data_set<-one_data_set[,c('Fz','F1','F2','Cz','Pz','times','epoch','trial_counter','ID','condition')]


return(one_data_set)

}

###--> LOAD Data of some electrodes for all participants ####

###reread - all eeg data - to correct condition and ID
list_eeg<-pblapply(file_list_set,fun_read_eeglab_output)

        # ### correct ID and condition --> has now also been fixed above
        # names(list_eeg)<-file_list_set
        #
        # condition<-substr(file_list_set,1,8)
        # ID<-substr(file_list_set,10,22)
        #
        # list_eeg<-pbmapply(function(x,y){
        #   x$condition<-y
        #   return(x)
        #   },x=list_eeg,y=condition,SIMPLIFY = F)
        #
        # list_eeg<-pbmapply(function(x,y){
        #   x$ID<-y
        #   return(x)
        #   },x=list_eeg,y=ID,SIMPLIFY = F)


save(list_eeg,file="C:/Users/nico/Desktop/leap_oddball_eeg_list.rdata")
df_eeg<-data.table::rbindlist(list_eeg) ##very fast compared to do.call(rbind)


###--> visualize Fz data

require(ggplot2)
require(gridExtra)
theme_set(theme_bw())


df_eeg$condition_recovered<-as.factor(substr(df_eeg$ID,1,8))
table(df_eeg$condition_recovered)

subsample<-sample(1:nrow(df_eeg),nrow(df_eeg)/100)

ggplot(df_eeg[subsample,],aes(x=times,y=Fz,group=condition_recovered,color=condition_recovered))+geom_smooth()

###---> calculate MMN

fun_estimate_mmn<-function(one_set){

  #testing
  #one_set<-list_eeg[[100]]

  ###restrict to timeframe - 35-350ms
  one_set<-one_set[one_set$times>50 & one_set$times<350,]

  mmn_trial<-with(one_set,by(Fz,trial_counter,function(channel_data){

      channel_data_smooth<-rollmean(channel_data,k=20) #rolling mean of 20ms
      min(channel_data_smooth)

  }))

  trial_counter<-names(mmn_trial)
  mmn<-as.numeric(mmn_trial)
  mmn_per_trial<-data.frame(trial_counter,mmn)
  return(mmn_per_trial)

}

list_mmn<-pblapply(list_eeg,fun_estimate_mmn)


###takes around 2 hours

###TESTING ####

      fun_bench_find_datafile<-function(files_to_read){

        print(paste0('processing: ',files_to_read))

        #read meta data
        one_set<-readMat(paste(data_path,files_to_read,sep='/'))

        # #name of data file in meta data
        one_data_file_path<-length(unlist(one_set[[1]][43])) #may has more than one entry
        # #one_set<-readMat(paste(data_path,"combined/196084871273_204_CombinedNB.set",sep='/')) # this has a weird one_data_file_path output
        # #alternative approach
        # one_data_file_path<-paste0(substr(files_to_read,1,nchar(files_to_read)-4),'.fdt')
        return(cbind(files_to_read,one_data_file_path))
      }

      list_bench<-pblapply(file_list_set,fun_bench_find_datafile)
      #--> takes two hours

      list_bench[[1]][2]
      index_of_files_that_fail<-which(sapply(list_bench,function(x){x[2]})!=1)
      failed_files<-sapply(list_bench[index_of_files_that_fail],function(x){x[1]})
      ###--> these fails fail as meta data field refers to wrong data file

      test_list_of_failed<-pblapply(failed_files,fun_read_eeglab_output)
      ###--> process all failed files with modified  function


sample_set<-sample(1:length(file_list_set),30)
test_list<-pblapply(file_list_set[sample_set],fun_read_eeglab_output)


hist(unlist(sapply(test_list,function(x){x['trial_counter']})),10)
summary(unlist(sapply(test_list,function(x){x['trial_counter']})),10)

#--> even distribution of trials

table(unlist(sapply(test_list,function(x){x['epoch']})))


###plot eeg

list_selected_electrodes<-pblapply(test_list,function(x){x[,c('Fz','F1','F2','Cz','Pz','times','epoch','trial_counter','ID','condition')]})
df_eeg<-plyr::rbind.fill(list_selected_electrodes)




ggplot(df_eeg,aes(x=times,y=Fz,group=condition,color=condition))+geom_smooth()


grid.arrange(ggplot(df_eeg,aes(x=times,y=Fz))+geom_smooth(),
             ggplot(df_eeg,aes(x=times,y=Pz))+geom_smooth(),
             ggplot(df_eeg,aes(x=times,y=Cz))+geom_smooth())



###testing

colnames(one_data_set)

table(one_data_set$timepoint)

ggplot(one_data_set,aes(x=times,y=FC2))+geom_smooth()

