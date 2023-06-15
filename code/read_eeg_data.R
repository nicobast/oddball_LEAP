
### INFO ####

# EEG data saved with matlab package EEGLAB - but separated as .set and .fdt file
# .fdt files are bianry and need to be provide with the correct data shape:
# https://www.mattcraddock.com/blog/2017/11/17/loading-eeglab-set-files-in-r-part-2/

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



### paths
data_path<-"Z:/nico/backup_data/data_AIMS/mismatch_negativity/EEG/05_segmented_baselinecorrected"
file_list<-list.files(data_path,recursive=T)

##files
file_list_set<-file_list[grepl(".set",file_list)] #.set files only contain the data structure
file_list_fdt<-file_list[grepl(".fdt",file_list)] #.fdt contains actual data - but as binary format

###analyze data structure

number_of_conditions<-4
length(file_list_set)/number_of_conditions
###--> 462 data sets


###READ DATA ####

#read meta data
one_set<-readMat(paste(data_path,file_list_set[30],sep='/'))

#meta data structure
one_set$EEG

#### data structure required to reconstruct correct matrix shape from binary file
var_names<-dimnames(one_set$EEG)[[1]] #get variable names of meta data
n_chans <- one_set$EEG[[which(var_names == 'nbchan')]] #number of channels
n_trials <- one_set$EEG[[which(var_names == 'trials')]] #number of retained trials
times <- one_set$EEG[[which(var_names == 'times')]] #time within trial (-100-500ms)

#name of data file in meta data
one_data_file_path<-unlist(one_set[[1]][43])


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


#other info
ID<-substr(one_data_file_path,1,12)
condition<-substr(one_data_file_path,14,16)
trial_info<-as.numeric(one_set$EEG[[which(var_names == 'trialskeep')]])

#TODO:
#wave?
#--> batch across participants

#### trial information ####

      #read meta data
      one_set<-readMat(paste(data_path,file_list_set[30],sep='/'))
      var_names<-dimnames(one_set$EEG)[[1]] #get variable names of meta data

###-->trial information in meta data? - trials keep - not succesful ####
var_names
test<-one_set$EEG[[which(var_names == 'trialskeep')]]
as.numeric(test)
cumsum(as.numeric(test))
plot(density(as.numeric(test)))
mean(test)
hist(test)
sum(test)

fun_trials_retained<-function(x){
  one_set<-readMat(paste(data_path,x,sep='/'))
  var_names<-dimnames(one_set$EEG)[[1]]
  trialskeep<-as.numeric(one_set$EEG[[which(var_names == 'trialskeep')]])
  return(trialskeep)
}

list_trialsretained<-pblapply(file_list_set[1:30],fun_trials_retained)

  sapply(list_trialsretained,function(x){cumsum(x)})
  hist(unlist(list_trialsretained),10)

  cumsum(unlist(list_trialsretained))
  unlist(list_trialsretained)

###--> might be triggers since last trial? --> not successful ####

#original event information

fun_trigger_seq<-function(x){

  #TESTING
  one_set<-readMat(paste(data_path,x,sep='/'))

  events_all<-one_set$EEG[[which(var_names == 'urevent')]] #event trigger in raw data
  events_all<-data.frame(matrix(unlist(events_all),ncol=4,byrow=T)) #bring matrix into right format
  trigger_sequence<-events_all[,1] #extract trigger

          # which(trigger_sequence=='T204')
          # which(trigger_sequence=='204')
          # which(trigger_sequence=='S204')

  #identify first oddball trigger
  first_oddball_trigger<-min(c(grep('201',trigger_sequence),
                               grep('202',trigger_sequence),
                               grep('203',trigger_sequence),
                               grep('204',trigger_sequence)))

  #identify last oddball trigger
  last_oddball_trigger<-max(c(grep('201',trigger_sequence),
                              grep('202',trigger_sequence),
                              grep('203',trigger_sequence),
                              grep('204',trigger_sequence)))

  #select oddball trigger sequence
  oddball_trigger_sequence<-trigger_sequence[first_oddball_trigger:last_oddball_trigger]

  return(oddball_trigger_sequence)

}

list_triggersequences<-pblapply(file_list_set[1:30],fun_trigger_seq)


    #checking
    sapply(list_triggersequences,table)
    list_clean_triggersequences<-list_triggersequences[sapply(list_triggersequences,function(x){length(table(x))})==4]

    list_trialsretained
    list_clean_trialsretained<-list_trialsretained[sapply(list_triggersequences,function(x){length(table(x))})==4]

    sapply(list_clean_triggersequences,table,simplify=F)

    hist(unlist(list_clean_trialsretained))
    hist(unlist(list_trialsretained))

    sapply(list_clean_trialsretained,cumsum)


    ##crossreferencing
    a_dataset<-1

    table(unlist(list_clean_triggersequences[a_dataset]))
    length(unlist(list_clean_trialsretained[a_dataset]))
    list_clean_trialsretained[a_dataset] #trial information?
    grep('204',unlist(list_clean_triggersequences[a_dataset])) #which are 204
    list_clean_triggersequences[a_dataset] #trigger sequence

#trial information in epoch data --> SUCCESFUL ####

    #read meta data
    one_set<-readMat(paste(data_path,file_list_set[30],sep='/'))
    var_names<-dimnames(one_set$EEG)[[1]] #get variable names of meta data

#information on retained triggers of trials is in there
epochs_retained<-one_set$EEG[[which(var_names == 'epoch')]]
epochs_retained<-data.frame(matrix(unlist(epochs_retained),ncol=6,byrow=T)) #bring matrix into right format

#information on ALL triggers of trials is in there
events_all<-one_set$EEG[[which(var_names == 'urevent')]] #event trigger in raw data
events_all<-data.frame(matrix(unlist(events_all),ncol=4,byrow=T)) #bring matrix into right format
trigger_sequence<-events_all[,1] #extract trigger

#compare
retained_triggers<-as.numeric(epochs_retained[,6]) ## information which trigger are retained
#trigger_positions_oftype<-grep('204',trigger_sequence) ##information on all triggers of that type

relevant_trigger_sequence<-sort(c(grep('201',trigger_sequence),
                                  grep('202',trigger_sequence),
                                  grep('203',trigger_sequence),
                                  grep('204',trigger_sequence)))

trial_counter<-seq_along(relevant_trigger_sequence)

retained_trigger_of_type<-relevant_trigger_sequence %in% retained_triggers

retained_trials<-trial_counter[retained_trigger_of_type]


##test
length_of_epoch<-table(one_data_set$epoch)[1]



one_data_set$trial_counter<-rep(retained_trials,each=length_of_epoch)

###--> batch for all participants ####

###testing

colnames(one_data_set)

ggplot(one_data_set,aes(x=times,y=T8))+geom_smooth()

