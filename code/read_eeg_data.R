
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
require(dplyr) #mutate function - in function

require(pbapply)
require(zoo) #rolling mean
require(data.table) #rbindlist

#visualization
require(wesanderson) #custom colors
require(ggplot2)
require(gridExtra)
require(RColorBrewer)

#plyr is implicitly called but not loaded as it overshadows some dplyr functions

project_path<-"C:/Users/nico/PowerFolders/project_oddball_LEAP"


##custom themes
#define color palette
theme_set(theme_bw())
custom_contrast_colors <- c(brewer.pal(6, "Blues")[4:6],brewer.pal(6, "Oranges")[4:6])



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
       one_set<-readMat(paste(data_path,file_list_set[1],sep='/'))
       one_data_file_path<-paste0(substr(file_list_set[1],1,nchar(file_list_set[104])-4),'.fdt')
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

#one_data_set<-one_data_set[,c('Fz','F1','F2','Cz','Pz','times','epoch','trial_counter','ID','condition')]


return(one_data_set)

}

###--> LOAD EEG Data of all participants ####

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

#save(list_eeg,file="C:/Users/nico/Desktop/leap_oddball_eeg_list.rdata")
save(list_eeg,file="C:/Users/nico/Desktop/leap_oddball_eeg_list_allchannels2.rdata")

###--> SCALP map - requires data of all EEG channels ####

load(file="C:/Users/nico/Desktop/leap_oddball_eeg_list_allchannels2.rdata")

#for standard coordinates
require(eegkit)
data(eegcoord) #load sample set with standard coordinates
rownames(eegcoord)

#for scalp topographical map
# install.packages("remotes")
# remotes::install_github("craddm/eegUtils")
require(eegUtils)

###load a dataset
#load(file="C:/Users/nico/Desktop/leap_oddball_eeg_list_allchannels2.rdata")

#select 30 random from each condition - same participants
combined_sampled<-sample(1:(length(list_eeg)/4),30)
duration_sampled<-combined_sampled+length(list_eeg)/4
frequency_sampled<-combined_sampled+2*length(list_eeg)/4
standard_sampled<-combined_sampled+3*length(list_eeg)/4

select_samples<-c(combined_sampled,
                  duration_sampled,
                  frequency_sampled,
                  standard_sampled)

#select thirty random participants
list_eeg_selected<-list_eeg[select_samples]

fun_prepare_channeldata_plotting<-function(one_set,sample_coordinates){

  #match names of sample coordnisates and data set
  names(one_set)<-toupper(names(one_set)) #capitalize names - naming convention of eegkit
  sample_coordinates<-sample_coordinates[rownames(sample_coordinates) %in% names(one_set),]
  sample_coordinates$variable<-rownames(sample_coordinates)

  one_set$TIMES<-as.factor(one_set$TIMES) #preserve in melting
  one_set<-reshape2::melt(one_set,value.name='amplitude') #convert to long format

  one_set<-merge(one_set,sample_coordinates,by='variable') #merge with sample coordniates
  names(one_set)[1]<-'electrode' #set name for later plotting function
  one_set$TIMES<-as.numeric(one_set$TIMES) #return to numeric
  return(one_set)

      }

list_eeg_selected<-pblapply(list_eeg_selected,fun_prepare_channeldata_plotting,sample_coordinates=eegcoord)
df_plot<-data.table::rbindlist(list_eeg_selected)

#scalp topography separated for condition
g_tm1<-topoplot(df_plot[df_plot$CONDITION=='standard' & df_plot$TIMES>250 & df_plot$TIMES<350,], #select time
         interp_limit='head',limits=c(-3,3))+labs(title='standard') #scale coordniates to head

g_tm2<-topoplot(df_plot[df_plot$CONDITION=='duration' & df_plot$TIMES>250 & df_plot$TIMES<350,], #select time
         interp_limit='head',limits=c(-3,3))+labs(title='length oddball') #scale coordniates to head

g_tm3<-topoplot(df_plot[df_plot$CONDITION=='frequenc' & df_plot$TIMES>250 & df_plot$TIMES<350,], #select time
         interp_limit='head',limits=c(-3,3))+labs(title='pitch oddball') #scale coordniates to head

g_tm4<-topoplot(df_plot[df_plot$CONDITION=='combined' & df_plot$TIMES>250 & df_plot$TIMES<350,], #select time
         interp_limit='head',limits=c(-3,3))+labs(title='pitch & length oddball') #scale coordniates to head


#prepare data for plotting
df_plot$times<-df_plot$TIMES-100
custom_condition_colors <- rev(wes_palette('FantasticFox1',5,type='discrete')[2:5]) #reverse custom colors to match color coding in other figures


##TODO: replace with section below
g_main<-ggplot(df_plot[df_plot$electrode=="FZ",],aes(x=times,y=amplitude,group=CONDITION,color=CONDITION,fill=CONDITION))+geom_smooth()+
  labs(y='electrode Fz - amplitute (uV)',x='trial duration (ms)')+
  scale_fill_manual(values = custom_condition_colors, labels=c("standard" = "standard", "frequenc" = "pitch oddball", "duration" = "length oddball ", "combined" = "pitch & length oddball"))+
  scale_color_manual(values = custom_condition_colors, labels=c("standard" = "standard", "frequenc" = "pitch oddball", "duration" = "length oddball ", "combined" = "pitch & length oddball"))


#ggplot(df_plot[df_plot$electrode %in% c('FZ','F1','F2','FC1','FC2'),],aes(x=times,y=amplitude,group=CONDITION,color=CONDITION,fill=CONDITION))+geom_smooth()

table(df_plot$electrode)

#extract legend
g_legend<-function(x){
  tmp <- ggplot_gtable(ggplot_build(x))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(g_tm1)


# #grid arrange plot
# grid.arrange(
#   g_main,
#   arrangeGrob(
#     g_tm1 + theme(legend.position="none"),
#     g_tm2 + theme(legend.position="none"),
#     g_tm3 + theme(legend.position="none"),
#     g_tm4 + theme(legend.position="none"),nrow=2),
#   mylegend, ncol=3, widths=c(7,5,1))


#save to file
tiff(file=paste0(project_path,"/output/figures/figure_mmn_eeg3.tiff"), # create a file in tiff format in current working directory
     width=12, height=5, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure

grid.arrange(
  g_main,
  arrangeGrob(
    g_tm1 + theme(legend.position="none"),
    g_tm2 + theme(legend.position="none"),
    g_tm3 + theme(legend.position="none"),
    g_tm4 + theme(legend.position="none"),nrow=2),
  mylegend, ncol=3, widths=c(7,5,1))

dev.off() #close operation and save file





###--> visualize Fz data ####

#load(file="C:/Users/nico/Desktop/leap_oddball_eeg_list.rdata")
load(paste0(project_path,"data/leap_oddball_eeg_list.rdata"))

df_eeg<-data.table::rbindlist(list_eeg) ##very fast compared to do.call(rbind)

subsample<-sample(1:nrow(df_eeg),nrow(df_eeg)/10)
custom_condition_colors <- rev(wes_palette('FantasticFox1',5,type='discrete')[2:5]) #reverse custom colors to match color coding in other figures

ggplot(df_eeg[subsample,],aes(x=times,y=Fz,group=condition,color=condition,fill=condition))+geom_smooth()+
  labs(y='electrode Fz - amplitute (uV)',x='trial duration (ms)')+
  scale_fill_manual(values = custom_condition_colors, labels=c("standard" = "standard", "frequenc" = "pitch oddball", "duration" = "length oddball ", "combined" = "pitch & length oddball"))+
  scale_color_manual(values = custom_condition_colors, labels=c("standard" = "standard", "frequenc" = "pitch oddball", "duration" = "length oddball ", "combined" = "pitch & length oddball"))


###---> calculate MMN ####

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

##--> create data frame
ID<-substr(names(list_mmn),10,22)
condition<-substr(names(list_mmn),1,8)

list_mmn<-pbmapply(function(x,y,z){
  x$condition<-y
  x$ID<-z
  return(x)
  },x=list_mmn,y=condition,z=ID,SIMPLIFY = F)

df_mmn<-data.table::rbindlist(list_mmn)

df_mmn$trial_counter<-as.numeric(df_mmn$trial_counter)
df_mmn$condition<-as.factor(df_mmn$condition)

### investigate mmn ####

require(lme4)
require(lmerTest)
require(emmeans)

##-investigate
with(df_mmn,by(mmn,condition,summary))
with(df_mmn,hist(mmn))
with(df_mmn,table(condition))

ggplot(df_mmn[df_mmn$mmn<10 & df_mmn$mmn>-20,],aes(x=mmn,group=condition,color=condition,fill=condition))+geom_density(alpha=0.2)


lmm<-lmer(scale(mmn)~condition*trial_counter+(1|ID),data=df_mmn)

anova(lmm)
plot(contrast(emmeans(lmm,~condition),'pairwise'))
fixef(lmm)['trial_counter']
###-->

    #create plot --> of contrasts
    plot_task_effect<-plot(contrast(emmeans(lmm,~condition),'pairwise'))[['data']]

    #modify plot
    ggplot(plot_task_effect,aes(x=contrast,y=the.emmean))+
      #conventional error bar
      geom_errorbar(aes(min=asymp.LCL,max=asymp.UCL),width=0.2)+
      geom_boxplot(aes(fill=contrast,
                       middle=the.emmean,
                       lower=the.emmean-1.5*SE,
                       upper=the.emmean+1.5*SE,
                       ymin=asymp.LCL,
                       ymax=asymp.UCL),stat = "identity",alpha=0.7)+
      #scale_fill_manual(values = custom_contrast_colors,labels = c('standard - pitch','standard - length','standard - pitch & length','pitch - length','pitch - pitch & length','length - pitch & length'))+
      scale_x_discrete(labels = NULL, breaks = NULL)+ #remove x-axis tick labels
      #coord_cartesian(ylim = c(-0.1, 0.1))+
      labs(x='contrast of task condition',y='MMN difference (z)')

### MERGE WITH df_trial ####

load("C:/Users/nico/PowerFolders/project_oddball_LEAP/data/mmn_leap_pd_final_dfs_21022023")

##prepare set
df_mmn$ID<-sub('/','',df_mmn$ID)
df_mmn$ID<-sub('_','',df_mmn$ID)

    table(df_mmn$ID)
    table(df_trial$id)
    table(df_trial$subjects)

    with(df_timepoint,table(subjects,wave))

    #compare data frame to merge
    hist(df_trial$EventCounter)
    hist(df_mmn$trial_counter)

    #eye-tracking data with MMN data
    table(df_timepoint$subjects %in% df_mmn$ID)
    table(unique(df_mmn$ID) %in% df_timepoint$subjects)


df_mmn$merge_id<-interaction(df_mmn$ID,df_mmn$trial_counter)
df_trial$merge_id<-interaction(df_trial$subjects,df_trial$EventCounter)

df_trial<-merge(df_trial,df_mmn,by='merge_id',all.x=T)

    table(is.na(df_mmn$mmn))
    table(is.na(df_trial$mmn))
    table(is.na(df_trial$rpd_auc))

with(df_trial,table(condition,EventData))
##some are not correctly matched
df_trial_mmn<-df_trial[(df_trial$EventData=='201' & df_trial$condition=='standard'|
                 df_trial$EventData=='203' & df_trial$condition=='duration'|
                 df_trial$EventData=='202' & df_trial$condition=='frequenc'|
                 df_trial$EventData=='204' & df_trial$condition=='combined'),]

###aggregate to df_timepoint

df_mmn_timepoint<-aggregate(mmn~id+EventData,data=df_trial_mmn,FUN=mean,na.rm=T)
df_mmn_timepoint<-reshape(df_mmn_timepoint, idvar = "id", timevar = "EventData", direction = "wide")

df_timepoint<-merge(df_timepoint,df_mmn_timepoint,by='id',all.x = T)

###which participants have MMN data
table(is.na(df_timepoint$mmn.201),df_timepoint$t1_diagnosis) ##

###calculate MMN as difference measures
df_timepoint$mmn_202_diff<-with(df_timepoint,mmn.202-mmn.201)
df_timepoint$mmn_203_diff<-with(df_timepoint,mmn.203-mmn.201)
df_timepoint$mmn_204_diff<-with(df_timepoint,mmn.204-mmn.201)

###ANALYSIS OF MMN ####
##### -- correlations on per trial level ####
with(df_trial_mmn[df_trial_mmn$EventData=='201',],cor.test(rpd_auc,mmn))
with(df_trial_mmn[df_trial_mmn$EventData=='202',],cor.test(rpd_auc,mmn))
with(df_trial_mmn[df_trial_mmn$EventData=='203',],cor.test(rpd_auc,mmn))
with(df_trial_mmn[df_trial_mmn$EventData=='204',],cor.test(rpd_auc,mmn))

with(df_trial_mmn[df_trial_mmn$EventData=='201',],cor.test(pd,mmn))
with(df_trial_mmn[df_trial_mmn$EventData=='202',],cor.test(pd,mmn))
with(df_trial_mmn[df_trial_mmn$EventData=='203',],cor.test(pd,mmn))
with(df_trial_mmn[df_trial_mmn$EventData=='204',],cor.test(pd,mmn))


#### -- correlations on per-participant level ####

#SEPR -PD
with(df_timepoint,cor.test(rpd_auc.201,pd))
with(df_timepoint,cor.test(rpd_auc.202,pd))
with(df_timepoint,cor.test(rpd_auc.203,pd)) #positive
with(df_timepoint,cor.test(rpd_auc.204,pd)) #negative

#SEPR - MMN --> uncorrelated
with(df_timepoint,cor.test(rpd_auc.201,mmn.201))
with(df_timepoint,cor.test(rpd_auc.202,mmn.202))
with(df_timepoint,cor.test(rpd_auc.203,mmn.203))
with(df_timepoint,cor.test(rpd_auc.204,mmn.204))

#BPS - MMN --> substantial negative correlations
with(df_timepoint,cor.test(mmn.201,pd))
with(df_timepoint,cor.test(mmn.202,pd))
with(df_timepoint,cor.test(mmn.203,pd))
with(df_timepoint,cor.test(mmn.204,pd))

# NG - BPS
with(df_timepoint,cor.test(gain_pcdm,pd))
with(df_timepoint,cor.test(gain_pcdm,pd))
with(df_timepoint,cor.test(gain_pcdm,pd))
with(df_timepoint,cor.test(gain_pcdm,pd))

# NG - SEPR
with(df_timepoint,cor.test(gain_pcdm,rpd_auc.201))
with(df_timepoint,cor.test(gain_pcdm,rpd_auc.202))
with(df_timepoint,cor.test(gain_pcdm,rpd_auc.203))
with(df_timepoint,cor.test(gain_pcdm,rpd_auc.204))

#NG - MMN
with(df_timepoint,cor.test(gain_pcdm,mmn.201))
with(df_timepoint,cor.test(gain_pcdm,mmn.202))
with(df_timepoint,cor.test(gain_pcdm,mmn.203))
with(df_timepoint,cor.test(gain_pcdm,mmn.204))

#Bonferoni correction for intitially significant correlations
number_of_comparisons<-24
as.numeric(with(df_timepoint,cor.test(rpd_auc.203,pd))['p.value'])*number_of_comparisons
as.numeric(with(df_timepoint,cor.test(rpd_auc.204,pd))['p.value'])*number_of_comparisons
as.numeric(with(df_timepoint,cor.test(mmn.201,pd))['p.value'])*number_of_comparisons
as.numeric(with(df_timepoint,cor.test(mmn.202,pd))['p.value'])*number_of_comparisons
as.numeric(with(df_timepoint,cor.test(mmn.203,pd))['p.value'])*number_of_comparisons
as.numeric(with(df_timepoint,cor.test(mmn.204,pd))['p.value'])*number_of_comparisons




#### -- linear models in aggregated data ####

#MMN
lm_mmn<-lm(scale(mmn_202_diff)~t1_diagnosis,df_timepoint)
anova(lm_mmn)

lm_mmn<-lm(scale(mmn_203_diff)~t1_diagnosis,df_timepoint)
anova(lm_mmn)

lm_mmn<-lm(scale(mmn_204_diff)~t1_diagnosis,df_timepoint)
anova(lm_mmn)


with(df_timepoint,by(mmn_202_diff,t1_diagnosis,summary))
with(df_timepoint,by(mmn_203_diff,t1_diagnosis,summary))

effectsize::cohens_d(x=scale(mmn_202_diff)~t1_diagnosis,data=df_timepoint)
effectsize::cohens_d(x=scale(mmn_203_diff)~t1_diagnosis,data=df_timepoint)
### rather stronger MMN in ASD


#### --- task effects on MMN unaggregated ####
lmm<-lmer(scale(mmn)~EventData*EventCounter+(1|subjects),df_trial_mmn)
anova(lmm)

fixef(lmm)['EventCounter'] ### increase in amplitude with trial counter
contrast(emmeans(lmm,~EventData),'pairwise') ### oddballs associated with lower MMN

      #create plot --> of contrasts
      plot_task_effect<-plot(contrast(emmeans(lmm,~EventData),'pairwise'))[['data']]

      #modify plot
      g3<-ggplot(plot_task_effect,aes(x=contrast,y=the.emmean))+
        #conventional error bar
        geom_errorbar(aes(min=asymp.LCL,max=asymp.UCL),width=0.2)+
        geom_boxplot(aes(fill=contrast,
                         middle=the.emmean,
                         lower=the.emmean-1.5*SE,
                         upper=the.emmean+1.5*SE,
                         ymin=asymp.LCL,
                         ymax=asymp.UCL),stat = "identity",alpha=0.7)+
        scale_fill_manual(values = custom_contrast_colors,labels = c('standard - pitch','standard - length','standard - pitch & length','pitch - length','pitch - pitch & length','length - pitch & length'))+
        scale_x_discrete(labels = NULL, breaks = NULL)+ #remove x-axis tick labels
        coord_cartesian(ylim = c(-0.1, 0.2))+
        labs(x='contrast of task conditions',y='MMN difference (z)')


              #SEPR -pupillary response - technical model
              lmm<-lmer(scale(rpd_auc)~EventData+
                          (1|subjects)+(1|wave),data=df_trial)

              #create plot --> of contrasts
              plot_task_effect<-plot(contrast(emmeans(lmm,~EventData),'pairwise'))[['data']]

              #modify plot
              g2<-ggplot(plot_task_effect,aes(x=contrast,y=the.emmean))+
                #conventional error bar
                geom_errorbar(aes(min=asymp.LCL,max=asymp.UCL),width=0.2)+

                #boxplot
                geom_boxplot(aes(fill=contrast,
                                 middle=the.emmean,
                                 lower=the.emmean-1.5*SE,
                                 upper=the.emmean+1.5*SE,
                                 ymin=asymp.LCL,
                                 ymax=asymp.UCL),stat = "identity",alpha=0.7)+
                scale_fill_manual(values = custom_contrast_colors,labels = c('standard - pitch','standard - length','standard - pitch & length','pitch - length','pitch - pitch & length','length - pitch & length'))+
                scale_x_discrete(labels = NULL, breaks = NULL)+ #remove x-axis tick labels
                coord_cartesian(ylim = c(-0.1, 0.2))+
                labs(x='contrast of task conditions',y='SEPR difference (z)')

              #BPS
              lmm<-lmer(scale(pd_baseline)~EventData+
                          (1|subjects)+(1|wave),data=df_trial)

              #create plot --> of contrasts
              plot_task_effect<-plot(contrast(emmeans(lmm,~EventData),'pairwise'))[['data']]

              #modify plot
              g1<-ggplot(plot_task_effect,aes(x=contrast,y=the.emmean))+
                #conventional error bar
                geom_errorbar(aes(min=asymp.LCL,max=asymp.UCL),width=0.2)+

                # #overplot predicted values
                # geom_jitter(data=df_plot_predicted_bps[complete.cases(df_plot_predicted_bps),],
                #             aes(x=predicted_bps_eventdata,
                #                 y=predicted_bps_mean,color=predicted_bps_eventdata),alpha=0.3,width=0.4,show.legend=F)+
                #boxplot
                geom_boxplot(aes(fill=contrast,
                                 middle=the.emmean,
                                 lower=the.emmean-1.5*SE,
                                 upper=the.emmean+1.5*SE,
                                 ymin=asymp.LCL,
                                 ymax=asymp.UCL),stat = "identity",alpha=0.7)+
                scale_fill_manual(values = custom_contrast_colors,labels = c('standard - pitch','standard - length','standard - pitch & length','pitch - length','pitch - pitch & length','length - pitch & length'))+
                scale_x_discrete(labels = NULL, breaks = NULL)+ #remove x-axis tick labels
                coord_cartesian(ylim = c(-0.1, 0.2))+
                labs(x='contrast of task conditions',y='BPS difference (z)')

              ###save figure

              #plot in markdown
              grid.arrange(g1+theme(legend.position="none"),g2+theme(legend.position="none"),g3,ncol=3,widths=c(1.2,1.2,2))

              tiff(file=paste0(project_path,"/output/figures/figure_taskcondition_effects_BPS_SEPR_MMN.tiff"), # create a file in tiff format in current working directory
                   width=10, height=4, units="in", res=300, compression='lzw') #define size and resolution of the resulting figure

              grid.arrange(g1+theme(legend.position="none"),g2+theme(legend.position="none"),g3,ncol=3,widths=c(1.2,1.2,2))

              dev.off() #close operation and save file


#### -- unaggregated effects of BPS and SEPR on MMN #####

lmm<-lmer(scale(mmn)~scale(pd_baseline)*EventData*EventCounter+(1|subjects),df_trial_mmn)
anova(lmm)

fixef(lmm)
###higer baseline with lower MMN amplitude

lmm<-lmer(scale(mmn)~scale(rpd_auc)*EventData*EventCounter+(1|subjects),df_trial_mmn)
anova(lmm)

fixef(lmm)
###higer respone with higher MMN amplitude


#### -- group differences ####
lmm<-lmer(scale(mmn)~EventData*t1_diagnosis*EventCounter+(1|subjects),df_trial_mmn)
anova(lmm)


confint(contrast(emtrends(lmm,~t1_diagnosis|EventData,var='EventCounter'),'pairwise'))
plot(emtrends(lmm,~t1_diagnosis|EventData,var='EventCounter'))
####--> effect directions follow rpd effects

lmm<-lmer(scale(mmn)~EventData*t1_diagnosis*sequence_position+(1|subjects),df_trial_mmn)
anova(lmm)



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

