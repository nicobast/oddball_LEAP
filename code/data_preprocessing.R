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


#detect system
ifelse(Sys.info()['sysname']=='Linux',
       home_path<-'~',
       home_path<-'C:/Users/Nico')

##search for the data
data_path<-'Z:/nico/backup_data/data_AIMS/mismatch_negativity/ET_data'
individual_data_folders<-list.dirs(data_path,recursive=T) #list data folders
individual_data_folders<-individual_data_folders[which(nchar(individual_data_folders)>(nchar(data_path)+6))] # remove topdirectories without data

#test matching of gaze and event data with one participant

test_folder<-individual_data_folders[1]

dir(test_folder)

#define id
length_pic<-12
id<-substr(test_folder,nchar(test_folder)-length_pic+1,nchar(test_folder))
wave<-substr(test_folder,nchar(test_folder)-length_pic-5,nchar(test_folder)-length_pic-1)
individual_id<-paste(id,wave,sep='_')

#read_data
events<-read(paste(test_folder,'session_events.csv',sep='/'))



