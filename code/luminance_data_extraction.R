require(av) # split video to images
require(png)


video_location<-"C:/Users/nico/Desktop/mmnpingu_silent.m4v"
image_target_location<-"C:/Users/nico/Desktop/pingu_single_frames"
images_per_second<-3
mean_trial_duration<-0.55

#convert video to frames
av_media_info(video_location)
av_video_images(video_location,
                destdir = image_target_location,
                format = 'png',
                fps= images_per_second)


#load RGB images - required files
rgb_data_paths<-list.files(path=image_target_location, full.names=T)
rgb_data_names<-list.files(path=image_target_location)

  rgb_data<-readPNG(rgb_data_paths[1],native=F,info=T)
  hist(unlist(rgb_data))


#loop across all images of a video
luminance_data_of_one_video<-as.list(0)
for(j in 1:length(rgb_data_paths)){
  #read one image
  rgb_data<-readPNG(rgb_data_paths[j],native=F,info=F)

  #step2 - convert relative values to a linear value (RGB values are gamma-encoded with a power curve):
  # see: https://en.wikipedia.org/wiki/Relative_luminance
  rgb_linear_values<-rgb_data^2.2 # 2.2 power curve

  #step3 - caculate LUMINANCE by sRGB coefficient:
  luminance_image<-(rgb_linear_values[,,1]*0.2126)+(rgb_linear_values[,,2]*0.7152)+(rgb_linear_values[,,3]*0.0722)

  #mean luminance of an image:
  luminance_image<-c(mean(luminance_image),sd(luminance_image))
  luminance_data_of_one_video[[j]]<-luminance_image #complete information
  print(paste0(' image: ',j))

}

##convert to data.frame
luminance_m<-sapply(luminance_data_of_one_video, function(x){x<-x[1]})
luminance_sd<-sapply(luminance_data_of_one_video, function(x){x<-x[2]})
df_luminance<-data.frame(luminance_m,luminance_sd)
df_luminance$video_time<-seq(1:nrow(df_luminance))/images_per_second
df_luminance$trials<-df_luminance$video_time/mean_trial_duration

ggplot(df_luminance,aes(x=trials,y=luminance_m,ymin=luminance_m-1.96*luminance_sd,ymax=luminance_m+1.96*luminance_sd),
       color='darkblue')+geom_smooth(fill='lightblue',stat = "identity")+
  labs(x='trial number',y='relative luminance')

ggplot(df_luminance,aes(x=trials,y=luminance_m),
       color='darkblue')+geom_smooth(fill='lightblue')+
  labs(x='trial number',y='relative luminance')


saveRDS(df_luminance,'C:/Users/nico/PowerFolders/project_oddball_LEAP/data/stimulus_luminance.rds')
