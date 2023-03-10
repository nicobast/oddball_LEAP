install.packages('workflowr')
require(workflowr)

# apply git information for version control
wflow_git_config(user.name = "nicobast", user.email = "nico.bast@kgu.de")

#detect system
ifelse(Sys.info()['sysname']=='Linux',
       home_path<-'~',
       home_path<-'C:/Users/Nico')

#create initial files
wflow_start(paste(home_path,"/PowerFolders/project_oddball_LEAP/"))

#build --> docs ar ebuild and can be viewed on local machine
wflow_build()

wflow_view()

# #publish --> can be viewed by others
wflow_publish(c("analysis/index.Rmd", "analysis/about.Rmd", "analysis/license.Rmd",
                "analysis/extract_salience.Rmd",
                "analysis/data_preprocessing.Rmd",
                "analysis/data_analysis_salience.Rmd"),
              "Publish the initial files for myproject")

# #deploy to github
#wflow_use_github("nicobast") #initiate project on github
wflow_git_push() #deploy updates