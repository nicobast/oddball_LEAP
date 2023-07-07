###per participant models

###--> table with t value, df, p, and R²

###t-test ####
#BPS
summary(lm<-lm(scale(pd)~t1_diagnosis,df_timepoint))
cohens_d(x=scale(pd)~t1_diagnosis,data=df_timepoint)


#NG
summary(lm(scale(gain_pcdm)~t1_diagnosis,df_timepoint[df_timepoint$Rsq_pcdm>0.2,]))
cohens_d(x=scale(pd)~t1_diagnosis,data=df_timepoint[df_timepoint$Rsq_pcdm>0.2,])


#SEPR
summary(lm(scale(rpd_auc.201)~t1_diagnosis,df_timepoint))
summary(lm(scale(rpd_auc.202)~t1_diagnosis,df_timepoint))
summary(lm(scale(rpd_auc.203)~t1_diagnosis,df_timepoint))
summary(lm(scale(rpd_auc.204)~t1_diagnosis,df_timepoint))

#MMN
summary(lm(scale(mmn.201)~t1_diagnosis,df_timepoint))
summary(lm(scale(mmn.202)~t1_diagnosis,df_timepoint))
summary(lm(scale(mmn.203)~t1_diagnosis,df_timepoint))
summary(lm(scale(mmn.204)~t1_diagnosis,df_timepoint))


###t-ANCOVA - SQ SUM type 3 ####
#BPS
lm<-lm(scale(pd)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
         ,df_timepoint)
car::Anova(lm,type=3)

#SEPR
lm<-lm(scale(rpd_auc.201)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)

lm<-lm(scale(rpd_auc.202)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)

lm<-lm(scale(rpd_auc.203)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)

lm<-lm(scale(rpd_auc.204)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)

#MMN
lm<-lm(scale(mmn.201)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)

lm<-lm(scale(mmn.202)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)

lm<-lm(scale(mmn.203)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)

lm<-lm(scale(mmn.204)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)
car::Anova(lm,type=3)


# significant covariate models ####
#full model - PD
lm<-lm(scale(pd)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint)

  ##supplements - full table
  row_names<-c('intercept','group',
               'age','perceptual IQ','sex',
               'sampling rate','gaze center deviation','data quality',
               'residuals')


  table_data_model<-round(Anova(lm,type=3),3)
  table_data_model<-cbind(row_names,table_data_model)

  table_data_model

  table_model<-table_data_model %>%
    kbl(caption = "Per participant: Full Linear model of baseline pupil size (BPS) with all potential covariates",
        col.names = c('','Sum Sq','df','F','p'),
        row.names = F) %>%
    kable_classic(full_width = F, html_font = "Cambria")

  table_model

  save_kable(table_model, file= paste0(project_path,'/output/supplements/table_linearmodel_BPS_allcovariates.html'))

#significant covariate model - PD
lm<-lm(scale(pd)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)
       ,df_timepoint)

  ##supplements - full table
  row_names<-c('intercept','group',
               'age','perceptual IQ','residuals')


  table_data_model<-round(Anova(lm,type=3),3)
  table_data_model<-cbind(row_names,table_data_model)

  table_data_model

  table_model<-table_data_model %>%
    kbl(caption = "Per participant: Reduced Linear model of baseline pupil size (BPS) with significant covariates",
        col.names = c('','Sum Sq','df','F','p'),
        row.names = F) %>%
    kable_classic(full_width = F, html_font = "Cambria")

  table_model

  save_kable(table_model, file= paste0(project_path,'/output/supplements/table_linearmodel_BPS_signcovariates.html'))



# full modelNG
lm<-lm(scale(gain_pcdm)~t1_diagnosis+
         scale(ageyrs)+scale(t1_piq)+sex+
         as.factor(sampling_rate)+scale(center_dev)+scale(missing_data_trial)
       ,df_timepoint[df_timepoint$Rsq_pcdm>0.2,])

  ##supplements - full table
  row_names<-c('intercept','group',
               'age','perceptual IQ','sex',
               'sampling rate','gaze center deviation','data quality',
               'residuals')


  table_data_model<-round(Anova(lm,type=3),3)
  table_data_model<-cbind(row_names,table_data_model)

  table_data_model

  table_model<-table_data_model %>%
    kbl(caption = "Per participant: Full Linear model of neural gain (NG) with all potential covariates",
        col.names = c('','Sum Sq','df','F','p'),
        row.names = F) %>%
    kable_classic(full_width = F, html_font = "Cambria")

  table_model

  save_kable(table_model, file= paste0(project_path,'/output/supplements/table_linearmodel_NG_allcovariates.html'))



#significant covariate model - NG
lm<-lm(scale(gain_pcdm)~t1_diagnosis
       ,df_timepoint[df_timepoint$Rsq_pcdm>0.2,])

    ##supplements - full table
    row_names<-c('intercept','group','residuals')


    table_data_model<-round(Anova(lm,type=3),3)
    table_data_model<-cbind(row_names,table_data_model)

    table_data_model

    table_model<-table_data_model %>%
      kbl(caption = "Per participant: Reduced Linear model of neural gain (NG) with significant covariates",
          col.names = c('','Sum Sq','df','F','p'),
          row.names = F) %>%
      kable_classic(full_width = F, html_font = "Cambria")

    table_model

    save_kable(table_model, file= paste0(project_path,'/output/supplements/table_linearmodel_NG_signcovariates.html'))




### table ####

BPS_group<-c(summary(lm<-lm(scale(pd)~t1_diagnosis,df_timepoint))[['fstatistic']],
  summary(lm<-lm(scale(pd)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])

NG_group<-c(summary(lm<-lm(scale(gain_pcdm)~t1_diagnosis,df_timepoint[df_timepoint$Rsq_pcdm>0.2,]))[['fstatistic']],
             summary(lm<-lm(scale(gain_pcdm)~t1_diagnosis,df_timepoint[df_timepoint$Rsq_pcdm>0.2,]))[['coefficients']][2,c(1,2,4)])

SEPR_201_group<-c(summary(lm<-lm(scale(rpd_auc.201)~t1_diagnosis,df_timepoint))[['fstatistic']],
             summary(lm<-lm(scale(rpd_auc.201)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])
SEPR_202_group<-c(summary(lm<-lm(scale(rpd_auc.202)~t1_diagnosis,df_timepoint))[['fstatistic']],
                  summary(lm<-lm(scale(rpd_auc.202)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])
SEPR_203_group<-c(summary(lm<-lm(scale(rpd_auc.203)~t1_diagnosis,df_timepoint))[['fstatistic']],
                  summary(lm<-lm(scale(rpd_auc.203)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])
SEPR_204_group<-c(summary(lm<-lm(scale(rpd_auc.204)~t1_diagnosis,df_timepoint))[['fstatistic']],
                  summary(lm<-lm(scale(rpd_auc.204)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])

MMN_201_group<-c(summary(lm<-lm(scale(mmn.201)~t1_diagnosis,df_timepoint))[['fstatistic']],
                  summary(lm<-lm(scale(mmn.201)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])
MMN_202_group<-c(summary(lm<-lm(scale(mmn.202)~t1_diagnosis,df_timepoint))[['fstatistic']],
                  summary(lm<-lm(scale(mmn.202)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])
MMN_203_group<-c(summary(lm<-lm(scale(mmn.203)~t1_diagnosis,df_timepoint))[['fstatistic']],
                  summary(lm<-lm(scale(mmn.203)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])
MMN_204_group<-c(summary(lm<-lm(scale(mmn.204)~t1_diagnosis,df_timepoint))[['fstatistic']],
                  summary(lm<-lm(scale(mmn.204)~t1_diagnosis,df_timepoint))[['coefficients']][2,c(1,2,4)])

group_comparison_table<-round(rbind(BPS_group,NG_group,
      SEPR_201_group,SEPR_202_group,SEPR_203_group,SEPR_204_group,
      MMN_201_group,MMN_202_group,MMN_203_group,MMN_204_group),3)


row_names<-c('BPS - across stimuli','NG - standards',
            'SEPR - standards','SEPR - pitch oddball','SEPR - length oddball','SEPR - pitch+length oddball',
            'MMN - standards','MMN - pitch oddball','MMN - length oddball','MMN - pitch+length oddball')

group_comparison_table<-cbind(row_names,group_comparison_table)

#produce table
group_comparison_table_format<-group_comparison_table %>%
  kbl(caption = "Per-participant level: Group comparisons of baseline pupil size (BPS), stimulus-evoked pupillary response (SEPR), mismatch negativity (MMN), and neural gain (NG)",
      row.names = F, col.names = c('','F','df1','df2','b','std. error','p')) %>%
  kable_classic(full_width = F, html_font = "Cambria")

group_comparison_table_format

save_kable(group_comparison_table_format,file=paste0(project_path,'/output/supplements/table_groupcomparison_perparticipant.html'))



