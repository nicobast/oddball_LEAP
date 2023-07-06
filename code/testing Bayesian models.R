##define parameters of posterior sampling in Bayesian modeling
number_of_warmups <- 4000
number_of_iterations <- 8000

#BPS Bayeisan ####

# A. translate data - make greta arrays
y <- as_data(scale(df_trial$pd_baseline))
design <- as.matrix(cbind(ifelse(df_trial$t1_diagnosis=='ASD',T,F),df_trial$EventCounter))
subject_id <- as.integer(factor(df_trial$subjects)) #random intercept (in LMM)

#define priors
int <- normal(0, 1)
design_coef<- normal(0, 1, dim=ncol(design))
sd <- cauchy(0, 3, truncation = c(0, Inf))
#sd <- student(3, 0, 1, truncation = c(0, Inf))
a_subject <- normal(0, 1, dim = max(subject_id))

#OPERATIONS - MODEL
mean <- int + design %*% design_coef + a_subject[subject_id]

# likelihood - distribution over data
distribution(y) <- normal(mean, sd)

# defining the model
m <- model(design_coef)

# sampling
start <-Sys.time()
draws <- mcmc(m, n_samples = number_of_iterations, warmup = number_of_warmups,) ###--> test
Sys.time() - start

#stability of estimates
mcmc_trace(draws) #see whether chains converge and are stationary
mcmc_dens(draws) # + stable estimates

#analyse
summary(draws)
plot(draws)

#posterior predictive checking
coda::gelman.diag(draws) #Rhat <= 1.01 is cutoff for reliable estimator
coda::effectiveSize(draws) #ESS

#save model
draws_BPS_twoway<-draws
save(draws_BPS_twoway,file=paste0(project_path,'/data/Bayesian_model_BPS_twoway.Rdata'))

#BPS habituation Bayesian ####

df_bps_bayesian_habituation<-df_trial[df_trial$EventData==201 & df_trial$sequence_position<10,]

# A. translate data - make greta arrays
y <- as_data(scale(df_bps_bayesian_habituation$pd_baseline))
design <- as.matrix(cbind(ifelse(df_bps_bayesian_habituation$t1_diagnosis=='ASD',T,F),df_bps_bayesian_habituation$sequence_position))
subject_id <- as.integer(factor(df_bps_bayesian_habituation$subjects)) #random intercept (in LMM)

#define priors
int <- normal(0, 1)
design_coef<- normal(0, 1, dim=ncol(design))
sd <- cauchy(0, 3, truncation = c(0, Inf))
#sd <- student(3, 0, 1, truncation = c(0, Inf))
a_subject <- normal(0, 1, dim = max(subject_id))

#OPERATIONS - MODEL
mean <- int + design %*% design_coef + a_subject[subject_id]

# likelihood - distribution over data
distribution(y) <- normal(mean, sd)

# defining the model
m <- model(design_coef)

# sampling
start <-Sys.time()
draws <- mcmc(m, n_samples = number_of_iterations, warmup = number_of_warmups,) ###--> test
Sys.time() - start

#stability of estimates
mcmc_trace(draws) #see whether chains converge and are stationary
mcmc_dens(draws) # + stable estimates

#analyse
summary(draws)
plot(draws)

#posterior predictive checking
coda::gelman.diag(draws) #Rhat <= 1.01 is cutoff for reliable estimator
coda::effectiveSize(draws) #ESS

#save model
draws_BPS_supp<-draws
save(draws_BPS_supp,file=paste0(project_path,'/data/Bayesian_model_BPS_supp.Rdata'))

#SEPR Bayesian ####

# A. translate data - make greta arrays
y <- as_data(scale(df_trial$rpd_auc))
trials <- as_data(scale(df_trial$EventCounter))
subject_id <- as.integer(factor(df_trial$subjects)) #random intercept (in LMM)
events <- as_data(model.matrix(~ EventData - 1, df_trial)) #alternative to dummy coding?
interaction_variable<-with(df_trial,interaction(t1_diagnosis,EventData))
group_by_events <- as_data(model.matrix(~interaction_variable-1, df_trial)) #alternative to dummy coding - remove intercept (-1)
#table(group_by_events[,1],group_by_events[,9])

#define priors
event_coefs<- normal(0, 1, dim=ncol(events))
group_event_coefs <- normal(0,1, dim=ncol(group_by_events))

sd <- cauchy(0, 3, truncation = c(0, Inf))
a_subject <- normal(0, 1, dim = max(subject_id))

#OPERATIONS - MODEL
mean <- trials * group_by_events %*% group_event_coefs + a_subject[subject_id]

# likelihood - distribution over data
distribution(y) <- normal(mean, sd)

# defining the model
#m <- model(int, coef_201, coef_202, coef_203, coef_204, sd)
m <- model(group_event_coefs)

# sampling
start <-Sys.time()
draws <- mcmc(m, n_samples = number_of_iterations, warmup = number_of_warmups) ###--> test
Sys.time() - start


#estimator convergence
mcmc_trace(draws) #see whether chains converge and are stationary
mcmc_dens(draws) # + stable estimates

##results
mcmc_areas(draws, prob= .66, prob_outer = .95)+scale_y_discrete(labels=c('ASD-standard','TD-standard ','ASD-pitch','TD-pitch','ASD-length ','TD-length','ASD-pitch+length','TD-pitch+length'))+coord_flip()+theme_bw()

#analyse
summary(draws)


#posterior predictive checking
coda::gelman.diag(draws) #Rhat <= 1.01 is cutoff for reliable estimator
coda::effectiveSize(draws) #ESS

summary(draws_SEPR_threeway)

#save model
draws_SEPR_threeway<-draws
save(draws_SEPR_threeway,file=paste0(project_path,'/data/Bayesian_model_SEPR_threewayinteraction.Rdata'))


#MMN Bayesian ####

#remove NA
df_bayesian_mmn<-df_trial[!is.na(df_trial$mmn),]

# A. translate data - make greta arrays
y <- as_data(scale(df_bayesian_mmn$mmn))
trials <- as_data(scale(df_bayesian_mmn$EventCounter))
subject_id <- as.integer(factor(df_bayesian_mmn$subjects)) #random intercept (in LMM)
events <- as_data(model.matrix(~ EventData - 1, df_bayesian_mmn)) #alternative to dummy coding?
interaction_variable<-with(df_bayesian_mmn,interaction(t1_diagnosis,EventData))
group_by_events <- as_data(model.matrix(~interaction_variable-1, df_bayesian_mmn)) #alternative to dummy coding - remove intercept (-1)
#table(group_by_events[,1],group_by_events[,9])

#define priors
event_coefs<- normal(0, 1, dim=ncol(events))
group_event_coefs <- normal(0,1, dim=ncol(group_by_events))

sd <- cauchy(0, 3, truncation = c(0, Inf))
a_subject <- normal(0, 1, dim = max(subject_id))

#OPERATIONS - MODEL
mean <- trials * group_by_events %*% group_event_coefs + a_subject[subject_id]

# likelihood - distribution over data
distribution(y) <- normal(mean, sd)

# defining the model
#m <- model(int, coef_201, coef_202, coef_203, coef_204, sd)
m <- model(group_event_coefs)

# sampling
start <-Sys.time()
draws <- mcmc(m, n_samples = number_of_iterations, warmup = number_of_warmups) ###--> test
Sys.time() - start


#estimator convergence
mcmc_trace(draws) #see whether chains converge and are stationary
mcmc_dens(draws) # + stable estimates

##results
mcmc_areas(draws, prob= .66, prob_outer = .95)+scale_y_discrete(labels=c('ASD-standard','TD-standard ','ASD-pitch','TD-pitch','ASD-length ','TD-length','ASD-pitch+length','TD-pitch+length'))+coord_flip()+theme_bw()

#analyse
summary(draws)


#posterior predictive checking
coda::gelman.diag(draws) #Rhat <= 1.01 is cutoff for reliable estimator
coda::effectiveSize(draws) #ESS - the higher the better

#save model
draws_MMN_threeway<-draws
save(draws_MMN_threeway,file=paste0(project_path,'/data/Bayesian_model_MMN_threewayinteraction.Rdata'))
