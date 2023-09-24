# This code is for comparison between MCMC and INLA
library(dplyr)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores()-1)
library(INLA)
#work_dir
setwd(work_dir)
source("./code/BiAR1 INLA model.R")
source("./code/simulation/test BiAR INLA STAN/Test BiAR1_LCM INLA STAN_user_define.R")
stan_file_path = "./code/simulation/test BiAR INLA STAN/stan code"
BiAR_LCM_stanmodel = stan_model(file =paste0(stan_file_path,"/BiAR_LCM.stan"))
n.sim = 50
N_list=c(10,20,50) ## number list of stocks
T_list=c(50,100,200) ## length list of time interval
# Data size setup
N_id=1
T_id=1
##################################################################
### hyperparameters for AR and LCM
##################################################################
# AR hyper
##################################################################
ar1_prec1 = 6; ar1_prec2 = 7
ar1_phi1 = 0.5; ar1_phi2 = 0.8; ar1_rho = 0.6 ## correlation between two ARs
##################################################################
# LCM hyper
alpha_prec1= 10;alpha_prec2 = 10; alpha_rho = 0.6
##################################################################
# Stan priors
wishart_prior_mat = diag(1,2) # Wishart prior scale matrix
wishart_nu =4 # wishar prior df
ar_rho_prior_sigma =2 # sigma for ar_rho prior distribution (normal)
beta0_prior_sigma = 10 # sigma for beta0 prior distribution (normal)



#### Simulation begins at here
N = N_list[N_id] # number of stock
T = T_list[T_id] # length of time interval
sim_idx =1
if(sim_idx==1){
  set.seed(05012023)
}
if(sim_idx==2){
  set.seed(05022023)
}
if(sim_idx==3){
  set.seed(05032023)
}
if(sim_idx==4){
  set.seed(05042023)
}

Intercepts = matrix(rnorm(2*N,mean=2,sd=0.1),ncol = 2) ## different intercepts for n in 1:N

CPU_INLA_vec=c()
CPU_STAN_vec =c()

Param_recover_INLA_vec =c()
Param_recover_STAN_vec =c()

Metrics_INLA_vec=c()
Metrics_STAN_vec =c()

for(n.sim_i in 1){
  list_temp = INLA_STAN_result(
    N,T,
    ar1_prec1 , ar1_prec2 ,
    ar1_phi1, ar1_phi2, ar1_rho , ## correlation between two ARs
    ##################################################################
    # LCM hyper
    alpha_prec1,alpha_prec2 , alpha_rho, 
    ##################################################################
    BiAR_LCM_stanmodel)
  
  CPU_INLA_temp = list_temp$list_INLA$CPU
  CPU_INLA_vec = c(CPU_INLA_vec,CPU_INLA_temp)
  
  Param_recover_INLA_temp = list_temp$list_INLA$Param_recover
  Param_recover_INLA_vec = rbind(Param_recover_INLA_vec,Param_recover_INLA_temp)
  
  Metrics_INLA_temp =list_temp$list_INLA$Metrics
  Metrics_INLA_vec  = rbind(Metrics_INLA_vec,Metrics_INLA_temp)

  CPU_STAN_temp = list_temp$list_STAN$CPU
  CPU_STAN_vec = c(CPU_STAN_vec,CPU_STAN_temp)
  
  Param_recover_STAN_temp = list_temp$list_STAN$Param_recover
  Param_recover_STAN_vec = rbind(Param_recover_STAN_vec,Param_recover_STAN_temp)
  
  Metrics_STAN_temp =list_temp$list_STAN$Metrics
  Metrics_STAN_vec  = rbind(Metrics_STAN_vec,Metrics_STAN_temp)
  }

CPU_INLA_perfomance = mean(CPU_INLA_vec)
Param_recover_INLA_performance = colMeans(Param_recover_INLA_vec)
Metrics_INLA_performance =colMeans(Metrics_INLA_vec)

CPU_STAN_perfomance = mean(CPU_STAN_vec)
Param_recover_STAN_performance = colMeans(Param_recover_STAN_vec)
Metrics_STAN_performance =colMeans(Metrics_STAN_vec)


INLA_STAN_performance_df=
  data.frame(INLA=c(CPU_STAN_perfomance,Param_recover_STAN_performance,Metrics_STAN_performance),
             STAN =c(CPU_INLA_perfomance,
                     Param_recover_INLA_performance,
                     Metrics_INLA_performance))
rownames(INLA_STAN_performance_df) =
  c("CPU","Intercepts","ar1_prec1","ar1_prec2","ar1_rho",
    "ar1_phi1","ar1_phi2","lcm_prec1","lcm_prec2","lcm_rho",
    "MAE","MSE","WMAPE")

 