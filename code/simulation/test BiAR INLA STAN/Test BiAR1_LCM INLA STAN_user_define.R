# This code is for comparison between MCMC and INLA
# Use-defined functions for iteration
# 11/09/2022
# Include INLA_only_result
############################################################################################
# This code is for comparison between MCMC and INLA
INLA_STAN_result = function(N,T,
                            ar1_prec1 , ar1_prec2 ,
                            ar1_phi1, ar1_phi2, ar1_rho , ## correlation between two ARs
                            ##################################################################
                            # LCM hyper
                            alpha_prec1,alpha_prec2 , alpha_rho, 
                            ##################################################################
                            BiAR_LCM_stanmodel){
  Intercepts = matrix(rnorm(2*N,mean=2,sd=0.1),ncol = 2) ## different intercepts for n in 1:N
  
  
  # Simulation for the bi-variate AR1
  stopifnot(T>1) # has to be greater than 1
  ## Starting point
  noise2d_s = mvtnorm::rmvnorm(n=1,
                               mean = rep(0,2),
                               sigma = matrix(
                                 c((ar1_prec1*(1-ar1_phi1**2))**-1,ar1_rho*(ar1_prec1*ar1_prec2)**-.5/(1-ar1_phi1*ar1_phi2),
                                   ar1_rho*(ar1_prec1*ar1_prec2)**-.5/(1-ar1_phi1*ar1_phi2),(ar1_prec2*(1-ar1_phi2**2))**-1),
                                 ncol=2
                               )
  ) 
  ## The rest of the noises
  noise2d_c=mvtnorm::rmvnorm(n=T-1,
                             mean = rep(0,2),
                             sigma = matrix(
                               c(ar1_prec1**-1,ar1_rho*ar1_prec1**-.5*ar1_prec2**-.5,
                                 ar1_rho*ar1_prec1**-.5*ar1_prec2**-.5,ar1_prec2**-1),ncol=2))
  
  noise2d = rbind(noise2d_s,noise2d_c)
  Phi_mat1=toeplitz(ar1_phi1**(0:(T-1)))
  Phi_mat1[upper.tri(Phi_mat1)]=0
  Phi_mat2=toeplitz(ar1_phi2**(0:(T-1)))
  Phi_mat2[upper.tri(Phi_mat2)]=0
  x1t =Phi_mat1 %*%noise2d[,1]
  x2t =Phi_mat2 %*%noise2d[,2]
  
  
  
  
  
  Sigma_alpha = diag(c(alpha_prec1,alpha_prec2)**-0.5)%*%
    matrix(c(1,alpha_rho,alpha_rho,1),ncol=2)%*%
    diag(c(alpha_prec1,alpha_prec2)**-0.5)
  alpha_st = mvtnorm::rmvnorm(N*T,mean = rep(0,2),sigma = Sigma_alpha)
  
  
  Y_para_tbl = dplyr::tibble(
    Intercept =c(apply(Intercepts,2,rep,each=T)),
    alpha = c(alpha_st),
    ar1 = c(rep(x1t,N),rep(x2t,N)),
    lambda = exp(Intercept+alpha+ar1),
    Y = rpois(2*N*T,lambda = lambda),
    t_idx = rep(1:T,2*N),
    J_idx = rep(1:2,each  =N*T),
    s_idx =rep(rep(1:N,each=T),2),
    alpha_idx = 1:(N*2*T)
  )
  
  
  # data frame restructuring
  Y_for_inla = matrix(NA,ncol = 2,nrow = 2*N*T)
  Y_for_inla[1:(N*T),1] = Y_para_tbl$Y[1:(N*T)]
  Y_for_inla[(1+N*T):(2*N*T),2] = Y_para_tbl$Y[(1+N*T):(2*N*T)]
  
  X_mat = matrix(NA,ncol = 2*N, nrow = 2*N*T) ## Separate the intercepts by N and alpha's dimension
  for(J in 1:2){
    for(s in 1:N){
      X_mat[(1:T)+((J-1)*N+s-1)*T,(J-1)*N+s]=1
    }
  }
  colnames(X_mat) = paste0(rep(c('intercept_J1','intercept_J2'),each = N),"_S",rep(1:N,2))
  
  t_idx_mat = matrix(NA,ncol = 2,nrow = 2*N*T)
  t_idx_mat[1:(N*T),1] = rep(1:T,N)
  t_idx_mat[(N*T+1):(2*N*T),2] = rep(1:T,N)
  colnames(t_idx_mat) = paste0('t_idx_J',1:2)
  
  dt_tbl = dplyr::as_tibble(cbind(X_mat,t_idx_mat))
  dt_tbl$alpha_idx = 1:(2*N*T)
  dt_tbl$s_idx = rep(rep(1:N,each=T),2)
  dt_tbl$BiAR1_time_idx  = c(rep(1:T,N),rep(1:T+T,N))
  
  
  formula_LCM = as.formula(paste0("Y_for_inla~-1+",paste0(colnames(X_mat),collapse = "+"),
                                  "+","f(BiAR1_time_idx,model = BiAR1.model)",
                                  "+","f(alpha_idx,model = 'iid2d',n=2*N*T,
                                hyper=list(prec1=list(param=c(4,1,1,0))))"))
  
  BiAR1.model = inla.rgeneric.define(inla.rgeneric.BiAR1.model,T=T)
  lcm_inla = INLA::inla(formula_LCM,data= dt_tbl,family = c('poisson',"poisson"),
                        control.compute=list(config=TRUE))
  #################################################################################################
  # STAN fitting
  data0= list(
    N=N,
    T=T,
    Y= Y_para_tbl$Y,
    wishart_prior_mat = wishart_prior_mat,
    wishart_nu = wishart_nu,
    ar_rho_prior_sigma =ar_rho_prior_sigma,
    beta0_prior_sigma = beta0_prior_sigma
  )
  
  stan_sampling_time = system.time({
    samples_temp = sampling(BiAR_LCM_stanmodel,data = data0,
                            pars =c('beta0',
                                    'Bi_AR_prec1',"Bi_AR_prec2","Bi_AR_rho",
                                    'ar_rho','LCM_prec1',"LCM_prec2","LCM_rho","Y_fitted"),chains=1, iter = 2000,warmup = 1000)
  })
  
  #################################################################################################
  # Begin to create INLA and STAN metrics here
  
  # INLA cpu
  CPU_INLA = lcm_inla$cpu.used[4]
  
  # INLA parameters recover
  ## Intercept
  fixed_contain = lcm_inla$summary.fixed
  fixed_contain$true = c(Intercepts)
  fixed_contain$contain = fixed_contain$`0.025quant`<fixed_contain$true & fixed_contain$`0.975quant`>fixed_contain$true
  fixed_contain_rate = mean(fixed_contain$contain) 
  
  ## Hypers
  alpha_m = c()
  fun_list = list(f1=exp,
                  f2=exp,
                  f3= function(x) 1-2/(exp(x)+1),
                  f4=function(x) 1-2/(exp(x)+1),
                  f5=function(x) 1-2/(exp(x)+1)
  )
  for(i in 1:5){
    zmargin_temp =  inla.zmarginal(inla.tmarginal(fun = fun_list[[i]],lcm_inla$marginals.hyperpar[[i]]),silent = TRUE)
    alpha_temp = c(zmargin_temp$mean,zmargin_temp$sd,
                   zmargin_temp$quant0.025,zmargin_temp$quant0.975)
    alpha_m = rbind(alpha_m,alpha_temp)
  }
  
  param_test=cbind(c(ar1_prec1,ar1_prec2,ar1_rho,ar1_phi1,ar1_phi2),
                   alpha_m)
  param_tbl_BiAR = data.frame(
    true = param_test[,1],
    est = param_test[,2],
    sd = param_test[,3],
    L = param_test[,4],
    U = param_test[,5]
  )
  
  param_tbl_LCM = lcm_inla$summary.hyperpar[-(1:5),c(1,2,3,5)]
  param_tbl_LCM = cbind(c(alpha_prec1,alpha_prec2,alpha_rho),param_tbl_LCM)
  colnames(param_tbl_LCM) = c("true",'est',"sd","L","U")
  param_tbl_overall = rbind(param_tbl_BiAR,param_tbl_LCM)
  
  param_tbl_overall$contain = (param_tbl_overall$true>param_tbl_overall$L)&(param_tbl_overall$true<param_tbl_overall$U)
  
  Param_recover_INLA = c(fixed_contain_rate,param_tbl_overall$contain)
  names(Param_recover_INLA) =c("Intercepts","ar1_prec1","ar1_prec2","ar1_rho",
                               "ar1_phi1","ar1_phi2","lcm_prec1","lcm_prec2","lcm_rho")
  
  ## MAE, MSE, WMAPE of fitted values
  Metrics_INLA =c(mean(abs(lcm_inla$summary.fitted.values$mean - Y_para_tbl$Y)),
                  mean((lcm_inla$summary.fitted.values$mean - Y_para_tbl$Y)^2),
                  sum(abs(Y_para_tbl$Y-lcm_inla$summary.fitted.values$mean))/sum(abs(Y_para_tbl$Y)))
  names(Metrics_INLA) =c("MAE","MSE","WMAPE")
  
  ## Parameters inference metrics for INLA without Intercept
  param_metrics_INLA=data.frame(
    param_name = c('ar1_prec1','ar1_prec2','ar1_rho','ar1_phi1','ar1_phi2',
                   'alpha_prec1','alpha_prec2','alpha_rho'),
    BIAS_squared = (param_tbl_overall$true-param_tbl_overall$est)^2,
    VAR = (param_tbl_overall$sd)^2,
    MSE = sqrt((param_tbl_overall$true-param_tbl_overall$est)^2 + (param_tbl_overall$sd)^2),
    Para_Est = param_tbl_overall$est,
    Para_True =  param_tbl_overall$true,
    Para_Est_L = param_tbl_overall$L,
    Para_Est_U = param_tbl_overall$U
  )

  
  list_INLA = list(CPU = CPU_INLA,
                   Param_recover = Param_recover_INLA,
                   param_metrics_INLA = param_metrics_INLA,
                   Metrics = Metrics_INLA) 
  
  # STAN cpu
  CPU_STAN = stan_sampling_time[3]
  
  # STAN parameters recover
  ## Intercept
  fixed_contain_STAN = as.data.frame(summary(samples_temp)$summary[1:(N*2),])
  fixed_contain_STAN$true = c(Intercepts)
  fixed_contain_STAN$contain = fixed_contain_STAN$`2.5%`<fixed_contain_STAN$true & fixed_contain_STAN$`97.5%`>fixed_contain_STAN$true
  fixed_contain_rate_STAN = mean(fixed_contain_STAN$contain)
  
  ## Hypers
  hyper_contain_STAN = as.data.frame(summary(samples_temp)$summary[2*N+(1:8),])
  hyper_contain_STAN$true = c(ar1_prec1,ar1_prec2,ar1_rho,ar1_phi1,ar1_phi2,
                              alpha_prec1,alpha_prec2,alpha_rho)
  hyper_contain_STAN$contain = hyper_contain_STAN$`2.5%` < hyper_contain_STAN$true &
    hyper_contain_STAN$`97.5%`> hyper_contain_STAN$true
  
  Param_recover_STAN = c(fixed_contain_rate,hyper_contain_STAN$contain)
  names(Param_recover_STAN) =c("Intercepts","ar1_prec1","ar1_prec2","ar1_rho",
                               "ar1_phi1","ar1_phi2","lcm_prec1","lcm_prec2","lcm_rho")
  
  ## MAE, MSE, WMAPE of fitted
  Y_fitted_STAN = as.data.frame(summary(samples_temp)$summary[2*N+8+1:(2*N*T),])
  Metrics_STAN =c(mean(abs(Y_fitted_STAN$mean - Y_para_tbl$Y)),
                  mean((Y_fitted_STAN$mean - Y_para_tbl$Y)^2),
                  sum(abs(Y_para_tbl$Y-Y_fitted_STAN$mean))/sum(abs(Y_para_tbl$Y)))
  names(Metrics_STAN) =c("MAE","MSE","WMAPE")
  
  
  ## Parameters inference metrics for STAN without intercepts
  param_metrics_STAN=data.frame(
    param_name = c('ar1_prec1','ar1_prec2','ar1_rho','ar1_phi1','ar1_phi2',
                   'alpha_prec1','alpha_prec2','alpha_rho'),
    BIAS_squared = (hyper_contain_STAN$mean-hyper_contain_STAN$true)^2,
    VAR = (hyper_contain_STAN$sd)^2,
    MSE = sqrt((hyper_contain_STAN$mean-hyper_contain_STAN$true)^2 + (hyper_contain_STAN$sd)^2),
    Para_Est = hyper_contain_STAN$mean,
    Para_True =  hyper_contain_STAN$true,
    Para_Est_L = hyper_contain_STAN$`2.5%`,
    Para_Est_U = hyper_contain_STAN$`97.5%`
  )

  
  list_STAN = list(CPU = CPU_STAN,
                   Param_recover = Param_recover_STAN,
                   param_metrics_STAN = param_metrics_STAN,
                   Metrics = Metrics_STAN) 
  
  
  contain_list = list(list_INLA=list_INLA,list_STAN=list_STAN)
  return(contain_list)
}







