
Simulate_dataset_group_and_individual=function(J=2,N=13,T=200,rho=.8,sig11=.64,sig22=.64){
  
  # length(beta0)
  err_temporal =c(arima.sim(model = list(order=c(0,1,0),sd=.001),n=T-1),arima.sim(model = list(order=c(0,1,0),sd=.001),n=T-1))
  err_temporal=rep(err_temporal,each=N)# repeated
  # length(err_temporal)
  # matplot(1:T,err_temporal,type='l') # The random walk component
  
  
  # correlation matrix for LCM
  # rho=.8;sig11=sig22=.64
  mat_alpha=matrix(c(sig11,rho*sqrt(sig11*sig22),rho*sqrt(sig11*sig22),sig22),nrow=2)
  # mat_alpha
  
  ## Case 1 grouped level effect
  err_alpha_grouped=c(mvtnorm::rmvnorm(T,mean = rep(0,2),sigma = mat_alpha))
  err_alpha_grouped=rep(err_alpha_grouped,each=N) # repeated
  # length(err_alpha_grouped)
  Y_index = rep(1:N,T)
  time_index=rep(1:T,each=N)
  Y12_grouped=cbind(Y_index,time_index,matrix(rpois(N*T*J,lambda = exp(err_temporal+err_alpha_grouped)),ncol=J))
  
  
  
  
  ## Case 2 non-grouped level effect
  err_alpha_individual=c(mvtnorm::rmvnorm(T*N,mean = rep(0,2),sigma = mat_alpha))
  # length(err_alpha_individual)
  
  Y12_individual=cbind(Y_index,time_index,matrix(rpois(N*T*J,lambda = exp(err_temporal+err_alpha_individual)),ncol=J))
  return(list(group_data=Y12_grouped,individual_data=Y12_individual))
  }
############################################################################################################
############################################################################################################
#### Function for creating inla model version one, without collapsing alpha latent effect.
############################################################################################################
############################################################################################################

inla_data_frame_restructure=function(Selected,
                                     coeffs_pos = "b0_pos",
                                     coeffs_neg = "b0_neg",
                                     with_diurnal= TRUE){
  # coeffs_pos and # coeffs_neg are 
  # corresponding selected fixed effects
  
  
  Selected$logret_neg =  as.numeric(Selected$logret_neg)
  Selected$logret_pos =  as.numeric(Selected$logret_pos)
  
  Selected$dura_max_neg =  scale(as.numeric(Selected$dura_max_neg),center=TRUE,scale = TRUE)
  Selected$dura_max_pos =  scale(as.numeric(Selected$dura_max_pos),center = TRUE,scale=TRUE)
  
  Selected$dura_std_neg =  scale(as.numeric(Selected$dura_std_neg),center = TRUE,scale = TRUE)
  Selected$dura_std_pos =  scale(as.numeric(Selected$dura_std_pos),center = TRUE,scale=TRUE)
  
  Selected$logret_std_neg =  scale(as.numeric(Selected$logret_std_neg),center = TRUE,scale=TRUE)
  Selected$logret_std_pos =  scale(as.numeric(Selected$logret_std_pos),center = TRUE,scale=TRUE)
  
  Selected$time_neg =  as.numeric(Selected$time_neg)
  Selected$time_pos =  as.numeric(Selected$time_pos)
  
  stock.sample =unique(Selected$ROOT)
  
  S_m = length(unique(Selected$ROOT))
  Selected$s_index=factor(Selected$ROOT)
  levels(Selected$s_index) = 1:S_m
  S_m = length(unique(Selected$ROOT))
  T=length(unique(Selected$time_interval))
  N=nrow(Selected)
  s_index_temp=as.numeric(sort(unique(Selected$s_index)))
  
  
  Y_mat_pos=matrix(NA,nrow=T*S_m,ncol=S_m)
  Y_mat_pos_vec=matrix(NA,nrow=T*S_m,ncol=1)
  X_mat_pos=matrix(NA,nrow=T*S_m,ncol=S_m*length(coeffs_pos))
  colnames(X_mat_pos)=paste(coeffs_pos,rep(1:S_m,each=length(coeffs_pos)),sep='')
  for(s in seq_along(s_index_temp)){
    Y_temp= Selected%>%
      filter((response==1)&(s_index==s_index_temp[s]))%>%
      arrange((time_interval))
    Y_mat_pos[1:T+T*(s-1),s]=Y_temp$counts
    Y_mat_pos_vec[1:T+T*(s-1),1]=Y_temp$counts
    X_mat_pos[1:T+T*(s-1),1:length(coeffs_pos)+length(coeffs_pos)*(s-1)]=cbind(as.matrix(Y_temp[,coeffs_pos[-length(coeffs_pos)]]),rep(1,T))
  }
  
  
  
  Y_mat_neg=matrix(NA,nrow=T*S_m,ncol=S_m)
  Y_mat_neg_vec=matrix(NA,nrow=T*S_m,ncol=1)
  X_mat_neg=matrix(NA,nrow=T*S_m,ncol=S_m*length(coeffs_neg))
  colnames(X_mat_neg)=paste(coeffs_neg,rep(1:S_m,each=length(coeffs_neg)),sep='')
  for(s in seq_along(s_index_temp)){
    Y_temp= Selected%>%
      filter((response==2)&(s_index==s_index_temp[s]))%>%
      arrange((time_interval))
    Y_mat_neg[1:T+T*(s-1),s]=Y_temp$counts
    Y_mat_neg_vec[1:T+T*(s-1),1]=Y_temp$counts
    X_mat_neg[1:T+T*(s-1),1:length(coeffs_neg)+length(coeffs_neg)*(s-1)]=cbind(as.matrix(Y_temp[,coeffs_neg[-length(coeffs_neg)]]),rep(1,T))
  }
  
  Y_mat=cbind(rbind(Y_mat_pos,matrix(NA,nrow = nrow(Y_mat_neg),ncol=ncol(Y_mat_neg))),
              rbind(matrix(NA,nrow = nrow(Y_mat_pos),ncol=ncol(Y_mat_pos)),Y_mat_neg))
  Y_mat_vec=cbind(rbind(Y_mat_pos_vec,matrix(NA,nrow = nrow(Y_mat_neg_vec),ncol=ncol(Y_mat_neg_vec))),
                  rbind(matrix(NA,nrow = nrow(Y_mat_pos_vec),ncol=ncol(Y_mat_pos_vec)),Y_mat_neg_vec))
  # View(Y_mat_vec)
  X_mat=cbind(rbind(X_mat_pos,matrix(NA,nrow = nrow(X_mat_neg),ncol=ncol(X_mat_neg))),
              rbind(matrix(NA,nrow = nrow(X_mat_pos),ncol=ncol(X_mat_pos)),X_mat_neg))
  Selected_new=data.frame(X_mat)
  Selected_new$logret_std_neg=Selected$logret_std_neg
  Selected_new$logret_std_pos=Selected$logret_std_pos
  Selected_new$dura_std_pos=Selected$dura_std_pos
  Selected_new$dura_std_neg=Selected$dura_std_neg
  Selected_new$dura_max_pos=Selected$dura_max_pos
  Selected_new$dura_max_neg=Selected$dura_max_neg
  Selected_new$b0_neg=rep(c(NA,1),each=N/2)
  Selected_new$b0_pos=rep(c(1,NA),each=N/2)
  
  Selected_new$response=rep(1:2,each=nrow(Y_mat)/2)
  Selected_new$time_pos=NA;Selected$time_neg=NA
  Selected_new$time_pos[Selected$response==1]=rep(1:T,length(unique(Selected$s_index)))
  Selected_new$time_neg[Selected$response==2]=rep(1:T,length(unique(Selected$s_index)))
  Selected_new$time_interval= rep(1:T,2*length(unique(Selected$s_index)))
  Selected_new$BiAR1_time_idx = c(rep(1:T,length(unique(Selected$s_index))),
                                  rep((T+1):(2*T),length(unique(Selected$s_index))))
  Selected_new$alpha=1:nrow(Y_mat)# alpha_{j,st}
  # Selected_new$alpha2=rep(1:(2*T),each=length(s_index_temp))
  Selected_new$alpha2=c(rep(1:T,length(s_index_temp)),rep(1:T+T,length(s_index_temp)))
  
  alpha_s=matrix(NA,nrow=N,ncol=S_m)
  for(S_mi in 1:S_m){
    alpha_s[c(1:T+(S_mi-1)*T,(1:T)+(S_m+S_mi-1)*T),S_mi]=1:(2*T)
  }
  colnames(alpha_s)=paste("alpha_",1:S_m,sep='')
  
  Selected_new=cbind(Selected_new,alpha_s)
  s_index_temp_vec= rep(rep(s_index_temp,each = T),2)
  Selected_new$ROOT = s_index_temp_vec
  
  return(list(Y_mat_vec,Selected_new))
}

running_inla_version=function(Selected,Selected_new,version=1,pr=0.5){ # Selected is the name of the original data frame.
  S_m = length(unique(Selected$ROOT))
  T=length(unique(Selected$time_interval))
  N=nrow(Selected)
  s_index_temp=as.numeric(sort(unique(Selected$s_index)))

  if(version==1){
    ## This is model 3.2 ## Stock-specific
    formula_new_str=paste('Y_mat_vec','~-1+',paste(
      c('dura_max_neg','dura_std_neg','logret_std_neg','b0_neg','dura_max_pos','dura_std_pos','logret_std_pos','b0_pos'),rep(s_index_temp,each=8),sep='',collapse = '+'),
      "+f(time_pos,model='rw1')+f(time_neg,model='rw1')+f(alpha,model='iid2d',n=N,hyper=list(theta1=list(param=c(4,",pr,",",pr,",0))))")
    formula_new = as.formula(formula_new_str)
    # formula_new
    res= inla(formula_new,family =rep("poisson",ncol(Y_mat_vec)),
              data = Selected_new,
              control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),
              control.fixed=list(prec=list(default=1)),
              control.predictor=list(compute=TRUE,link = 1),
              num.thread=1)
  }
  else if(version==2){
    ## This is model 3.1 and 3.3 ## individual and common coef
    formula_new_str=paste('Y_mat_vec','~-1+',paste(
      c('dura_max_neg','dura_std_neg','logret_std_neg','b0_neg','dura_max_pos','dura_std_pos','logret_std_pos','b0_pos'),sep='',collapse = '+'),
      "+f(time_pos,model='rw1')+f(time_neg,model='rw1')+f(alpha,model='iid2d',n=N,hyper=list(theta1=list(param=c(4,",pr,",",pr,",0))))")
    formula_new = as.formula(formula_new_str)
    res= inla(formula_new,family =rep("poisson",ncol(Y_mat_vec)),
              data = Selected_new,
              control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),
              control.fixed=list(prec=list(default=1)),
              control.predictor=list(compute=TRUE,link = 1),
              num.thread=1)
  }
  # else if(version=='3'){
  #   formula_new_str=paste('Y_mat_vec','~-1+',paste(
  #     c('dura_max_neg','dura_std_neg','logret_std_neg','b0_neg','dura_max_pos','dura_std_pos','logret_std_pos','b0_pos'),rep(s_index_temp,each=8),sep='',collapse = '+'),
  #     "+f(time_pos,model='rw1')+f(time_neg,model='rw1')+f(alpha,model='iid2d',n=N,hyper=list(theta1=list(param=c(4,",pr,",",pr,",0))))")
  #   formula_new = as.formula(formula_new_str)
  #   # formula_new
  #   res= inla(formula_new,family =rep("poisson",ncol(Y_mat_vec)),
  #             data = Selected_new,
  #             control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),
  #             control.fixed=list(prec=list(default=1)),
  #             control.predictor=list(compute=TRUE,link = 1),
  #             num.thread=1)
  # else if(version=="1a"){
  #   formula_new_str=paste('Y_mat_vec','~-1+',paste(
  #     c('dura_max_neg','dura_std_neg','logret_std_neg','b0_neg','dura_max_pos','dura_std_pos','logret_std_pos','b0_pos'),rep(s_index_temp,each=8),sep='',collapse = '+'),
  #     "+f(time_pos,model='rw1')+f(time_neg,model='rw1')+",paste("f(alpha_",s_index_temp,",model='iid2d',n=2*T,hyper=list(theta1=list(param=c(4,",pr,",",pr,",","0))))",sep='',collapse = "+"))
  #   formula_new = as.formula(formula_new_str)
  #   # formula_new
  #   res= inla(formula_new,family =rep("poisson",ncol(Y_mat_vec)),
  #             data = Selected_new,
  #             control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),
  #             control.fixed=list(prec=list(default=1)),
  #             control.predictor=list(compute=TRUE,link = 1),
  #             num.thread=1)
  # }

  # else if(version=="2a"){
  #   formula_new_str=paste('Y_mat_vec','~-1+',paste(
  #     c('dura_max_neg','dura_std_neg','logret_std_neg','b0_neg','dura_max_pos','dura_std_pos','logret_std_pos','b0_pos'),sep='',collapse = '+'),
  #     "+f(time_pos,model='rw1')+f(time_neg,model='rw1')+",paste("f(alpha_",s_index_temp,",model='iid2d',n=2*T,hyper=list(theta1=list(param=c(4,",pr,",",pr,",","0))))",sep='',collapse = "+"))
  #   formula_new = as.formula(formula_new_str)
  #   # formula_new
  #   res= inla(formula_new,family =rep("poisson",ncol(Y_mat_vec)),
  #             data = Selected_new,
  #             control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE,config=TRUE),
  #             control.fixed=list(prec=list(default=1)),
  #             control.predictor=list(compute=TRUE,link = 1),
  #             num.thread=1)
  # }
  
  else{
   stopifnot(is.element(version,c(1:2)))
  }

  return(list(inla_res=(res),data_restr=Selected_new))
}


ggplot_for_coef=function(inla_obj,coeffs,n_sample,S_m,filefolder,mon){
  for(coeff_i in coeffs){
    coef_mat=rbind(rep(1,S_m))
    coef_list=split(coef_mat,col(coef_mat))
    names(coef_list)=paste(coeff_i,1:S_m,sep='')
    
    sample_temp=inla.posterior.sample(inla_list_temp$inla_res,n=n_sample,selection =coef_list)
    sample_mat=c()
    for(i in 1:n_sample){
      sample_mat=rbind(sample_mat,(sample_temp[[i]]$latent))
    }
    
    sample_data_frame=data.frame(Coef=factor(rep(stock.sample,n_sample),levels = stock.sample),samples=sample_mat[,1])
    # View(sample_data_frame)
    sample_data_frame%>%
      ggplot()+geom_boxplot(aes(Coef,samples,group=Coef,color=Coef),width=0.1)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      ggtitle(paste(coeff_i,' of different stocks'))+theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste(filefolder,"\\",coeff_i,'_',mon,'.pdf',sep=''))
  }
}
