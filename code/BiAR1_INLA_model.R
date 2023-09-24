'inla.rgeneric.BiAR1.model' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  #Internal function
  interpret.theta <- function() {
    return(
      list(
        prec_1 = exp(theta[1L]),
           prec_2 = exp(theta[2L]),
        rho =(exp(theta[3L])-1) / (1 + exp(theta[3L])),
           phi_1 = (exp(theta[4L])-1) / (1 + exp(theta[4L])),
           phi_2 = (exp(theta[5L])-1) / (1 + exp(theta[5L]))
        )
      
    )
  }
  
  graph <- function(){
    require(Matrix)
    
    Lambda = sparseMatrix(
      i=c(1L:(2L*T),2L*(1L:T)-1L),
      j=c(1L:(2L*T),2L*(1L:T)),
      x=1,
      repr = "C",symmetric = TRUE
    )
    
    # Linear transformation matrix from
    # (x_1,x_2,...x_t, y_1,y_2,...y_t) to
    # (w_1,u_1,w_2,u_2,...,w_t,u_t), t>1
    A = sparseMatrix(
      i=c(1L,2L,rep(3L:(2L*T),each=2L)),
      j=c(1L,T+1L,T+1L+cumsum(rep(c(-T,1L,T-1L,1L),T-1L))),
      x =c(1)) 
    

Prec_mat_Graph = t(A)%*%Lambda%*%A 
    return(((Prec_mat_Graph!=0)*1))
  }
  
  Q <- function() {

    param0 <- interpret.theta()
    
    Lambda = sparseMatrix(
      i=c(1L:(2L*T),2L*(1L:T)-1L),
      j=c(1L:(2L*T),2L*(1L:T)),
      x=c(rep(1,2L*T),rep(-param0$rho,T)),
      repr = "C",symmetric = TRUE
    )
    
    # Linear transformation matrix from
    # (x_1,x_2,...x_t, y_1,y_2,...y_t) to
    # (w_1,u_1,w_2,u_2,...,w_t,u_t), t>1
    A = sparseMatrix(
      i=c(1L,2L,rep(3L:(2L*T),each=2L)),
      j=c(1L,T+1L,T+1L+cumsum(rep(c(-T,1L,T-1L,1L),T-1L))),
      x =c((1-param0$phi_1**2)**0.5*param0$prec_1**0.5,
           (1-param0$phi_2**2)**0.5*param0$prec_2**0.5,
           rep(c(-param0$phi_1*param0$prec_1**0.5,param0$prec_1**0.5,
                 -param0$phi_2*param0$prec_2**0.5,param0$prec_2**0.5),T-1L)) 
    )
    Prec_mat = t(A)%*%Lambda%*%A/(1-param0$rho^2) 

    return(Prec_mat)
  }
  
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  log.prior <- function() {
    require(CholWishart)
    param = interpret.theta()
    
    res <-dnorm(log(param$phi_1+1)-log(1-param$phi_1),0,.15)+
      dnorm(log(param$phi_2+1)-log(1-param$phi_2),0,.15)+
      dWishart((1-param$rho**2)/param$prec_1/param$prec_2*matrix(
        c(param$prec_2^-1,-param$rho*param$prec_1^-0.5*param$prec_2^-.5,
          -param$rho*param$prec_1^-0.5*param$prec_2^-.5,param$prec_1^-1),
        ncol=2
      ),df=7,Sigma = diag(rep(1,2)),log = TRUE)+log(param$prec_1)+
      log(param$prec_2)+log(0.5)+log(1+param$rho)+log(1-param$rho)

    
    return(res)
  }
  
  initial <- function() {
    return(c(4,4,4,0,0))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}


'inla.rgeneric.BiRW1.model' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {
  
  #Internal function
  interpret.theta <- function() {
    return(
      list(
        prec_1 = exp(theta[1L]),
        prec_2 = exp(theta[2L]),
        rho =(exp(theta[3L])-1) / (1 + exp(theta[3L]))
      )
      
    )
  }
  
  graph <- function(){
    require(Matrix)
    
    Lambda = sparseMatrix(
      i=c(1L:(2L*T),2L*(1L:T)-1L),
      j=c(1L:(2L*T),2L*(1L:T)),
      x=1,
      repr = "C",symmetric = TRUE
    )
    
    # Linear transformation matrix from
    # (x_1,x_2,...x_t, y_1,y_2,...y_t) to
    # (w_1,u_1,w_2,u_2,...,w_t,u_t), t>1
    A = sparseMatrix(
      i=c(1L,2L,rep(3L:(2L*T),each=2L)),
      j=c(1L,T+1L,T+1L+cumsum(rep(c(-T,1L,T-1L,1L),T-1L))),
      x =c(1)) 
    
    
    Prec_mat_Graph = t(A)%*%Lambda%*%A 
    return(((Prec_mat_Graph!=0)*1))
  }
  
  Q <- function() {
    
    param0 <- interpret.theta()
    
    Lambda = sparseMatrix(
      i=c(1L:(2L*T),2L*(1L:T)-1L),
      j=c(1L:(2L*T),2L*(1L:T)),
      x=c(rep(1,2L*T),rep(-param0$rho,T)),
      repr = "C",symmetric = TRUE
    )
    
    # Linear transformation matrix from
    # (x_1,x_2,...x_t, y_1,y_2,...y_t) to
    # (w_1,u_1,w_2,u_2,...,w_t,u_t), t>1
    A = sparseMatrix(
      i=c(1L,2L,rep(3L:(2L*T),each=2L)),
      j=c(1L,T+1L,T+1L+cumsum(rep(c(-T,1L,T-1L,1L),T-1L))),
      x =c(param0$prec_1**0.5,
           param0$prec_2**0.5,
           rep(c(-param0$prec_1**0.5,param0$prec_1**0.5,
                 -param0$prec_2**0.5,param0$prec_2**0.5),T-1L)) 
    )
    Prec_mat = t(A)%*%Lambda%*%A/(1-param0$rho^2) 
    
    return(Prec_mat)
  }
  
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
    
  }
  
  log.prior <- function() {
    require(CholWishart)
    param = interpret.theta()
    
    res <-dWishart((1-param$rho**2)/param$prec_1/param$prec_2*matrix(
      c(param$prec_2^-1,-param$rho*param$prec_1^-0.5*param$prec_2^-.5,
        -param$rho*param$prec_1^-0.5*param$prec_2^-.5,param$prec_1^-1),
      ncol=2
    ),df=7,Sigma = diag(rep(1,2)),log = TRUE)+log(param$prec_1)+
      log(param$prec_2)+log(0.5)+log(1+param$rho)+log(1-param$rho)
    
    
    return(res)
  }
  
  initial <- function() {
    return(c(4,4,4))
  }
  
  quit <- function() {
    return(invisible())
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}


