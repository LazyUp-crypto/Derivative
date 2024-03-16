BSMPrice <- function(S0,K,T,r,sigma,Type){
  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1-sigma*sqrt(T)
  if(Type=="Call"){
    return(S0*pnorm(d1)-K*exp(-r*T)*pnorm(d2))
  }else{
    return(K*exp(-r*T)*pnorm(-d2)-S0*pnorm(-d1))
  }
}
BSMPrice(36,40,0.5,0.06,0.2,"Put")
