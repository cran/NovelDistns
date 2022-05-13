# sample function for simulation

reggap <- function(n, dist, param)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!= "weibull" & dist!= "lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(n %% 1 !=0){stop("Parameter n must be an integer")}
  if(dist== "exp")
  {
    quant <- function(par,u)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      lambda = par[4]
      A <- ((1-u^(1/b))^(1/a))
      B <- -1+A
      C <- (B*log(alpha))/alpha
      D <- gsl::lambert_W0(C)
      E <- log(alpha) + D
      G <- log(alpha)/E
      values <- (1/lambda)*log(G)
    }
  }
  if(dist=="rayleigh")
  {
    quant <- function(par,u)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      delta = par[4]
      A <- log(alpha)^2
      B <- log(alpha)
      C <- (-1+(1-u^(1/b))^(1/a))*log(alpha)
      D <- gsl::lambert_Wm1(C/alpha)
      E <- (B + D)^2
      G <- delta * sqrt(log((A/E)))
    }
  }
  if(dist =="weibull")
  {
    quant <- function(par,u)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      theta = par[4]
      beta = par[5]
      A <- log(alpha)
      B <- (-1+(1-u^(1/b))^(1/a))*log(alpha)
      D <- A + gsl::lambert_Wm1(B/alpha)
      E <- (log(A/D))^(1/beta)
    }
  }
  if(dist =="lomax")
  {
    quant <- function(par,u)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      theta = par[4]
      beta = par[5]
      A <- log(alpha)
      B <- (-1+(1-u^(1/b))^(1/a))*log(alpha)
      D <- -1+(A/(A + gsl::lambert_Wm1(B/alpha)))^(1/beta)
      E <- (1/theta)*D
    }
  }
  return(quant(param, stats::runif(n)))
}


# quantile function

qeggap <- function(p, dist, param, log.p = FALSE, lower.tail = TRUE)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(p < 0 & p > 1){stop("Parameter p must be between 0 and 1")}
  if(dist== "exp")
  {
    quant <- function(par,p)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      lambda = par[4]
      A <- ((1-p^(1/b))^(1/a))
      B <- -1+A
      C <- (B*log(alpha))/alpha
      D <- gsl::lambert_W0(C)
      E <- log(alpha) + D
      G <- log(alpha)/E
      values <- (1/lambda)*log(G)
    }
  }
  if(dist=="rayleigh")
  {
    quant <- function(par,p)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      delta = par[4]
      A <- log(alpha)^2
      B <- log(alpha)
      C <- (-1+(1-p^(1/b))^(1/a))*log(alpha)
      D <- gsl::lambert_Wm1(C/alpha)
      E <- (B + D)^2
      G <- delta * sqrt(log((A/E)))
    }
  }

  if(dist =="weibull")
  {
    quant <- function(par,p)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      theta = par[4]
      beta = par[5]
      A <- log(alpha)
      B <- (-1+(1-p^(1/b))^(1/a))*log(alpha)
      D <- A + gsl::lambert_Wm1(B/alpha)
      E <- (log(A/D))^(1/beta)
    }
  }
  if(dist =="lomax")
  {
    quant <- function(par,p)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      theta = par[4]
      beta = par[5]
      A <- log(alpha)
      B <- (-1+(1-p^(1/b))^(1/a))*log(alpha)
      D <- -1+(A/(A + gsl::lambert_Wm1(B/alpha)))^(1/beta)
      E <- (1/theta)*D
    }
  }
  output <- quant(param,p)
  if(log.p==TRUE & lower.tail==FALSE) out <- quant(param,exp(p-1))
  if(log.p==TRUE & lower.tail==TRUE) out <- quant(param,exp(-p))
  if(log.p==FALSE & lower.tail==FALSE) out <- quant(param,exp(1-p))
  return(output)
}


# CDF

peggap <- function(data, dist, param, log.p = FALSE, lower.tail = TRUE)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    den = function(par,x){rate = par[4]; stats::dexp(x,rate)}
    cum = function(par,x){rate = par[4]; stats::pexp(x,rate)}
  }
  if(dist == "rayleigh")
  {
    den = function(par, x) {delta= par[4];(x/(delta^2))*(exp((-x^2)/(2*(delta^2))))}
    cum = function(par,x) {delta = par[4]; 1- exp((-x^2)/(2*(delta^2)))}
  }
  if(dist=="weibull")
  {
    den=function(par,x){shape=par[4]; scale=par[5]; stats::dweibull(x,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5];  stats::pweibull(x,shape,scale)}
  }
  if(dist=="lomax"){
    den=function(par,x){theta=par[4]; beta=par[5]; (theta*beta)/((1+theta*(x))^(beta+1))}
    cum=function(par,x){theta=par[4]; beta=par[5]; 1-(1+theta*(x))^(-beta)}
  }


    cdffinal <- function(par,x)
  {
    s0 = den(par,x)
    d0 = cum(par,x)
    a = par[1]
    b = par[2]
    alpha = par[3]
    (1-(1-(alpha*d0)/(alpha^d0))^a)^b
  }
  cdfoutput <- cdffinal(param,data)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdfoutput)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdfoutput)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdfoutput
  return(cdfoutput)
}


# density function

deggap <- function(data, dist, param, log = FALSE)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!="weibull" & dist!="lomax" )
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    den = function(par,x){rate = par[4]; stats::dexp(x,rate)}
    cum = function(par,x){rate = par[4]; stats::pexp(x,rate)}
  }
  if(dist == "rayleigh")
  {
    den = function(par, x) {delta= par[4];(x/(delta^2))*(exp((-x^2)/(2*(delta^2))))}
    cum = function(par,x) {delta = par[4]; 1- exp((-x^2)/(2*(delta^2)))}
  }

  if(dist=="weibull")
  {
    den=function(par,x){shape=par[4]; scale=par[5]; stats::dweibull(x,shape,scale)}
    cum=function(par,x){shape=par[4]; scale=par[5];  stats::pweibull(x,shape,scale)}
  }
  if(dist=="lomax"){
    den=function(par,x){theta=par[4]; beta=par[5]; (theta*beta)/((1+theta*(x))^(beta+1))}
    cum=function(par,x){theta=par[4]; beta=par[5]; 1-(1+theta*(x))^(-beta)}
  }

  pdffinal <- function(par,x)
  {
    s0 = den(par,x)
    d0 = cum(par,x)
    a = par[1]
    b = par[2]
    alpha = par[3]
    (a*b*s0*alpha^d0)*(log(alpha)*d0 + 1)*((1-(alpha*d0)/(alpha^d0))^(a-1))*((1-(1-(alpha*d0)/(alpha^d0))^a)^(b-1))
  }
  pdfoutput <-  pdffinal(param,data)
  if(log==TRUE){pdf<-log(pdfoutput)}
  return(pdfoutput)
}


# estimation

mleggap <- function(data, dist, starts, method ="SANN")
{
  if(dist!="exp"& dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    PDF0 <- function(par,x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      lambda = par[4]
      (((a*b*alpha^(exp(-lambda*x)))*(lambda*exp(-lambda*x))) - (a*b*alpha^(exp(-lambda*x))) * log(alpha)*lambda*(1-exp(-lambda*x))*exp(-lambda*x))*(((1-(1-(alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x))))^a)^(b-1))*(1-(alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x))))^(a-1))
    }
    CDF0 <- function(par,x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      lambda = par[4]
      (1-(1-(alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x))))^a)^b
    }
  }

  if(dist=="weibull"){

    PDF0 <- function(par,x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      lambda = par[4]
      beta = par[5]
      A <- a*b*lambda*beta*(alpha^(exp(-lambda*x^beta)))*(exp(-lambda*(x^beta)))*(1-(log(alpha)*(1-exp(-lambda*(x^beta)))))
      B <- (1-((alpha*(1-exp(-lambda*(x^beta))))/(alpha^(1-exp(-lambda*(x^beta))))))^(a-1)
      C <- (1-B)^(b-1)
      pdfegw <- A*B*C
    }

    CDF0 <- function(par,x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      lambda = par[4]
      beta = par[5]
      cdfegw <- (1-(1-((alpha*(1-exp(-(lambda)*x^(beta))))/(alpha^(1-exp(-((lambda)*x^(beta)))))))^a)^b
    }

  }

  if(dist=="rayleigh"){
    PDF0 <- function(par, x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      delta = par[4]
      A <- 1- exp((-x^2)/(2*(delta^2)))
      B <- (x/(delta^2))*(exp((-x^2)/(2*(delta^2))))
      C <- (a*b*B*(alpha^(1-A))) - (A*B*a*b*log(alpha)*(alpha^(1-A)))
      D <- (1-((alpha*A)/(alpha^A)))^(a-1)
      E <- (1-((1-((alpha*A)/(alpha^A)))^(a)))^(b-1)
      G <- C*D*E
    }
    CDF0 <- function(par,x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      delta = par[4]
      A <- 1- exp((-x^2)/(2*(delta^2)))
      B <- (alpha*A)/(alpha^A)
      C <- (1-(1-B)^a)^b
    }
  }

  if(dist=="lomax"){
    PDF0 <- function(par, x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      beta = par[4]
      theta = par[5]
      A <- (theta*beta)/(1+theta*x)^(beta+1)
      B <- 1-(1+theta*x)^(-beta)
      C <- (a*b*A*(alpha^(1-B))) - (A*B*a*b*log(alpha)*(alpha^(1-B)))
      D <- (1-((alpha*B)/(alpha^B)))^(a-1)
      E <- (1-((1-((alpha*B)/(alpha^B)))^(a)))^(b-1)
      G <- C*D*E
    }
    CDF0 <- function(par,x)
    {
      a = par[1]
      b = par[2]
      alpha = par[3]
      beta = par[4]
      theta = par[5]
      A <- 1-(1+theta*x)^(-beta)
      B <- (alpha*A)/(alpha^A)
      C <- (1-(1-B)^a)^b
    }
  }


  ans<-suppressWarnings(AdequacyModel::goodness.fit(pdf=PDF0,
                                                    cdf=CDF0,
                                                    starts = starts, data = data,
                                                    method= method,mle=NULL))

  aux=cbind(ans$mle,ans$Erro,ans$mle+stats::qnorm(0.025)*ans$Erro,ans$mle+stats::qnorm(0.975)*ans$Erro)

  colnames(aux)=c("MLE","Std. Error.","Lower. 95% CI","Upper 95% CI")

  aux1=cbind(ans$AIC, ans$CAIC, ans$BIC, ans$HQIC, ans$W, ans$A, ans$Value)
  colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "-log-Likelihood")
  rownames(aux1)=c("")

  aux2=cbind(ans$KS$statistic,ans$KS$p.value)
  colnames(aux2)=c("KS Statistic","KS p-value")
  rownames(aux2)=c("")

  aux3=cbind(if(ans$Convergence==0){"Algorithm Converged"} else{"Algorithm Not Converged"})
  colnames(aux3)=c("")
  rownames(aux3)=c("")

  list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}




regap <- function(n, dist, param)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!= "weibull" & dist!= "lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(n %% 1 !=0){stop("Parameter n must be an integer")}
  if(dist== "exp")
  {
    quant <- function(par,u)
    {
      a = par[1]
      alpha = par[2]
      lambda = par[3]
      A <- (-log(alpha)*exp(log(u)/alpha))/a
      B <- log(alpha)+gsl::lambert_Wm1(A)
      C <- log(B/(log(alpha)))
      values <- (-1/lambda)*C
    }
  }
  if(dist =="rayleigh")
  {
    y <- NA
    quant <- function(par,u)
    {
      a = par[1]
      alpha = par[2]
      delta = par[3]
      cdf_R <- function(par,x)
      {
        a =par[1]
        alpha = par[2]
        delta = par[3]
        ((alpha*(1- exp((-x^2)/(2*(delta^2)))))/(alpha^(1- exp((-x^2)/(2*(delta^2))))))^a
      }
      for(m in 1:n)
      {
        f=function(x){
          cdf_R(par,x)-u[m]
        }
        y[m] = min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
      }
      return(y)
    }
  }
  if(dist =="lomax")
  {
    y <- NA
    quant <- function(par,u)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      cdf_R <- function(par,x)
      {
        a =par[1]
        alpha = par[2]
        theta = par[3]
        beta = par[4]
        ((alpha*(1-(1+theta*x)^(-beta)))/(alpha^(1-(1+theta*x)^(-beta))))^a
      }
      for(m in 1:n)
      {
        f=function(x){
          cdf_R(par,x)-u[m]
        }
        y[m] = min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
      }
      return(y)
    }
  }
  if(dist =="weibull")
  {
    quant <- function(par,u)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      A <- (-log(alpha)*(u^(1/a)))/alpha
      B <- log(alpha)+gsl::lambert_Wm1(A)
      C <- -log(B/log(alpha))
      D <- theta*C^(1/beta)
    }
  }
  return(quant(param, stats::runif(n)))
}


qegap <- function(p, dist, param, log.p = FALSE, lower.tail = TRUE)
{
  if(dist!="exp" & dist!="rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(p < 0 & p > 1){stop("Parameter p must be between 0 and 1")}
  if(dist== "exp")
  {
    quant <- function(par,p)
    {
      a = par[1]
      alpha = par[2]
      lambda = par[3]
      A <- (-log(alpha)*exp(log(p)/alpha))/a
      B <- log(alpha)+gsl::lambert_Wm1(A)
      C <- log(B/(log(alpha)))
      values <- (-1/lambda)*C
    }
  }

  if(dist =="rayleigh")
  {
    quant <- function(par,p)
    {
      a = par[1]
      alpha = par[2]
      delta = par[3]
      cdf_R <- function(par,x)
      {
        a =par[1]
        alpha = par[2]
        delta = par[3]
        ((alpha*(1- exp((-x^2)/(2*(delta^2)))))/(alpha^(1- exp((-x^2)/(2*(delta^2))))))^a
      }
      f=function(x){
        cdf_R(par,x)-p
      }
      min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
    }
  }
  if(dist =="weibull")
  {
    quant <- function(par,p)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      A <- (-log(alpha)*(p^(1/a)))/alpha
      B <- log(alpha)+gsl::lambert_Wm1(A)
      C <- -log(B/log(alpha))
      D <- theta*C^(1/beta)
    }
  }

  if(dist =="lomax")
  {
    quant <- function(par,p)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      cdf_R <- function(par,x)
      {
        a = par[1]
        alpha = par[2]
        theta = par[3]
        beta = par[4]
        ((alpha*(1-(1+theta*x)^(-beta)))/(alpha^(1-(1+theta*x)^(-beta))))^a
      }
      f=function(x){
        cdf_R(par,x)-p
      }
      min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
    }
  }
  output <- quant(param,p)
  if(log.p==TRUE & lower.tail==FALSE) out <- quant(param,exp(p-1))
  if(log.p==TRUE & lower.tail==TRUE) out <- quant(param,exp(-p))
  if(log.p==FALSE & lower.tail==FALSE) out <- quant(param,exp(1-p))
  return(output)
}



pegap <- function(data, dist, param, log.p = FALSE, lower.tail = TRUE)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    den = function(par,x){rate = par[3]; stats::dexp(x,rate)}
    cum = function(par,x){rate = par[3]; stats::pexp(x,rate)}
  }
  if(dist == "rayleigh")
  {
    den = function(par, x) {delta= par[3];(x/(delta^2))*(exp((-x^2)/(2*(delta^2))))}
    cum = function(par,x) {delta = par[3]; 1- exp((-x^2)/(2*(delta^2)))}
  }
  if(dist=="weibull")
  {
    den=function(par,x){shape=par[3]; scale=par[4]; stats::dweibull(x,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4];  stats::pweibull(x,shape,scale)}
  }
  if(dist=="lomax"){
    den=function(par,x){theta=par[3]; beta=par[4]; (theta*beta)/((1+theta*(x))^(beta+1))}
    cum=function(par,x){theta=par[3]; beta=par[4]; 1-(1+theta*(x))^(-beta)}
  }
  cdffinal <- function(par,x)
  {
    s0 = den(par,x)
    d0 = cum(par,x)
    a = par[1]
    alpha = par[2]
    ((alpha*d0)/(alpha^d0))^a
  }
  cdfoutput <- cdffinal(param,data)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdfoutput)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdfoutput)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdfoutput
  return(cdfoutput)
}



degap <- function(data, dist, param, log = FALSE)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    den = function(par,x){rate = par[3]; stats::dexp(x,rate)}
    cum = function(par,x){rate = par[3]; stats::pexp(x,rate)}
  }
  if(dist == "rayleigh")
  {
    den = function(par, x) {delta= par[3];(x/(delta^2))*(exp((-x^2)/(2*(delta^2))))}
    cum = function(par,x) {delta = par[3]; 1- exp((-x^2)/(2*(delta^2)))}
  }

  if(dist=="weibull")
  {
    den=function(par,x){shape=par[3]; scale=par[4]; stats::dweibull(x,shape,scale)}
    cum=function(par,x){shape=par[3]; scale=par[4];  stats::pweibull(x,shape,scale)}
  }

  if(dist=="lomax"){
    den=function(par,x){theta=par[3]; beta=par[4]; (theta*beta)/((1+theta*(x))^(beta+1))}
    cum=function(par,x){theta=par[3]; beta=par[4]; 1-(1+theta*(x))^(-beta)}
  }

  pdffinal <- function(par,x)
  {
    s0 = den(par,x)
    d0 = cum(par,x)
    a = par[1]
    alpha = par[2]
    a*(((alpha*d0)/(alpha^d0))^(a-1))*s0*((alpha^(1-d0))-(log(alpha)*(alpha^(1-d0))))
  }
  pdfoutput <-  pdffinal(param,data)
  if(log==TRUE){pdf<-log(pdfoutput)}
  return(pdfoutput)
}



mlegap <- function(data, dist, starts, method ="SANN")
{
  if(dist!="exp"& dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    PDF0 <- function(par,x)
    {
      a = par[1]
      alpha = par[2]
      lambda = par[3]
      (((a*alpha^(exp(-lambda*x)))*(lambda*exp(-lambda*x))) - (a*alpha^(exp(-lambda*x))) * log(alpha)*lambda*(1-exp(-lambda*x))*exp(-lambda*x))*(((1-(1-(alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x)))))^(a-1)))
    }
    CDF0 <- function(par,x)
    {
      a = par[1]
      alpha = par[2]
      lambda = par[3]
      (1-(1-(alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x)))))^a
    }
  }

  if(dist=="weibull"){

    PDF0 <- function(par,x)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      A <- 1- exp(-(x/theta)^beta)
      B <- (theta/beta)*((x/beta)^(theta-1))*exp(-(x/beta)^theta)
      C <- (alpha*A)/(alpha^A)
      D <- a*(C^(a-1))*B*((alpha^(1-A))-log(alpha)*(alpha^(1-A)))
    }

    CDF0 <- function(par,x)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      A <- 1- exp(-(x/theta)^beta)
      B <- ((alpha*A)/(alpha^A))^a
    }

  }

  if(dist=="rayleigh"){
    PDF0 <- function(par, x)
    {
      a = par[1]
      alpha = par[2]
      delta = par[3]
      A <- 1- exp((-x^2)/(2*(delta^2)))
      B <- (x/(delta^2))*(exp((-x^2)/(2*(delta^2))))
      C <- (a*B*(alpha^(1-A))) - (A*B*a*log(alpha)*(alpha^(1-A)))
      E <- ((alpha*A)/(alpha^A))^(a-1)
      G <- C*E
    }
    CDF0 <- function(par,x)
    {
      a = par[1]
      alpha = par[2]
      delta = par[3]
      A <- 1- exp((-x^2)/(2*(delta^2)))
      B <- (alpha*A)/(alpha^A)
      C <- B^a
    }
  }
  if(dist=="lomax"){
    PDF0 <- function(par, x)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      A <- 1- (1+theta*x)^(-beta)
      B <- (theta*beta)*(1+theta*x)^(-beta+1)
      C <- (a*B*(alpha^(1-A))) - (A*B*a*log(alpha)*(alpha^(1-A)))
      E <- ((alpha*A)/(alpha^A))^(a-1)
      G <- C*E
    }
    CDF0 <- function(par,x)
    {
      a = par[1]
      alpha = par[2]
      theta = par[3]
      beta = par[4]
      A <- 1- (1+theta*x)^(-beta)
      B <- (alpha*A)/(alpha^A)
      C <- B^a
    }
  }

  ans<-suppressWarnings(AdequacyModel::goodness.fit(pdf=PDF0,
                                                    cdf=CDF0,
                                                    starts = starts, data = data,
                                                    method= method,mle=NULL))

  aux=cbind(ans$mle,ans$Erro,ans$mle+stats::qnorm(0.025)*ans$Erro,ans$mle+stats::qnorm(0.975)*ans$Erro)

  colnames(aux)=c("MLE","Std. Error.","Lower. 95% CI","Upper 95% CI")

  aux1=cbind(ans$AIC, ans$CAIC, ans$BIC, ans$HQIC, ans$W, ans$A, ans$Value)
  colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "-log-Likelihood")
  rownames(aux1)=c("")

  aux2=cbind(ans$KS$statistic,ans$KS$p.value)
  colnames(aux2)=c("KS Statistic","KS p-value")
  rownames(aux2)=c("")

  aux3=cbind(if(ans$Convergence==0){"Algorithm Converged"} else{"Algorithm Not Converged"})
  colnames(aux3)=c("")
  rownames(aux3)=c("")

  list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}



rgap <- function(n, dist, param)
{
  if(dist!="exp" & dist!="rayleigh" & dist!= "weibull" & dist!= "lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(n %% 1 !=0){stop("Parameter n must be an integer")}
  if(dist== "exp")
  {
    quant <- function(par,u)
    {
      alpha = par[1]
      lambda = par[2]
      A <- gsl::lambert_Wm1((-u*log(alpha))/alpha)
      B <- log((A+log(alpha))/log(alpha))
      values <- -(1/lambda)*B
    }
  }

  if(dist== "rayleigh")
  {
    y <- NA
    quant <- function(par,u)
    {
      cdf_R <- function(par,x)
      {

        alpha = par[1]
        delta = par[2]
        (alpha*(1- exp((-x^2)/(2*(delta^2)))))/(alpha^(1- exp((-x^2)/(2*(delta^2)))))
      }

      for(m in 1:n)
      {
        f=function(x){
          alpha = par[1]
          delta = par[2]
          cdf_R(par,x)-u[m]
        }
        y[m]=min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
      }
      return(y)
    }
  }

  if(dist =="lomax")
  {
    y <- NA
    quant <- function(par,u)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      cdf_R <- function(par,x)
      {
        alpha = par[1]
        theta = par[2]
        beta = par[3]
        (alpha*(1-(1+theta*x)^(-beta)))/(alpha^(1-(1+theta*x)^(-beta)))
      }
      for(m in 1:n)
      {
        f=function(x){
          cdf_R(par,x)-u[m]
        }
        y[m] = min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
      }
      return(y)
    }
  }

  if(dist =="weibull")
  {
    quant <- function(par,u)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      A <- gsl::lambert_Wm1((-u*log(alpha))/alpha)
      B <- -log((A+log(alpha))/log(alpha))
      C <- theta*(B^(1/beta))
    }
  }
  return(quant(param, stats::runif(n)))
}




qgap <- function(p, dist, param, log.p = FALSE, lower.tail = TRUE)
{
  if(dist!="exp" & dist!="rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(p < 0 & p > 1){stop("Parameter p must be between 0 and 1")}
  if(dist== "exp")
  {
    quant <- function(par,p)
    {
      alpha = par[1]
      lambda = par[2]
      A <- gsl::lambert_Wm1((-p*log(alpha))/alpha)
      B <- log((A+log(alpha))/log(alpha))
      values <- -(1/lambda)*B
    }
  }

  if(dist== "rayleigh")
  {
    quant <- function(par,p)
    {
      cdf_R <- function(par,x)
      {

        alpha = par[1]
        delta = par[2]
        (alpha*(1- exp((-x^2)/(2*(delta^2)))))/(alpha^(1- exp((-x^2)/(2*(delta^2)))))
      }

      alpha = par[1]
      delta = par[2]
      f=function(x){
        cdf_R(par,x)-p
      }
      x=min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
    }
  }
  if(dist =="weibull")
  {
    quant <- function(par,p)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      A <- gsl::lambert_Wm1((-p*log(alpha))/alpha)
      B <- -log((A+log(alpha))/log(alpha))
      C <- theta*(B^(1/beta))
    }
  }
  if(dist =="lomax")
  {
    quant <- function(par,p)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      cdf_R <- function(par,x)
      {
        alpha = par[1]
        theta = par[2]
        beta = par[3]
        (alpha*(1-(1+theta*x)^(-beta)))/(alpha^(1-(1+theta*x)^(-beta)))
      }
      f=function(x){
        cdf_R(par,x)-p
      }
      min(rootSolve::uniroot.all(f,lower=0,upper=100000000000,tol=0.0001))
    }
  }
  output <- quant(param,p)
  if(log.p==TRUE & lower.tail==FALSE) out <- quant(param,exp(p-1))
  if(log.p==TRUE & lower.tail==TRUE) out <- quant(param,exp(-p))
  if(log.p==FALSE & lower.tail==FALSE) out <- quant(param,exp(1-p))
  return(output)
}


pgap <- function(data, dist, param, log.p = FALSE, lower.tail = TRUE)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    den = function(par,x){rate = par[2]; stats::dexp(x,rate)}
    cum = function(par,x){rate = par[2]; stats::pexp(x,rate)}
  }
  if(dist == "rayleigh")
  {
    den = function(par, x) {delta= par[2];(x/(delta^2))*(exp((-x^2)/(2*(delta^2))))}
    cum = function(par,x) {delta = par[2]; 1- exp((-x^2)/(2*(delta^2)))}
  }
  if(dist=="weibull")
  {
    den=function(par,x){shape=par[2]; scale=par[3]; stats::dweibull(x,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3];  stats::pweibull(x,shape,scale)}
  }
  if(dist=="lomax"){
    den=function(par,x){theta=par[2]; beta=par[3]; (theta*beta)/((1+theta*(x))^(beta+1))}
    cum=function(par,x){theta=par[2]; beta=par[3]; 1-(1+theta*(x))^(-beta)}
  }
  cdffinal <- function(par,x)
  {
    s0 = den(par,x)
    d0 = cum(par,x)
    alpha = par[1]
    (alpha*d0)/(alpha^d0)
  }
  cdfoutput <- cdffinal(param,data)
  if(log.p==TRUE & lower.tail == FALSE) cdf<-log(1-cdfoutput)
  if(log.p==TRUE & lower.tail == TRUE) cdf<-log(cdfoutput)
  if(log.p==FALSE & lower.tail == FALSE) cdf<-1-cdfoutput
  return(cdfoutput)
}

dgap <- function(data, dist, param, log = FALSE)
{
  if(dist!="exp" & dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    den = function(par,x){rate = par[2]; stats::dexp(x,rate)}
    cum = function(par,x){rate = par[2]; stats::pexp(x,rate)}
  }
  if(dist == "rayleigh")
  {
    den = function(par, x) {delta= par[2];(x/(delta^2))*(exp((-x^2)/(2*(delta^2))))}
    cum = function(par,x) {delta = par[2]; 1- exp((-x^2)/(2*(delta^2)))}
  }

  if(dist=="weibull")
  {
    den=function(par,x){shape=par[2]; scale=par[3]; stats::dweibull(x,shape,scale)}
    cum=function(par,x){shape=par[2]; scale=par[3];  stats::pweibull(x,shape,scale)}
  }

  if(dist=="lomax"){
    den=function(par,x){theta=par[2]; beta=par[3]; (theta*beta)/((1+theta*(x))^(beta+1))}
    cum=function(par,x){theta=par[2]; beta=par[3]; 1-(1+theta*(x))^(-beta)}
  }

  pdffinal <- function(par,x)
  {
    s0 = den(par,x)
    d0 = cum(par,x)
    alpha = par[1]
    s0*(alpha^(1-d0))*(1-log(alpha)*d0)
  }
  pdfoutput <-  pdffinal(param,data)
  if(log==TRUE){pdf<-log(pdfoutput)}
  return(pdfoutput)
}


mlgap <- function(data, dist, starts, method ="SANN")
{
  if(dist!="exp"& dist!= "rayleigh" & dist!="weibull" & dist!="lomax")
  {stop("Baseline distribution not implemented or misspelled.")}
  if(dist == "exp")
  {
    PDF0 <- function(par,x)
    {
      alpha = par[1]
      lambda = par[2]
      (((alpha^(exp(-lambda*x)))*(lambda*exp(-lambda*x))) - (alpha^(exp(-lambda*x))) * log(alpha)*lambda*(1-exp(-lambda*x))*exp(-lambda*x))

    }
    CDF0 <- function(par,x)
    {
      alpha = par[1]
      lambda = par[2]
      (alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x)))
    }
  }

  if(dist=="weibull"){

    PDF0 <- function(par,x)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      A <- 1- exp(-(x/theta)^beta)
      B <- (theta/beta)*((x/beta)^(theta-1))*exp(-(x/beta)^theta)
      C <- (B*alpha^(1-A))*(1-log(alpha)*A)
    }

    CDF0 <- function(par,x)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      A <- 1- exp(-(x/theta)^beta)
      B <- (alpha*A)/(alpha^A)
    }

  }

  if(dist=="rayleigh"){
    PDF0 <- function(par, x)
    {
      alpha = par[1]
      delta = par[2]
      A <- 1- exp((-x^2)/(2*(delta^2)))
      B <- (x/(delta^2))*(exp((-x^2)/(2*(delta^2))))
      C <- (B*(alpha^(1-A))) - (A*B*log(alpha)*(alpha^(1-A)))
    }
    CDF0 <- function(par,x)
    {
      alpha = par[1]
      delta = par[2]
      A <- 1- exp((-x^2)/(2*(delta^2)))
      B <- (alpha*A)/(alpha^A)
    }
  }
  if(dist=="lomax"){
    PDF0 <- function(par, x)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      A <- 1- (1+theta*x)^(-beta)
      B <- (theta*beta)*(1+theta*x)^(-beta+1)
      C <- (B*(alpha^(1-A))) - (A*B*log(alpha)*(alpha^(1-A)))
      E <- (alpha*A)/(alpha^A)
      G <- C*E
    }
    CDF0 <- function(par,x)
    {
      alpha = par[1]
      theta = par[2]
      beta = par[3]
      A <- 1- (1+theta*x)^(-beta)
      B <- (alpha*A)/(alpha^A)
    }
  }


  ans<-suppressWarnings(AdequacyModel::goodness.fit(pdf=PDF0,
                                                    cdf=CDF0,
                                                    starts = starts, data = data,
                                                    method= method,mle=NULL))

  aux=cbind(ans$mle,ans$Erro,ans$mle+stats::qnorm(0.025)*ans$Erro,ans$mle+stats::qnorm(0.975)*ans$Erro)

  colnames(aux)=c("MLE","Std. Error.","Lower. 95% CI","Upper 95% CI")

  aux1=cbind(ans$AIC, ans$CAIC, ans$BIC, ans$HQIC, ans$W, ans$A, ans$Value)
  colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "-log-Likelihood")
  rownames(aux1)=c("")

  aux2=cbind(ans$KS$statistic,ans$KS$p.value)
  colnames(aux2)=c("KS Statistic","KS p-value")
  rownames(aux2)=c("")

  aux3=cbind(if(ans$Convergence==0){"Algorithm Converged"} else{"Algorithm Not Converged"})
  colnames(aux3)=c("")
  rownames(aux3)=c("")

  list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}




