simulateSA<-function(
    N, 
    TT, 
    r, 
    tr.threshold = 0.5, 
    tr.start = 10, 
    p = 10, 
    beta = c(4,4,2,2,rep(0, 6)),    
    alpha = NULL, 
    xi = NULL,    
    mu=0, 
    force = 3, 
    factor = NULL, 
    Rtype = "n", 
    Ftype = "drift", 
    Fsize = 1, 
    eff.size = 1,
    eff.noise = 1,
    tr.noise = 0.3,
    time.invariant = TRUE, 
    ar.coef = 0.7,
    error.type = "n", 
    seed=NULL
) {
  if (!is_count(N, 1) || !is_count(TT, 2) || !is_count(r, 2) || !is_count(p))
    stop("Require N >= 1, TT >= 2, r >= 2 and p >= 0. Treatment assignment uses two loadings.")
  if (!is_count(force) || force > 3 || !is_switch(time.invariant))
    stop("force must be 0, 1, 2 or 3; time.invariant must be a logical flag.")
  if (!length(tr.threshold) || length(tr.threshold) != length(tr.start) ||
      any(!is.finite(tr.threshold) | tr.threshold <= 0 | tr.threshold >= 1) ||
      any(!is.finite(tr.start) | tr.start != floor(tr.start) | tr.start < 1 | tr.start > TT))
    stop("Supply one valid adoption period per probability threshold.")
  if (length(beta) < p || any(!is.finite(beta[seq_len(p)])))
    stop("Supply a finite fixed coefficient for every covariate.")
  if (!is.null(alpha) && !identical(dim(alpha), as.integer(c(N, p)))) stop("alpha must be N by p.")
  if (!is.null(xi) && !identical(dim(xi), as.integer(c(TT, p)))) stop("xi must be TT by p.")
  if (!is.null(factor) && !identical(dim(factor), as.integer(c(TT, r)))) stop("factor must be TT by r.")
  if (!length(Ftype) || !length(Fsize) || !(length(Ftype) %in% c(1, r)) ||
      !(length(Fsize) %in% c(1, r)) || any(!is.finite(Fsize)))
    stop("Ftype and Fsize must have length one or r; Fsize must be finite.")
  error.type <- match.arg(error.type, c("n", "g", "t"))
  
  if (is.null(seed)==FALSE) {set.seed(seed)}
  
  
  getTS <- function(type = "ar1", TT, ar.coef) {
    if (type %in% c("ar1","drift")) {
      ts <- arima.sim(list(order = c(1,0,0), ar = ar.coef), n = TT)               
    } else {
      ts <- rnorm(TT)
    }
    if (type == "drift") {
      ts <- ts + seq(from=0, length.out = TT, by=1)
    }
    return(ts)
  }
  ntr <- length(tr.threshold) + 1
  if (r > 0) {
    if (Rtype %in% c("u","uniform")) {
      bound <- sqrt(3)
      lambda<-matrix(runif(N*r, min = 0, max = 2*bound), N, r, byrow=TRUE)
    } else {
      lambda<-matrix(rnorm(N*r, 1, 1), N, r, byrow=TRUE) # mean = 1
    }      
    
    if (is.null(factor) == TRUE) {
      factor<-matrix(NA,TT,r)   
      if (length(Ftype)==1) {Ftype <- rep(Ftype, r)}
      for (i in 1:r) {
        type <- Ftype[i]
        factor[,i] <- getTS(type, TT, ar.coef)
      }  
    }      
    if (length(Fsize)==1)  {Fsize <- rep(Fsize,r)}
    for (i in 1:r) {
      factor[,i] <- factor[,i]*Fsize[i]
    }  
  }  
  
  if (p>0) {
    
    if (time.invariant == FALSE) {
      X <- array(rnorm(N*TT*p), dim = c(TT, N, p))
    } else {
      W <- matrix(rnorm(N*p), nrow = N, ncol = p)
      X <- array(NA, dim = c(TT, N, p))
      for (i in 1:TT) {
        X[i,,] <- W
      }
    }
    
    Xfit <- matrix(0, TT, N)
    for (i in 1:p) {
      Xfit <- Xfit + X[,,i] * beta[i]
    }
 
    Zfit <- matrix(0, TT, N)
    if (is.null(alpha)==FALSE) {            
      for (i in 1:p) {
        Zfit <- Zfit + X[,,i] * matrix(rep(alpha[,i], each = TT), TT, N)
      }
    } else {
      alpha <- matrix(0, N, p)
    }
    
    if (is.null(xi)==TRUE) {  
      xi <- matrix(0, TT, p)
      for (i in c(1:p)) {
        xi[,i] <- getTS(type = "ar1", TT, ar.coef)
      }
    }
    Afit <- matrix(0, TT, N)
    for (i in 1:p) {
      Afit <- Afit + X[,,i] * matrix(rep(xi[,i], N), TT, N)
    }         
  } 
  
  if (force==1|force==3) { 
    unitFE <- matrix(rep(rnorm(N,sd=1),each=TT),TT,N)
  }
  if (force==2|force==3) { 
    ts <- arima.sim(list(order = c(1,0,0), ar = ar.coef), TT) 
    timeFE <- matrix(rep(ts,N),TT,N)
  }
  

  tr.raw <- 0.7*lambda[,1] + 0.3*lambda[,2] + rnorm(N, 0, tr.noise)
  D <- matrix(0, TT, N)
  treat <- rep(0, N) 
  T0 <- rep(TT, N)
  for (i in 1:(ntr-1)) {
    sd.adj <- sqrt(0.7^2 + 0.3^2 + tr.noise^2)
    threshold <- qnorm(tr.threshold[i], mean = 1, sd = sd.adj) 
    tr.id <- which(tr.raw * 10000 >= threshold * 10000 + 1) 
    
    treat[tr.id] <- i             
    T0[tr.id] <- tr.start[i]-1
    for (j in tr.id) {
      D[tr.start[i]:TT,j] <- 1  
    }
  } 
  
  if (error.type == "n") {
    e <- rnorm(TT*N,sd=1)
  } else if (error.type == "g") {
    e <- (rgamma(TT*N, shape = 2, scale = 2) - 4)/sqrt(8)
  }else if (error.type == "t") {
    e <- (rt(TT*N,3))
  }
  e <- matrix(e, TT, N)
  
  
  Y <- matrix(mu, TT, N)
  
  if (r>0) {
    Y <- Y + factor%*%t(lambda)
  }
  if (p>0) {
    Y <- Y + Xfit + Zfit + Afit
  }
 
  if (force==1|force==3) { 
    Y <- Y + 1*unitFE
  }
  if (force==2|force==3) { 
    Y <- Y + 1*timeFE
  }  
  
 
  eff <- matrix(0,TT,N)
  for (i in 1:N){
    eff[, i] <- rnorm(TT, mean = cumsum(D[, i]) * eff.size, sd = eff.noise)
  }
  eff <- eff*D
  Y <- Y + eff
  
 
  Y <- Y + e
  
  

  panel<-as.data.frame(cbind(
    rep(101:(100+N),each=TT),
    rep(1:TT,N),
    rep(1:TT,N),
    c(Y),
    c(e),
    c(eff),
    rep(T0, each = TT)))
  cname<-c("id","time","t","Y","error","eff","T0")
  

  treat <- rep(D[nrow(D),], each = TT)
  panel <- cbind(panel,c(D), treat)
  cname <- c(cname,"D","treat")
  

  if (p>0) {
    for (i in 1:p) {
      panel<-cbind(panel,c(X[,,i]))
      cname<-c(cname,paste("X",i,sep=""))
    }
  }    
  
  if (force==1|force==3) {
    panel <- cbind(panel,c(unitFE))
    cname <- c(cname,"unitFE")
  }
  if (force==2|force==3) {
    panel <- cbind(panel,c(timeFE))
    cname <- c(cname,"timeFE")
  }
  
  if (r>0) {
    for (i in 1:r) {
      panel<-cbind(panel,rep(factor[,i],N))
      cname<-c(cname,paste("F",i,sep=""))
    }
    for (i in 1:r) {
      panel<-cbind(panel,rep(lambda[,i],each=TT))
      cname<-c(cname,paste("L",i,sep=""))
    }
  }
  colnames(panel)<-cname
  
  if (r>0) for (i in 1:r) {
    panel[,paste("FL",i,sep="")]<-panel[,paste("F",i,sep="")]*panel[,paste("L",i,sep="")]
  }
  
  return(panel)
}
