crmh <- function(x)  exp(-x^2/2/s1p^2) * vcrm(x) #posterior
crmht <- function(x)  x * exp(-x^2/2/s1p^2) * vcrm(x) #posterior*x
vcrm <- function(a) { #likleihood
  v <- 1
  for (i in 1:length(x1p))
    v <- v*((x1p[i]^exp(a))^y1p[i])*(((1-w1p[i]*x1p[i]^exp(a))^(1-y1p[i])))
  return(v)
}

titecrm <- function(prior, target, tox, level, n=length(level),
                    weights=NULL, followup=NULL, obswin=NULL, scheme="linear",
                    dosename=NULL, include=1:n, pid=1:n, method="bayes",
                    scale=sqrt(1.34), model.detail=TRUE, patient.detail=TRUE) {
  if (is.null(weights)) {
    if (scheme=="linear") { weights <- followup/obswin; }
    if (scheme=="KM") { } ### Work in progress; }
  }
  if (any(weights>1) | any(weights<0)) stop(" Weights have to be between 0 and 1!")
  if (is.null(pid)) {
    if (! (length(tox)==length(level) & length(tox)==length(weights)))
      stop(" tox, level, and weights are of different lengths!")
  }
  else {
    if (! (length(tox)==length(level) & length(tox)==length(weights) & length(tox)==length(pid)) )
      stop(" pid, tox, level, and weights are of different lengths!")
  }
  weights[tox==1] <- 1
  
  x1p <<- c(prior[level[include]])
  y1p <<- tox[include]
  w1p <<- weights[include]
  s1p <<- scale
  est <- integrate(crmht,-Inf,Inf)[[1]]/integrate(crmh,-Inf,Inf)[[1]]
  ptox <- prior^exp(est)
  rec <- order(abs(ptox-target))[1]
  foo <- list(prior=prior, target=target, tox=tox, level=level, weights=weights,
              followup=followup, obswin=obswin, scheme=scheme, dosename=dosename,
              subset=pid[include], estimate=est,ptox=ptox,
              recommend=rec,method=method,scale=scale,patient.detail=patient.detail,
              model.detail=model.detail,pid=pid,include=include)
  class(foo) <- "mtd"
  foo
}

print.mtd <- function(x,dgt=3,model.detail=x$model.detail,patient.detail=x$patient.detail, ...) {
  cat("Today:", date(),"\n")
  n <- length(x$pid)
  used <- rep(0,n)
  used[x$include] <- 1
  wts <- round(x$weights,digits=dgt)
  ptox <- round(x$ptox, digits=dgt)
  if (is.null(x$followup)) { followup <- rep("N/A",n); }
  else { followup <- round(x$followup,digits=dgt); }
  if (patient.detail) {
    cat("DATA SUMMARY \n")
    cat("PID","\t","Level","\t","Tox.event","\t","f/u","\t","Weight","\t","Included","\n")
    for (i in 1:n)
      cat(x$pid[i],"\t",x$level[i],"\t",x$tox[i],"\t\t",followup[i],"\t",wts[i],"\t\t",used[i],"\n")
  }
  cat("\nToxicity probability update:\n")
  if (is.null(x$dosename)) {
    cat("Level","\t","Prior","\t","n","\t","total.wts","\t","total.tox","\t","Ptox","\n")
    K <- length(prior)
    for (k in 1:K) {
      expt <- which(x$level==k & used==1)
      cat(k,"\t",x$prior[k],"\t",length(expt), "\t",round(sum(x$weights[expt]),digits=dgt),"\t\t",sum(x$tox[expt]),"\t\t",ptox[k],"\n")
    }
    cat("Next recommended dose level:",x$recommend,"\n")
  }
  else {
    cat("Dose","\t\t","Level","\t","Prior","\t","n","\t","total.wts","\t","total.tox","\t","Ptox","\n")
    K <- length(prior)
    for (k in 1:K) {
      expt <- which(x$level==k & used==1)
      cat(x$dosename[k],"\t",k,"\t",x$prior[k],"\t",length(expt), "\t",signif(sum(x$weights[expt]),digits=dgt),"\t\t",sum(x$tox[expt]),"\t\t",ptox[k],"\n")
    }
    cat("Next recommended dose:",x$dosename[x$recommend],"( Level",x$recommend,")","\n")
  }
  cat("Recommendation is based on TITE-CRM for a target at pT =",x$target,"\n")

  if (model.detail) {
    cat("\nEstimation details:\n")
    cat("Empiric dose-toxicity model: p = alpha^{exp(beta)}\n")
    cat("alpha =",prior,"\n")
    if (x$method=="bayes") {
      cat("Normal prior on beta with mean 0 and variance",x$scale^2,"\n")
      cat("Posterior mean of beta:",x$estimate,"\n")
    }
    else if (x$method=="mle") {
	cat("Mle of beta:",x$estimate,"\n")
    }
  }
}

plot.mtd <- function(x,ask=FALSE,den=40, ...) {
  p <- x$ptox
  p0 <- x$prior
  target <- x$target
  K <- length(p)
  status <- x$tox
  followup <- x$followup
  obswin <- x$obswin
  scheme <- x$scheme
  
  if (!(is.null(obswin) | is.null(followup)))  {
    hor <- seq(0,obswin,len=den)
    ver <- hor/obswin
    if (scheme=="KM") {
      support <- sort(followup[status==1])
      z <- length(support)
      if (z) {
        for (j in 1:den) {
          m <- length(support[support<=hor[j]])
          if (!m)  ver[j] <- hor[j] / support[1] / (z+1)
          else if (m==z)
            ver[j] <- (z + (hor[j]-support[z])/(obswin-support[z])) / (z+1)
          else
            ver[j] <- (m + (hor[j]-support[m])/(support[(m+1)]-support[m])) / (z+1)
        }
      }
    }
  }
  
  choices <- c("All","Dose-response curves","Weight function")
  choices <- substring(choices,1,40)
  tmenu <- paste("plot:",choices)
  maxpick <- pick <- 2
  ask.now <- ask
  while (pick <= length(tmenu) + maxpick) {
    if(ask.now)
      pick <- menu(tmenu,title="\nMake a plot selection (or 0 to exit):\n") + 1
    switch(pick,return(invisible(x)),ask.now <- FALSE,
           {
             par(las=1,mgp=c(2.8,1,0))
             plot(1:K,p,type="l",xlab="Dose level",ylab="Probability of toxicity",ylim=c(0,1))
             points(1:K,p,pch="X")
             lines(1:K,p0,lty=2)
             points(1:K,p0,pch="O")
             abline(h=target,lwd=2)
             mtext("Prior (dotted) and updated (solid) dose-toxicity curves",line=0.5)
           }
           ,
           {
             if (! (is.null(followup) | is.null(obswin)) & scheme=="KM") {
               plot(hor,ver,type="n",xlab="follow-up time of patient",ylab="weight")
               lines(hor,ver,lty=5)
               mtext("The adaptive weight function",line=0.5)
               if (z) points(myjitter(support),rep(0,z),pch=16)
             }
             else {
               cat("No plot output of weight function:\n")
               cat("Weights are provided by user, or are calculated by linear weight function.\n")
             }
           }
           )
    if (!ask.now)  pick <- pick + 1
    if (pick == length(tmenu) + 2) ask.now <- ask
  }
  return(invisible(x))
}

myjitter <- function(x,factor=1) {
  z <- diff(range(x[!is.na(x)]))
  z <- factor * (z/50)
  if (!z)  z <- abs(mean(x)/50)
  x + runif(length(x),-z,z)
}

cohere <- function(prior,target,n,x0,method="bayes",scale=sqrt(1.34), detail=TRUE) {
  coh <- rep(TRUE,(n-1)); vlevel <- rep(NA,(n-1))
  for (i in 1:(n-1)) {
    level <- x0[1:i]
    y <- c(rep(0,(i-1)),1)
    obj <- titecrm(prior,target,y,level,weights=rep(1,i),method=method,scale=scale)
    est <- obj$estimate
    pest <- prior^exp(est)
    cur <- order(abs(pest-target))[1]
    if (cur > level[i]) { 
      coh[i] <- FALSE; vlevel[i] <- cur;
      if (detail) {
        cat("Incoherent escalation occurs after n =",i,"patients.\n")
        cat("Level","\t",level,"\n"); cat("tox","\t",y,"\n");
        cat("Recommended level:",cur,"\n\n")
      }
    }
  }
  if (all(coh)) { msg <- "Coherent" }
  else {
    ind <- which(!coh)
    m <- length(ind)
    msg <- "Incoherent! Take a less conservative initial design or a small target DLT rate." 
  }
  msg
}

titesim2 <- function(PI,prior,target,n,x0,obswin=1,tgrp=obswin,rate=1,accrual="fixed",surv="uniform",scheme="linear",method="bayes",scale=sqrt(1.34),seed=1099) {
	set.seed(seed)
	if (length(x0)!=n) stop("Sample size of initial design is not equal to n!")
 	est <- u <- y <- level <- arrival <- rep(NA,n)
  	if (accrual=="fixed") next.arrival <- obswin/rate
  	else if (accrual=="poisson") next.arrival <- rexp(1,rate/obswin)	
  	m <- 1
	while (TRUE) {
		level[m] <- cur <- x0[m]
		if (is.na(arrival[m])) arrival[m] <- next.arrival
		if (is.na(y[m])) {
			if (surv=="uniform") {
				y[m] <- ynew <- rbinom(1,1,PI[cur]); 
				if (ynew) unew <- runif(1,0,obswin)
				else unew <- Inf
				u[m] <- unew; utox <- u + arrival;			
			}
		}
		if (accrual=="fixed") next.arrival <- next.arrival + obswin/rate
		else if (accrual=="poisson") next.arrival <- next.arrival + rexp(1,rate/obswin)
		B <- rep(0,n);  B[utox<=next.arrival] <- 1;
		if (sum(B)>0 | m==(n-1)) break
		if (x0[m+1]==cur | (next.arrival-arrival[m])>=tgrp)  m <- m+1
	}
	if (m==(n-1)) {
		if (sum(B)==0) { cur <- x0[n]; }
		else {
			censor <- pmin(next.arrival,utox) - arrival;
			followup <- pmin(censor,obswin)
			obj <- titecrm(prior,target,B[1:m],level[1:m],followup=followup[1:m],obswin=obswin,scheme=scheme,method=method,scale=scale)
			cur <- obj$recommend; est[m] <- obj$estimate;
		}
		arrival[n] <- next.arrival; level[n] <- cur;
		if (surv=="uniform") {
			y[n] <- ynew <- rbinom(1,1,PI[cur])
			if (ynew) unew <- runif(1,0,obswin)
			else unew <- Inf
			u[n] <- unew; utox <- u + arrival
		}
	}
	else {
		censor <- pmin(next.arrival,utox)-arrival; followup <- pmin(censor,obswin);
		obj <- titecrm(prior,target,B[1:m],level[1:m],followup=followup[1:m],obswin=obswin,scheme=scheme,method=method,scale=scale)
		cur <- obj$recommend; est[m] <- obj$estimate;
		if (accrual=="fixed") next.arrival <- next.arrival + obswin/rate
		else next.arrival <- next.arrival + rexp(1,rate/obswin)
		
		for (i in (m+1):(n-1)) {
			arrival[i] <- next.arrival; level[i] <- cur;
			if (surv=="uniform") {
				y[i] <- ynew <- rbinom(1,1,PI[cur])
				if (ynew) unew <- runif(1,0,obswin)
				else unew <- Inf
				u[i] <- unew; utox <- u + arrival;
			}
			if (accrual=="fixed") next.arrival <- next.arrival + obswin/rate
		    else if (accrual=="poisson") next.arrival <- next.arrival + rexp(1,rate/obswin)
   		 	B <- rep(0,n);  B[utox<=next.arrival] <- 1;
    		censor <- pmin(next.arrival,utox) - arrival; followup <- pmin(censor,obswin);
    		obj <- titecrm(prior,target,B[1:i],level[1:i],followup=followup[1:i],obswin=obswin,scheme=scheme,method=method,scale=scale)
    		cur <- obj$recommend; est[i] <- obj$estimate
		}
		arrival[n] <- next.arrival; level[n] <- cur;
		if (surv=="uniform") {
			y[n] <- ynew <- rbinom(1,1,PI[cur])
			if (ynew) unew <- runif(1,0,obswin)
			else unew <- Inf
			u[n] <- unew; utox <- u + arrival
		}
	}
	obj <- titecrm(prior,target,y,level,weights=rep(1,n),scheme=scheme,method=method)
  	cur <- obj$recommend
  	est[n] <- obj$estimate

  	foo <- list(PI=PI,prior=prior,target=target,tox=y,level=level,arrival=arrival,ttox.pt=u,ttox.cal=utox,obswin=obswin,tgrp=tgrp,
              initD=x0,estimate=est,recommend=cur,scheme=scheme,method=method,scale=scale,accrual=accrual,rate=rate,surv=surv,nstage=2)
  	class(foo) <- "sim"
  	foo
} 

titesim1 <- function(PI,prior,target,n,obswin=1,rate=1,accrual="fixed",surv="uniform",scheme="linear",method="bayes",scale=sqrt(1.34),seed=1099) {
  set.seed(seed)
  est <- u <- y <- level <- arrival <- NULL
  if (accrual=="fixed") next.arrival <- obswin/rate
  else if (accrual=="poisson") next.arrival <- rexp(1,rate/obswin)
  cur <- order(abs(prior-target))[1]
  
  for (i in 1:(n-1)) {
    arrival <- c(arrival, next.arrival)
    level <- c(level, cur)
    if (surv=="uniform") {
      ynew <- rbinom(1,1,PI[cur])
      if (ynew) unew <- runif(1,0,obswin)
      else unew <- Inf
      y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
    }
    if (accrual=="fixed") next.arrival <- next.arrival + obswin/rate
    else if (accrual=="poisson") next.arrival <- next.arrival + rexp(1,rate/obswin)
    B <- rep(0,length(y));  B[utox<=next.arrival] <- 1;
    censor <- pmin(next.arrival,utox) - arrival
    followup <- pmin(censor,obswin)
    obj <- titecrm(prior,target,B,level,followup=followup,obswin=obswin,scheme=scheme,method=method,scale=scale)
    cur <- obj$recommend
    est <- c(est, obj$estimate)
  }
  arrival <- c(arrival, next.arrival)
  level <- c(level,cur)
  if (surv=="uniform") {
    ynew <- rbinom(1,1,PI[cur])
    if (ynew) unew <- runif(1,0,obswin)
    else unew <- Inf
    y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
  }
  obj <- titecrm(prior,target,y,level,weights=rep(1,n),scheme=scheme,method=method)
  cur <- obj$recommend
  est <- c(est, obj$estimate)

  foo <- list(PI=PI,prior=prior,target=target,tox=y,level=level,arrival=arrival,ttox.pt=u,ttox.cal=utox,obswin=obswin,
              estimate=est,recommend=cur,scheme=scheme,method=method,scale=scale,accrual=accrual,rate=rate,surv=surv,nstage=1)
  class(foo) <- "sim"
  foo
} 

plot.sim <- function(x, ...) {
	n <- length(x$level)
	tevent <- c(x$arrival, x$ttox.cal)
	pid <- rep(1:n,2)
	event <- c(rep("enrol",n), rep("TOX",n))
	level <- c(x$level, x$level)
	est <- c(0, x$estimate[1:(n-1)], rep(NA,n))
	o <- order(tevent)
	tevent <- tevent[o]; pid <- pid[o]; event <- event[o]; level <- level[o]; est <- est[o];
	ind <- which(tevent<Inf)
	tevent <- tevent[ind]; pid <- pid[ind]; event <- event[ind];	level <- level[ind]; est <- est[ind];
	m <- length(ind)
   
	xmax <- max( max(tevent), max(x$arrival)+x$obswin);
	xmax <- max(tevent)
	leveljit <- level; leveljit[event=="TOX"] <- leveljit[event=="TOX"] + 0.1
	plot(tevent,level,ylim=c(1,length(x$prior)),type="n",xlab="Study time",ylab="Dose level",xlim=c(0,xmax))
	text(tevent,(leveljit),as.character(pid),cex=0.7)
	circ <- rep(" ",m)
	circ[event=="TOX"] <- "O"
	text(tevent,(leveljit+0.02),circ,cex=1.6)
	mtext("Each number respresents a patient",line=2)
	mtext("An uncircled number indicates time of arrival, circled indicates time of toxicity",line=0.5)
}

print.sim <- function(x,dgt=3, ...) {
	n <- length(x$level)
	tevent <- signif( c(x$arrival, x$ttox.cal), digits=dgt)
	pid <- rep(1:n,2)
	event <- c(rep("enrol",n), rep("TOX",n))
	level <- c(x$level, x$level)
	est <- signif( c(0, x$estimate[1:(n-1)], rep(NA,n)), digits=dgt)
	o <- order(tevent)
	tevent <- tevent[o]; pid <- pid[o]; event <- event[o]; level <- level[o]; est <- est[o];
	ind <- which(tevent<Inf)
	tevent <- tevent[ind]; pid <- pid[ind]; event <- event[ind];	level <- level[ind]; est <- est[ind];
	m <- length(ind)
	cat("TRIAL SUMMARY on calendar time\n")
	cat("Time","\t","PID","\t","Event","\t","Level","\t","*Estimate","\n")
	for (j in 1:m) cat(tevent[j],"\t",pid[j],"\t",event[j],"\t",level[j],"\t",est[j],"\n")
	
	pid <- 1:n
	arrival <- signif( x$arrival, digits=dgt)
	level <- x$level
	tox <- x$tox
	ttox <- signif( x$ttox.pt, digits=dgt)
	est <- signif( c(0, x$estimate[1:(n-1)]), digits=dgt)
	cat("\nPATIENT SUMMARY on patient time\n")
	cat("PID","\t","Arrive","\t","Level","\t","Tox","\t","T.tox","\t","*Estimate","\n")
	for (i in 1:n) cat(pid[i],"\t",arrival[i],"\t\t",level[i],"\t",tox[i],"\t",ttox[i],"\t",est[i],"\n")
	cat("*Estimate used to calculate the dose for the patient; \n")
	cat(" Assume real-time update.\n")

	ntox <- expt <- rep(NA,length(x$prior))
	ptox <- signif( x$prior^exp(x$estimate[n]), digits=dgt)
	mtd <- order(abs(x$PI-x$target))[1]
	cat("\nDOSE and OVERALL SUMMARY\n")
	cat("Level","\t","Ptrue","\t","Prior","\t","n","\t","ntox","\t","ptox","\n")
	for (k in 1:length(x$prior)) {
		expt[k] <- length(which(x$level==k))
		ntox[k] <- length(which(x$level==k & x$tox==1))
		cat(k,"\t",x$PI[k],"\t",x$prior[k],"\t",expt[k],"\t",ntox[k],"\t",ptox[k],"\n")
	}
	cat("True MTD = level", mtd, "    Recommended MTD = level",x$rec,"\n")
	if (x$nstage==1) { cat("Recommendation is based on TITE-CRM with",n,"patients and target pT",x$target,"\n"); }
	else if (x$nstage==2) {
		cat("Recommendation is based on two-stage TITE-CRM with",n,"patients and target pT",x$target,"\n")
		cat(" with initial design:\n")
		cat("\t\t\t","Level","\t",names(table(x$initD)),"\n")
		cat("\t\t\t","Size","\t",as.numeric(table(x$initD)),"\n")
		cat("Before switching to TITE-CRM, escalation occurs only when\n")
		cat(" all patients in previous cohort have been followed for at least",x$tgrp,"time units\n")
	}
	cat("Final estimate of model parameter:",signif(x$estimate[n],digits=dgt),"\n")
	if (x$method=="bayes")	cat("using bayes estimation with prior mean 0 and variance",x$scale^2,"\n")
	cat("Time-to-event is modeled as",x$surv,"\n")
	cat("Patient arrival is modeled as",x$accrual,"\n")
	cat(" with rate",x$rate,"patients per",x$obswin,"time units (= observation window)\n")
}

drplot <- function(x,n=length(x$level),include=n) {
  p0 <- x$prior
  target <- x$target
  if (length(x$level) != length(x$estimate)) { pid <- 1; }
  else { pid <- (1:n)[include]; }
  est <- x$estimate
  K <- length(x$prior)
  plot(1:K, p0,type="n",xlab="Dose level",ylab="Probability of toxicity",ylim=c(0,1))
  lines(1:K,p0,lty=2)
  points(1:K, p0,pch="o")
  abline(h=target,lwd=2)
  for (j in pid) {
    ptox <- prior^exp(est[j])
    lines(1:K,ptox,col=j)
    sym <- as.character(j)
    text(1:K,ptox,sym)
  }	
}

crmsens <- function(prior,target,eps=1e-8,maxit=100,incr=0.1,LB=0,UB=Inf,detail=FALSE) {
  K <- length(prior)
  ini <- log(target)/log(prior)
  theta <- ini[1:(K-1)]   # equal to exp(beta)

  for (k in 1:(K-1)) {
    notdone <-  TRUE
    lb <- ub <- ini[k]
    lhs <- (prior[k]^theta[k] + prior[k+1]^theta[k])/2
    if (lhs < target) {
      while (notdone) {
        lb <- lb-incr
        lhs <- (prior[k]^lb + prior[k+1]^lb)/2
        if (lhs>=target)  notdone <-  FALSE
        theta[k] <- lb
      }
    }
    else if (lhs > target) {
      while (notdone) {
        ub <- ub+incr
        lhs <- (prior[k]^ub + prior[k+1]^ub)/2
        if (lhs<=target) notdone <- FALSE
        theta[k] <- ub
      }
    }
    
    if (abs(lhs-target)>eps) {
      notdone <- TRUE
      iter <- 0
      mid <- (lb+ub)/2
      lhs <- (prior[k]^mid + prior[k+1]^mid)/2
      while (notdone) {
        if (lhs<target)  ub <- mid
        else  lb <- mid
        mid <- (ub + lb) / 2
        lhs <- (prior[k]^mid + prior[k+1]^mid)/2
        if (iter>maxit || abs(lhs-target)<=eps)  notdone <- FALSE
        iter <- iter + 1
      }
      theta[k] <- mid
    }
  }

  val <- rep(NA,K*2)
  dim(val) <- c(K,2)
  val[1,1] <- LB
  val[K,2] <- UB
  for (k in 1:(K-1)) {
    val[k,2] <- val[(k+1),1] <- theta[k]
  }

  IRm <- IRp <- rep(NA,K*2)
  dim(IRm) <- dim(IRp) <- c(K,2)
  IRm[2:K,2] <- IRp[1:(K-1),1] <- target
  for (l in 1:K) {
    if (l>1) IRm[l,1] <- prior[l-1]^val[l,1]
    if (l<K)  IRp[l,2] <- prior[l+1]^val[l,2]
  }
  
  a <- list(homeset=val,iint=cbind(IRm,IRp),target=target,prior=prior,K=K,detail=detail)
  class(a) <- "homesets"
  a
}

print.homesets <- function(x,dgt=3, ...) {
  if (x$detail) {
    cat("\nHome sets for the model setup:")
    dimnames(x$homeset) <- list(NULL,c("",""))
    print(x$homeset)
  }

  cat("\nTrue","\t\t","Indifference interval","\n")
  cat("MTD","\t\t","Lower","\t\t","Upper","\n")
  for (k in 1:x$K) cat(k,"\t\t",round(x$iint[k,1],digits=dgt),"\t\t",round(x$iint[k,4],digits=dgt),"\n")
  cat("With this model, the CRM will eventually choose a dose with \n")
  cat(" ptox between", signif(min(x$iint[(2:x$K),1]),digits=dgt),"and",signif(max(x$iint[(1:(x$K-1)),4]),digits=dgt),
      "while targeting at",x$target,"\n\n")
  cat("Consistency will hold if\n")
  cat("(i) the dose below the MTD does not exceed the Lower limit\n")
  cat("(ii) the dose above the MTD is more toxic than the Upper limit\n\n")
}


#.First.lib <- function(lib,pkg) {
#  library.dynam("titecrm",pkg,lib)
#  provide(titecrm)
#}
