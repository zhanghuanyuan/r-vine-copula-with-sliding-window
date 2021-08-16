require("pracma")
require("fGarch")
require("rugarch")
require("VineCopula")
require("copula")
require(tseries)
require(FinTS)

RollCheck<-function(data,roll_long,roll_interval){
  #data，sliding window，step length
  #ensure these 3 variables are available
  
  if(length(data[1,])==1){
    data_set=0
    stop('Single column data cannot be calculated.')
  }
  else{
    data_set<-length(data[1,])
  }
  if(length(data[,1])<=roll_long){
    data_number
    stop('the length of sliding window cannot greater than the length of data.')
  }
  else{
    data_number=length(data[,1])
  }
  if(roll_long<=0){
    stop('the length of sliding window must greater than 0.')
  }
  else{
    a<-roll_long
  }
  
  if(roll_interval<=0){
    stop('the step length must greater than 0.')
  }
  else{
    b<-roll_interval
  }
  
  if(floor((data_number-roll_long)/roll_interval)==(data_number-roll_long)/roll_interval){
    interval_all<-(data_number-roll_long)/roll_interval+1
  }
  else{
  interval_all<-(floor((data_number-roll_long)/roll_interval))+2
  }

  return(list(data_set=data_set,
              data_number=data_number,
              roll_interval=b,
              roll_long=a,
              interval_all=interval_all))
}

panduan<-function(x){
  if(x<=0.01){
    a="***"
  }else if(x<=0.05){
    a="**"
  }else if(x<=0.1){
    a="*"
  }else{
    a=""
  }
  return(a)
}

market_simp_name<-function(x){
    name<-switch(x,
        "SH",
        "US",
        "FR",
        "UK",
        "HK",
        "JP",
        "CH",
        "SG",
        "TW")
}

summarise<-function(data,RC){
  data_set=as.numeric(RC$data_set)
  summ<-matrix(data=NA,nrow=data_set,ncol=10)
  colnames(summ)<-c("name","mean","max","min","sd","skewness","kurtosis","ks","adf","Arch-LM")
  
  for(i in 1:data_set){
    #name
    summ[i,1]<-market_simp_name(i)
    #data
    data_1<-as.numeric(c(unlist(data[,(i)])))
    #mean,max,min,sd
    a<-as.numeric(summary(data_1))
    summ[i,2]<-a[4]
    summ[i,3]<-a[6]
    summ[i,4]<-a[1]
    summ[i,5]<-sd(data_1)
    #skewness&kurtosis
    summ[i,6]<-skewness(data_1)
    summ[i,7]<-kurtosis(data_1)
    #ks_test&adf_test
    ks_result_test<-ks.test(jitter(data_1),"pnorm")
    adf_result_test<-adf.test(data_1)
    summ[i,8]<-paste(round(ks_result_test[["statistic"]],3),panduan(ks_result_test[["p.value"]]))
    summ[i,9]<-paste(round(adf_result_test[["statistic"]],3),panduan(adf_result_test[["p.value"]]))
    #ARCH-LM
    archlm<-ArchTest(data_1,1)
    summ[i,10]<-paste(round(archlm[["statistic"]],3),panduan(archlm[["p.value"]]))
    }

  summ
}

edf<-function(bzh){
  cdf<-matrix(data=NA,nrow=length(bzh),ncol=1)
  cdf<-rank(bzh)
  for(i in 1:length(bzh)){
    cdf[i]<-cdf[i]/length(bzh)
  }
  cdf
}

copula_jisuan<-function(data){#arma(0,1)garch(1,1) as default
  fit<-garchFit(~arma(0,1)+garch(1,1),data,cond.dist = "std",trace=FALSE)
  bzh<-residuals(fit)/volatility(fit)
  bzh<-as.numeric(bzh)
  jisuan_copula<-edf(bzh)
}

copula_gjr<-function(data){#arma(0,1)gjrgarch(1,1) as default
  spec<-ugarchspec(
    variance.model = list(model = "gjrGARCH", 
                          garchOrder = c(1, 1), #garch(1,1)
                          submodel = NULL, 
                          external.regressors = NULL, 
                          variance.targeting = FALSE),
    mean.model = list(armaOrder = c(0, 1), #arma(0,1)
                      include.mean = TRUE, 
                      archm = FALSE, 
                      archpow = 1, 
                      arfima = FALSE, 
                      external.regressors = NULL, 
                      archex = FALSE),
    distribution.model = "norm"
  )
  gjrfit<-ugarchfit(spec,data=data,solver="solnp")
  bzh<-gjrfit@fit[["residuals"]]/gjrfit@fit[["sigma"]]
  bzh<-as.numeric(bzh)
  jisuan_gjr<-edf(bzh)
}

CopulaData<-function(data,RC,t=1){#t=1,garch;t=2,gjr-garch
  data_set=as.numeric(RC$data_set)
  data_number<-as.numeric(RC$data_number)
  data_copula<-matrix(data=NA,nrow=data_number,ncol=data_set)
  if(t==1){
    for(i in 1:data_set){
      data_copula[,i]<-copula_jisuan(data[,i])
    }
  }
  if(t==2){
    data_copula<-matrix(data=NA,nrow=data_number,ncol=data_set)
    for(i in 1:data_set){
      data_copula[,i]<-copula_gjr(data[,i])
    }
  }
  else{
    stop('t must be 1 or 2.')
  }
  data_copula
}

RollRVM<-function(data_copula,RC){
  data_set=as.numeric(RC$data_set)
  data_number<-as.numeric(RC$data_number)
  roll_interval=as.numeric(RC$roll_interval)
  roll_long=as.numeric(RC$roll_long) 
  interval_all=as.numeric(RC$interval_all)
  for(j in 1:interval_all){

    if(j==1){
      roll_info<-matrix(data=NA,nrow=interval_all,ncol=3)
      data_copula_roll<-matrix(data=NA,nrow=roll_long,ncol=data_set)
      time1<-Sys.time()
      time_est<-"(Estimating time)"
    } 

    if(j!=1){
       ttl<-ttl_sec*(interval_all-j+1)
   	   h<-ttl%/%3600
       m<-(ttl%%3600)%/%60
       s<-(ttl%%3600)%%60
       time_est<-paste("(Estimated to take ",h,"h",m,"m",s,"s.)")
    }
  	
    nameprocess<-paste("Calculating ",j," of ",interval_all,"...",sep="")
    print(paste(nameprocess,time_est))
    
    if(j!=interval_all){
      data_copula_roll<-data_copula[(roll_interval*(j-1)+1):(roll_interval*(j-1)+roll_long),]
    }
    else{
      data_copula_roll<-data_copula[(data_number-roll_long+1):(data_number),]
    }

    RVCopula<-RVineStructureSelect(data_copula_roll,familyset=c(1:6))
    roll_info[j,]<-AIC_Loglike(RVCopula)
    
    if(j==interval_all){
      colnames(roll_info)<-c("aic","bic","loglik")
      print("Completed!")
    }

    if(j==1){
    	a<-TimeDiff(time1)
    	hour<-as.numeric(substr(a,1,2))
    	min<-as.numeric(substr(a,4,5))
    	sec<-as.numeric(substr(a,7,8))
    	ttl_sec<-(hour*60*60)+(min*60)+sec    	
    }
  }
  roll_info
}

RVCopula<-function(data_copula_roll){
  RVM<-RVineStructureSelect(data_copula_roll,c(1:6))
  RVM1<-RVineCopSelect(data_copula_roll,familyset = c(1:6),RVM$Matrix)
  RVM2<-RVineSeqEst(data_copula_roll,RVM1,method="mle")

  mle<-RVineMLE(data_copula_roll,RVM1)
  
  vine_matrix<-RVineMatrix(mle$RVM$Matrix,
                           family = mle$RVM$family,
                           par=mle$RVM$par,
                           par2=mle$RVM$par2,
                           names=mle$RVM$names)
  return(vine_matrix)
}

RollFixStru<-function(data_copula,RC){
  data_set=as.numeric(RC$data_set)
  data_number=as.numeric(RC$data_number)
  roll_interval=as.numeric(RC$roll_interval)
  roll_long=as.numeric(RC$roll_long)
  interval_all=as.numeric(RC$interval_all) 
  RVM<-RVineStructureSelect(data_copula,c(1:6))
  fami<-vector(mode="numeric",length=0)
  tauv<-vector(mode="numeric",length=0)

  for(j in 1:interval_all){
    nameprocess<-paste("Calculating ",j," of ",interval_all,"...",sep="")
    print(nameprocess)
    
    if(j==1){
      data_copula_roll<-matrix(data=NA,nrow=roll_long,ncol=data_set)
    }
    
    if(j!=interval_all){
      data_copula_roll<-data_copula[(roll_interval*(j-1)+1):(roll_interval*(j-1)+roll_long),]
    }
    else{
      data_copula_roll<-data_copula[(data_number-roll_long+1):(data_number),]
    }

    RVM1<-RVineCopSelect(data_copula_roll,familyset = c(1:6),RVM$Matrix)
    vine_matrix<-RVineMatrix(RVM$Matrix,
                           family = RVM1$family,
                           par=RVM1$par,
                           par2=RVM1$par2,
                           names=RVM1$names)

    fami<-append(fami,as.vector(extract_t(vine_matrix$family)))                
    tauv<-append(tauv,as.vector(extract_t(vine_matrix$tau)))

    if(j==interval_all){
      print("Completed!")
    }
  }

  z=0
  family_seq<-vector(mode="numeric",length=0)
  tau_seq<-vector(mode="numeric",length=0)
  for(x in 1:data_set-1){
    z=z+1

    for(y in 1:interval_all){
      family_seq=append(data,fami[data_set*(y-1)+z])
      tau_seq=append(data,fami[data_set*(y-1)+z])
    }
  }

  list(family=family_seq,tau=tau_seq)
}

AIC_Loglike<-function(RVM){
  roll_info<-matrix(data=NA,nrow=1,ncol=3)
  roll_info[1,3]<-RVM$logLik
  roll_info[1,1]<-RVM$AIC
  roll_info[1,2]<-RVM$BIC
  roll_info
}

RVMtree1<-function(RVM_matrix){
  RVMtree1_col<-matrix(data=NA,nrow=(length(RVM_matrix[,1])-1),ncol=2)
  RVMtree1_col[,1]<-t(RVM_matrix[(length(RVM_matrix[,1])),1:(length(RVM_matrix[,1])-1)])
  for (i in 1:(length(RVM_matrix[,1])-1)){
    RVMtree1_col[i,2]<-RVM_matrix[i,i]
  }

  RVMtree1_net<-matrix(data=0,
                       nrow=(length(RVM_matrix[,1])),
                       ncol=(length(RVM_matrix[,1])))
  for (i in 1:(length(RVMtree1_col[,1]))){
    xmin<-min(RVMtree1_col[i,1],RVMtree1_col[i,2])
    xmax<-max(RVMtree1_col[i,1],RVMtree1_col[i,2])
    RVMtree1_net[xmin,xmax]<-RVMtree1_net[xmin,xmax]+1
    RVMtree1_net[xmax,xmin]<-RVMtree1_net[xmax,xmin]+1
  }

  return(list(RVMtree1_col=RVMtree1_col,
              RVMtree1_net=RVMtree1_net))
}

extract_t<-function(matrix)
{
  extract_t<-t(matrix[(length(matrix[,1])),1:(length(matrix[,1])-1)])
}

Extract<-function(RVM,RC){
  data_set=as.numeric(RC$data_set)
  info<-matrix(data=NA,nrow=(data_set-1),ncol=7)
  info[,1:2]<-RVMtree1(RVM$Matrix)$RVMtree1_col  #knot
  info[,3]<-extract_t(RVM$family)                #family
  info[,5]<-extract_t(RVM$tau)                   #tau
  info[,6]<-extract_t(RVM$par)                   #par
  info[,7]<-extract_t(RVM$par2)                  #par2
  for(f in 1:(data_set-1)){
    if(as.numeric(info[f,3]) <= 14){
      info[f,4]<-switch(as.numeric(info[f,3]),
                      "Gaussian",
                      "Student",
                      "Clayton",
                      "Gumbel",
                      "Frank",
                      "Joe",
                      "BB1",
                      "BB6",
                      "BB7",
                      "BB8",
                      11,
                      12,
                      "RotClayton",
                      "RotGumbel")

    }
    else{
      info[f,4]=info[f,3]
    }
    
  }
  
  colnames(info)<-c("knot1","knot2","family","family name","tau","par","par2")
  info
}

Centrality<-function(RVMtree1_net){
  point_centrality_direct<-matrix(data=0,nrow=length(RVMtree1_net[,1]),ncol=2)
  point_centrality_indirect<-matrix(data=0,nrow=length(RVMtree1_net[,1]),ncol=2)
  for (i in 1:(length(RVMtree1_net[,1]))){
    point_centrality_direct[i,1]<-i
    point_centrality_indirect[i,1]<-i
    for (j in 1:(length(RVMtree1_net[1,]))){
      point_centrality_direct[i,2]<-point_centrality_direct[i,2]+RVMtree1_net[i,j]
    }
  }
  point_centrality_indirect[,2]<-point_centrality_direct[,2]/length(RVMtree1_net[,1])

  return(list(point_centrality_direct=point_centrality_direct,
              point_centrality_indirect=point_centrality_indirect))
}

CoVaRCalculate<-function(data,knot,par,par2,type){
    data_covar<-as.numeric(data[,knot])
    spec<-ugarchspec(
      variance.model = list(model = "gjrGARCH", 
                            garchOrder = c(1, 1), #garch(1,1)
                            submodel = NULL, 
                            external.regressors = NULL, 
                            variance.targeting = FALSE),
      mean.model = list(armaOrder = c(0, 1), #arma(0,1)
                        include.mean = TRUE, 
                        archm = FALSE, 
                        archpow = 1, 
                        arfima = FALSE, 
                        external.regressors = NULL, 
                        archex = FALSE),
      distribution.model = "norm"
    )
    gjrfit<-ugarchfit(spec,data=data_covar,solver="solnp") 
    covar_.05<-as.matrix(CoVaR(0.05,0.05,par=par,par2=par2,dof=0,gamma=0,
             cond.mean=mean(data_covar),cond.sigma=gjrfit@fit[["sigma"]],
             dist="gauss",type=type)$CoVaR)
    covar_.5<-as.matrix(CoVaR(0.5,0.5,par=par,par2=par2,dof=0,gamma=0,
             cond.mean=mean(data_covar),cond.sigma=gjrfit@fit[["sigma"]],
             dist="gauss",type=type)$CoVaR)
    delta<-matrix(data=NA,length(covar_.05),ncol=1)
    for(j in 1:length(covar_.05)){
        delta[j,]<-(covar_.05[j,]-covar_.5[j,])
    }
    delta
}

DeltaCoVaR<-function(data,info,outputpath){
  for(i in 1:length(info[,1])){
    nameprocess<-paste("Calculating ",i," of ",
                   length(info[,1]),
                   "...(",i,"/",
                   length(info[,1]),")",
                   sep="")
    print(nameprocess)

    knot1<-as.numeric(info[i,1])
    knot2<-as.numeric(info[i,2])
    par<-as.numeric(info[i,6])
    par2<-as.numeric(info[i,7])
    type<-as.character(info[i,4])
    
    #delta<-deltacovar(data,knot1,knot2,par,par2,type)
    delta_2_1<-CoVaRCalculate(data,knot1,par,par2,type)
    delta_1_2<-CoVaRCalculate(data,knot2,par,par2,type)
    
    covar_name1<-market_simp_name(knot1)
    covar_name2<-market_simp_name(knot2)
    
    covar2_name<-paste(covar_name1,covar_name2,sep="_to_")
    covar1_name<-paste(covar_name2,covar_name1,sep="_to_")
    
    delta<-data.frame(delta_1_2,delta_2_1)
    colnames(delta)<-c(covar1_name,covar2_name)
    covar_name<-paste(outputpath,"/CoVaR_",covar_name1,covar_name2,".csv",sep="")

    write.csv(delta,file=covar_name)
    if(i==length(info[,1])){
      print("Completed!")
    }
  }
}

RollCor<-function(knot1,knot2,RC){
  interval_all=as.numeric(RC$interval_all)
  roll_interval=as.numeric(RC$roll_interval) 
  roll_long=as.numeric(RC$roll_long) 
  data_number=as.numeric(RC$data_number)
  cor<-matrix(data=NA,nrow=interval_all,ncol=1)
  
  for(i in 1:interval_all){
    if(i!=interval_all){
      x<-as.numeric(data[(roll_interval*(i-1)+1):(roll_interval*(i-1)+roll_long),knot1])
      y<-as.numeric(data[(roll_interval*(i-1)+1):(roll_interval*(i-1)+roll_long),knot2])
    }
    else{
      x<-as.numeric(data[(data_number-roll_long+1):(data_number),knot1])
      y<-as.numeric(data[(data_number-roll_long+1):(data_number),knot2])
    }
    
    cor[i,]<-cor.test(x,y,method="kendall")[["estimate"]]
  }
  cor
}

Rollcopulafamily<-function(knot1,data_copula,RC){
  data_set=as.numeric(RC$data_set)
  data_number<-as.numeric(RC$data_number) 
  roll_interval=as.numeric(RC$roll_interval)
  roll_long=as.numeric(RC$roll_long) 
  interval_all=as.numeric(RC$interval_all) 

  for(j in 1:interval_all){
    nameprocess<-paste("Calculating ",j," of ",interval_all,"...(",j,"/",interval_all,")",sep="")
    print(nameprocess)
    
    if(j==1){
      roll_info<-matrix(data=NA,nrow=interval_all,ncol=3)
      data_copula_roll<-matrix(data=NA,nrow=roll_long,ncol=data_set)
    }
    
    if(j!=interval_all){
      data_copula_roll<-data_copula[(roll_interval*(j-1)+1):(roll_interval*(j-1)+roll_long),]
    }
    else{
      data_copula_roll<-data_copula[(data_number-roll_long+1):(data_number),]
    }
    RVCopula<-RVineStructureSelect(data_copula_roll,familyset=c(1:6)) 
    
    info<-Extract(RVCopula,data_set)

    for(i in 1:(data_set-1)){
      if(as.numeric(info[i,1])==knot1){
        info_SHHK[j,1]<-info[i,1]
        info_SHHK[j,2]<-info[i,2]
        info_SHHK[j,3]<-info[i,4]
        info_SHHK[j,4]<-info[i,5]
      }
      
      else if(as.numeric(info[i,2])==knot1){
        info_SHHK[j,1]<-info[i,2]
        info_SHHK[j,2]<-info[i,1]
        info_SHHK[j,3]<-info[i,4]
        info_SHHK[j,4]<-info[i,5]
      }
    }
    
    if(j==interval_all){
      print("Completed!")
    }
  }
  colnames(info_SHHK)<-c("knot1","knot2","family","tau")
  info_SHHK  
}

#绘图
#aic、bic、lik
plotaic<-function(roll_info){
  plot(roll_info[,1],type="l",col = "red", xlab = "Time", ylab = "Value", 
       main = "AIC & BIC & Lik")
  lines(roll_info[,2], type = "l", col = "blue")
  lines(roll_info[,3], type = "l", col = "green")
  legend("topright", c("aic", "bic","lik"),pch = c(1,1,1),col=c("red","blue","green"),bg ="white")
}

TimeDiff <- function(start_time) {
 start_time <- as.POSIXct(start_time)
 dt <- difftime(Sys.time(), start_time, units="secs")
 format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
