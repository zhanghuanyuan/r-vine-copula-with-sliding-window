"""输入
data           去除时间列
roll_long      区间长度 eg:200
roll_interval  间隔长度 eg:5
"""

"""输出
data_set      数据量     
data_number   时间长度   
roll_long     区间长度
interval_all  区间总数量
roll_interval 间隔长度
"""

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

interval_all<-(floor((data_number-roll_long)/roll_interval))+1

return(list(data_set=data_set,
    data_number=data_number,
    roll_interval=b,
    roll_long=a,
    interval_all=interval_all))
}


"""
描述性统计分析
"""
#判断
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

#名称
market_simp_name<-function(x){
    name<-switch(x,
        "SH",
        "NL",
        "AU",
        "AT",
        "BE",
        "US",
        "FR",
        "UK",
        "DE",
        "HK",
        "ID",
        "MY",
        "KR",
        "IT",
        "MX",
        "JP",
        "NZ",
        "NO",
        "PH",
        "RU",
        "IN",
        "CH",
        "SG",
        "TW")
}

#描述性统计分析
summarise<-function(data,data_set){
  summ<-matrix(data=NA,nrow=data_set+1,ncol=10)
  summ[1,]<-c("name","mean","max","min","sd","skewness","kurtosis","ks","adf","Arch-LM")
  
  for(i in 2:(data_set+1)){
    #name
    summ[i,1]<-i-1
    #data
    data_1<-as.numeric(c(unlist(data[,(i-1)])))
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



"""
标准化Copula数据
"""
#标准化的copula数据，用于计算 
##基于GARCH模型的标准化
copula_jisuan<-function(data){#arma(0,1)garch(1,1) as default
  fit<-garchFit(~arma(0,1)+garch(1,1),data,cond.dist = "std",trace=FALSE)
  bzh<-residuals(fit)/volatility(fit)
  bzh<-as.numeric(bzh)
  jisuan_copula<-ecdf(bzh)
}
##基于GJR-GARCH模型的标准化
copula_gjr<-function(data){#arma(0,1)gjrgarch(1,1) as default
  spec<-ugarchspec(#模型设置
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
  jisuan_gjr<-ecdf(bzh)
}

CopulaData<-function(data,data_set,data_number,t=1){#t=1,garch;t=2,gjr-garch
  data_copula<-matrix(data=NA,nrow=data_number,ncol=data_set)
  if(t==1){
    for(i in 1:data_set){#列
      data_copula[,i]<-copula_jisuan(data[,i])
    }
  }
  if(t==2){
    data_copula<-matrix(data=NA,nrow=data_number,ncol=data_set)
    for(i in 1:data_set){#列
      data_copula[,i]<-copula_gjr(data[,i])
    }
  }
  else{
    stop('t must be 1 or 2.')
  }
  data_copula
}

"""
R-Vine Copula相关计算
"""
#计算R-vine copula 
##MLE方法计算量巨大，谨慎使用
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

AIC_Loglike<-function(RVM,interval_all){
  roll_info<-matrix(data=NA,nrow=interval_all+1,ncol=3)
  roll_info[j+1,3]<-RVM$logLik
  roll_info[j+1,1]<-RVM$AIC
  roll_info[j+1,2]<-RVM$BIC
  colnames(roll_info)<-c("aic","bic","loglik")
  roll_info
}

#tree1矩阵提取
RVMtree1<-function(RVM_matrix){
  #both RVM$matrix and vine_matrix$matrix are available

  #输出竖排相关矩阵
  RVMtree1_col<-matrix(data=NA,nrow=(length(RVM_matrix[,1])-1),ncol=2)
  RVMtree1_col[,1]<-t(RVM_matrix[(length(RVM_matrix[,1])),1:(length(RVM_matrix[,1])-1)])
  for (i in 1:(length(RVM_matrix[,1])-1)){
    RVMtree1_col[i,2]<-RVM_matrix[i,i]
  }

  #输出社会网络型矩阵
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
"""
提取节点信息
"""
#提取treefamily、tau-value、par、par2
extract_t<-function(matrix){#
  extract_t<-t(matrix[(length(matrix[,1])),1:(length(matrix[,1])-1)])
}

Extract<-function(RVM,data_set){
  info<-matrix(data=NA,nrow=(data_set-1),ncol=11)
  info[,1:2]<-RVMtree1(RVM$Matrix)$RVMtree1_col  #knot
  info[,3]<-extract_t(RVM$family)                #family
  info[,5]<-extract_t(RVM$tau)                   #tau
  info[,6]<-extract_t(RVM$par)                   #par
  info[,7]<-extract_t(RVM$par2)                  #par2

  for(fnum in 1:length(info[,3])){
  family<-info[funm,3]
  info[fnum,4]<-switch(as.numeric(family),
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

  colnames(info)<-c("knot1","knot2","family","family name","tau","par","par2")
  info
}

"""
社会网络方法
"""
#点中心度计算
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

"""
CoVaR的计算
"""
#计算Copula-CoVaR
CoVaRCalculate<-function(data,knot,par,par2,type){#原始数据（去除时间列）与info的knot信息
    data_covar<-data[,knot]
    gjrfit<-ugarchfit(spec,data=data_covar,solver="solnp") #可以顺便提取边缘分布信息
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

CoVaRresult<-function(data,info,i){
  knot1<-as.numeric(info[i,1])
  knot2<-as.numeric(info[i,2])
  par<-as.numeric(info[i,6])
  par2<-as.numeric(info[i,7])
  type<-info[i,4]
  delta_2_1<-CoVaRCalculate(data,knot1,par,par2,type)
  delta_1_2<-CoVaRCalculate(data,knot2,par,par2,type)

  covar_name1<-market_simp_name(knot1)
  covar_name2<-market_simp_name(knot2)

  covar1_name<-paste(covar_name1,covar_name2,sep="_to_")
  covar2_name<-paste(covar_name2,covar_name1,sep="_to_")

  delta<-data.frame(delta_1_2,delta_2_1)
  colnames(delta)<-c(covar1_name,covar2_name)

  covar_name<-paste("CoVaR_",covar_name1,"_",covar_name2,".csv",sep="")
  comple<-assign(delta,covar_name)
  comple
  #write.csv(delta,file=covar_name)
}

deltacovar<-function(data,knot1,knot2,par,par2,type){
  delta_2_1<-CoVaRCalculate(data,knot1,par,par2,type)
  delta_1_2<-CoVaRCalculate(data,knot2,par,par2,type)
  
  covar_name1<-market_simp_name(knot1)
  covar_name2<-market_simp_name(knot2)
  
  covar2_name<-paste(covar_name1,covar_name2,sep="_to_")
  covar1_name<-paste(covar_name2,covar_name1,sep="_to_")
  
  delta<-data.frame(delta_1_2,delta_2_1)
  colnames(delta)<-c(covar1_name,covar2_name)
  list(delta,covar1_name,covar2_name)
}





