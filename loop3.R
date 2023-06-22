partstat<-function(nv,groups,ranks,center=TRUE){
  stat<-array(NA,length(nv))
  for(ii in seq(length(nv))){
    stat[ii]<-sum(ranks*(groups==ii))
    if(center) stat[ii]<-stat[ii]-nv[ii]*(sum(nv)+1)/2
  }
  return(stat)
}
loopme2<-function(nv){
  k3<-array(0,rep(length(nv),3))
  u4<-k4<-array(0,rep(length(nv),4))
  groups<-rep(seq(length(nv)),nv)
  dg<-diff(groups)
  #  cat("groups",groups,"dg",dg,"\n")
  ranks<-seq(sum(nv))
  count<-1
  stat<-partstat(nv,groups,ranks,TRUE)
  ustat<-partstat(nv,groups,ranks,FALSE)
  for(ii in 1:length(nv)){
    for(jj in 1:length(nv)){
      for(kk in 1:length(nv)){
        k3[ii,jj,kk]<-k3[ii,jj,kk]+ stat[ii]* stat[jj]* stat[kk]
        for(ll in 1:length(nv)){
          u4[ii,jj,kk,ll]<-u4[ii,jj,kk,ll]+ustat[ii]*ustat[jj]*ustat[kk]*ustat[ll]
          k4[ii,jj,kk,ll]<-k4[ii,jj,kk,ll]+ stat[ii]* stat[jj]* stat[kk]*stat[ll]
        }
      }
    }
  }
  ii<-0
  while(any(dg>0)){
    ii<-max(seq(length(dg))[dg>0]); # cat("ii",ii,"\n")
    tg<-groups[(ii+1):length(groups)]>groups[ii]
    newg<-min(groups[(ii+1):length(groups)][tg])
    kk<-min(((ii+1):length(groups))[groups[(ii+1):length(groups)]==newg]); # cat("kk",kk,"\n")
    #     cat("pre swap groups",groups,"\n")
    temp<-groups[ii]
    groups[ii]<-groups[kk]
    groups[kk]<-temp
    #     cat("pre sort groups",groups,"ii",ii,"\n")
    if(ii<=(length(groups)-2)){
      groups[(ii+1):length(groups)]<-sort(groups[(ii+1):length(groups)])
    }
    dg<-diff(groups)
    count<-count+1
    stat<-partstat(nv,groups,ranks,TRUE)
    ustat<-partstat(nv,groups,ranks,FALSE)
    #     cat("count",count,"groups",groups,"dg",dg,"ustat",ustat,"\n")
    for(ii in 1:length(nv)){
      for(jj in 1:length(nv)){
        for(kk in 1:length(nv)){
          k3[ii,jj,kk]<-k3[ii,jj,kk]+ stat[ii]* stat[jj]* stat[kk]
          for(ll in 1:length(nv)){
            u4[ii,jj,kk,ll]<-u4[ii,jj,kk,ll]+ustat[ii]*ustat[jj]*ustat[kk]*ustat[ll]
            k4[ii,jj,kk,ll]<-k4[ii,jj,kk,ll]+ stat[ii]* stat[jj]* stat[kk]* stat[ll]
          }
        }
      }
    }
    #     cat("ustat[1]*ustat[1]*ustat[1]*ustat[1]",ustat[1]*ustat[1]*ustat[1]*ustat[1],"u4[1,1,1,1]",u4[1,1,1,1],"\n")
    if(count>100000000) dg<--1
  }
  return(list(k1=nv*(sum(nv)+1)/2,k3=k3/count,k4=k4/count,u4=u4/count,flag=(length(dg)==1)))
}
sfcn<-function(n,hard=FALSE,center=FALSE){
  if(hard){
    S<-rep(0,5)
    if(center){
      for(ii in 1:n){
        S[1]<-S[1]+(ii-(n+1)/2)^4
      }
      for(ii in 1:n) for(jj in (1:n)){
        if(ii!=jj) S[2]<-S[2]+(ii-(n+1)/2)^3*(jj-(n+1)/2)
      }
      for(ii in 1:n) for(jj in (1:n)) {
        if(ii!=jj) S[3]<-S[3]+(ii-(n+1)/2)^2*(jj-(n+1)/2)^2
      }
      for(ii in 1:n) for(jj in (1:n)) for(kk in (1:n)){
        if((ii!=jj)&(ii!=kk)&(jj!=kk)) S[4]<-S[4]+(ii-(n+1)/2)^2*(jj-(n+1)/2)*(kk-(n+1)/2);
      }
      for(ii in 1:n) for(jj in (1:n)) for(kk in (1:n)) for(mm in (1:n)){
        if( (ii!=jj)&(ii!=kk)&(jj!=kk)& (mm!=jj)&(ii!=mm)&(mm!=kk)) S[5]<-S[5]+(ii-(n+1)/2)*(jj-(n+1)/2)*(kk-(n+1)/2)*(mm-(n+1)/2)
      }
      
    }else{
      for(ii in 1:n){
        S[1]<-S[1]+ii^4
      }
      for(ii in 1:n) for(jj in (1:n)){
        if(ii!=jj) S[2]<-S[2]+ii^3*jj
      }
      for(ii in 1:n) for(jj in (1:n)) {
        if(ii!=jj) S[3]<-S[3]+ii^2*jj^2
      }
      for(ii in 1:n) for(jj in (1:n)) for(kk in (1:n)){
        if((ii!=jj)&(ii!=kk)&(jj!=kk)) S[4]<-S[4]+ii^2*jj*kk;
      }
      for(ii in 1:n) for(jj in (1:n)) for(kk in (1:n)) for(mm in (1:n)){
        if( (ii!=jj)&(ii!=kk)&(jj!=kk)& (mm!=jj)&(ii!=mm)&(mm!=kk)) S[5]<-S[5]+ii*jj*kk*mm
      }
    }
  }
  else{
    if(center){
      S<-c((n*(7 - 10*n^2 + 3*n^4))/240.,
           (-7*n + 10*n^3- 3*n^5)/240.,
           (n*(-21 + 5*n + 30*n^2- 10*n^3- 9*n^4+ 5*n^5))/720.,
           (n*(42 - 5*n - 60*n^2+ 10*n^3+ 18*n^4- 5*n^5))/720.,
           (n*(-42 + 5*n + 60*n^2- 10*n^3- 18*n^4+ 5*n^5))/240.)
    }
    else{
      S<-c((n*(1 + n)*(1 + 2*n)*(-1 + 3*n + 3*n^2))/30,
           (n*(1 + n)*(4 - 4*n - 21*n^2 + 6*n^3 + 15*n^4 ))/120,
           (n*(1 + n)*(1 + 2*n)*(6 - 13*n - 3*n^2 + 10*n^3 ))/180,
           (n*(1 + n)*(-24 + 14*n + 91*n^2 - 56*n^3 - 55*n^4 + 30*n^5 ))/360.,
           (n*(1 + n)*(48 - 28*n - 152*n^2 + 127*n^3 + 65*n^4 - 75*n^5 + 15*n^6 ))/240.)
    }
  }
  return(S)
}
truetrue<-function(nv,center=FALSE){
  n<-sum(nv)
  Power<-function(n,j) n^j
  S<-sfcn(n,center=center)
  #cat("S[1]",S[1],"\n")
  #cat("S[2]",S[2],"\n")
  #cat("S[3]",S[3],"\n")
  #cat("S[4]",S[4],"\n")
  #cat("S[5]",S[5],"\n")
  d4<-array(NA,rep(length(nv),4));d2<-array(NA,rep(length(nv),2))
  for(aa in seq(length(nv))){
    d2[aa,aa]<- (sum(nv)-nv[aa])*nv[aa]*(sum(nv)+1)/12
    d4[aa,aa,aa,aa]<-S[1]*nv[aa]/n +4*S[2]*nv[aa]*(nv[aa]-1)/(n*(n-1)) +3*S[3]*nv[aa]*(nv[aa]-1)/(n*(n-1)) +6*S[4]*nv[aa]*(nv[aa]-1)*(nv[aa]-2)/(n*(n-1)*(n-2)) +S[5]*nv[aa]*(nv[aa]-1)*(nv[aa]-2)*(nv[aa]-3)/(n*(n-1)*(n-2)*(n-3))
    for(bb in seq(length(nv))){
      if(bb!=aa){
        d2[aa,bb]<- -nv[aa]*nv[bb]*(sum(nv)+1)/12
        d4[aa,bb,bb,bb]<-S[1]*0+S[2]*nv[aa]*nv[bb]/(n*(n-1))+S[3]*0+S[4]*3*nv[aa]*nv[bb]*(nv[bb]-1)/(n*(n-1)*(n-2)) +S[5]*nv[aa]*nv[bb]*(nv[bb]-1)*(nv[bb]-2)/(n*(n-1)*(n-2)*(n-3))
        d4[bb,aa,bb,bb]<-d4[aa,bb,bb,bb]
        d4[bb,bb,aa,bb]<-d4[aa,bb,bb,bb] 
        d4[bb,bb,bb,aa]<-d4[aa,bb,bb,bb]
        d4[aa,aa,bb,bb]<-S[1]*0+S[2]*0+S[3]*nv[aa]*nv[bb]/(n*(n-1))+S[4]*nv[aa]*nv[bb]*(nv[aa]+nv[bb]-2)/(n*(n-1)*(n-2)) +S[5]*nv[aa]*(nv[aa]-1)*nv[bb]*(nv[bb]-1)/(n*(n-1)*(n-2)*(n-3))
        d4[aa,bb,bb,aa]<-d4[aa,aa,bb,bb]
        d4[aa,bb,aa,bb]<-d4[aa,aa,bb,bb]
        for(cc in seq(length(nv))){
          if((cc!=bb)&(cc!=aa)){
            d4[aa,aa,bb,cc]<-S[1]*0+S[2]*0+S[3]*0+S[4]*nv[aa]*nv[bb]*nv[cc]/(n*(n-1)*(n-2)) +S[5]*nv[aa]*(nv[cc]-1)*nv[bb]*nv[cc]/(n*(n-1)*(n-2)*(n-3))
            d4[aa,bb,aa,cc]<- d4[aa,aa,bb,cc]
            d4[bb,aa,aa,cc]<- d4[aa,aa,bb,cc]
            d4[bb,cc,aa,aa]<- d4[aa,aa,bb,cc]
            d4[bb,aa,cc,aa]<- d4[aa,aa,bb,cc]
            d4[aa,bb,cc,aa]<- d4[aa,aa,bb,cc]
            for(dd in seq(length(nv))){
              #                    if((aa==1)&(bb==2)&(cc==3)) cat("dd=",dd)
              if((dd!=aa)& (dd!=bb)& (dd!=cc)) d4[aa,bb,cc,dd]<-S[1]*0+S[2]*0+S[3]*0+S[4]*0+S[5]*nv[aa]*nv[bb]*nv[cc]*nv[dd]/(n*(n-1)*(n-2)*(n-3))
            }
          }
        }
      }
    }
  }
  return(list(d2=d2,d4=d4))
}
central_moments_calc <- function(nv) {
  len <- length(nv)
  m2m4 <-truetrue(nv, center=TRUE)
  m1 <- array(0,rep(len,1))
  m2 <- m2m4$d2
  m3 <- array(0,rep(len,3))
  m4 <- m2m4$d4
  return(list(m1=m1,m2=m2,m3=m3,m4=m4))
}


# Testing...



