################ function:count ###################
countl<-function(center,nn,kk,rsq,a,v,leqcalc=T){
#                                      T T
# Counts points in an elipse given by x a a x,
#  browser()
   eps<-if(leqcalc) 1.0e-8 else -1.0e-8
   if(nn<kk){
      cat("kk=",kk," nn=",nn," kk should never exceed nn\n")
   }
   if(is.na(rsq)) rsq<- -1
   if(rsq<0){
      count<-NA
   }else{
      count<-0
      oc<-center[kk-1]
      #cat("center[kk]= ", center[kk], "rsq = ", rsq, "a[kk, kk] =", a[kk, kk], "leqcalc=", leqcalc)
      b<-round(.5+center[kk]-sqrt(rsq)/a[kk,kk]-eps)
      e<-round(-.5+center[kk]+sqrt(rsq)/a[kk,kk]+eps)
     #cat("kk=",kk,"b=",b,"e=",e,"\n")
#     browser()
      if(e>=b){
         for(j in b:e){
            nrsq<-rsq-a[kk,kk]**2*(j-center[kk])**2
            v[kk]<-j      
            if(kk>1){
               center[kk-1]<-oc-a[kk-1,kk]*(j-center[kk])/a[kk-1,kk-1]
               count<-count+countl(center,nn,kk-1,nrsq, a,v)
            }else{
               v[nn+1]<-nrsq
               count<-1+count
            }
         }
      }
      center[kk-1]<-oc
   }
   return(count)
}
