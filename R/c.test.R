c.test <-
function(matr,cla,cor=FALSE,net=NULL,nam=NULL){
      library("ICSNP")
      a=unique(cla)[1]
      b=unique(cla)[2]
      a1=t(matr[,which(cla==a)])
      b1=t(matr[,which(cla==b)])
      N1=dim(a1)[1]
      N2=dim(b1)[1]
      n=dim(matr)[1]
      
      KK=HotellingsT2(a1,b1)
      if(cor==FALSE){
          r=c(KK[[1]][1],KK[[2]])
          names(r)=c("C.f","p-value")
          r
      }else{
          K1=KK[[1]][1]*(N1+N2-2)*n/(N1+N2-1-n)
          xx=apply(net,1,function(I){sum(any(I[1]==nam),any(I[2]==nam))==2})
          net=net[which(xx),]
          
          if(is.null(net)){
              B1=sapply(1:(n-1),function(I){
                  sapply((I+1):n,function(J){
                      summary(lm(a1[,I+1]~a1[,I]))[[4]][2,1:2]
                  })
                  
              })
              BB1=sapply(1:(n-1),function(I){
                  sapply((I+1):n,function(J){
                      summary(lm(b1[,I+1]~b1[,I]))[[4]][2,1:2]
                  })
                  
              })
              k1=sapply(1:((n+1)*n/2),function(I){(B1[1,I]-BB1[1,I])^2/(B1[2,I]^2+BB1[2,I]^2)})
              k1=sum(k1)
              r1=sum(k1,K1)
              r2=1-pchisq(r1,(n+3)*n/2)
              r=c(r1,r2)
              names(r)=c("C.x","p-value")
              r
          }else{
              nn=dim(net)[1]
              B1=sapply(1:nn,function(I){
                  x=c(which(nam==net[I,1]),which(nam==net[I,2]))
                  summary(lm(a1[,x[1]]~a1[,x[2]]))[[4]][2,1:2]
                  
              })
              BB1=sapply(1:nn,function(I){
                  x=c(which(nam==net[I,1]),which(nam==net[I,2]))
                  summary(lm(b1[,I+1]~b1[,I]))[[4]][2,1:2]
                  
              })
              k1=sapply(1:nn,function(I){(B1[1,I]-BB1[1,I])^2/(B1[2,I]^2+BB1[2,I]^2)})
              k1=sum(k1)
              r1=sum(k1,K1)
              r2=1-pchisq(r1,n+nn)
              r=c(r1,r2)
              names(r)=c("C.x","p-value")
              r
          }
          
          }
  }
