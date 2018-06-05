scia <-
function(matr, cla, nam, cor=FALSE, net=HPRD, set1=0.02, set2=0.02,pvalue=0.05,perm=200){
    library(igraph)
    library("ICSNP")
    a=unique(cla)[1]
    b=unique(cla)[2]
    b1=apply(matr[,which(a==cla)],1,function(I){length(unique(I))==1})
    b2=apply(matr[,which(b==cla)],1,function(I){length(unique(I))==1})
    b=unique(c(which(b1),which(b2)))
    if(length(b)>0){
        matr=matr[-b,]
        nam=nam[-b]
    }
    a=unique(c(net[,1],net[,2]))
    name=intersect(a,nam)
    b=setdiff(a,name)
    if(length(b)>0){
        b1=sapply(net[,1],function(I){any(I==b)})
        b2=sapply(net[,2],function(I){any(I==b)})
        b1=unique(c(which(b1),which(b2)))
        net=net[-b1,]
    }
    name=unique(c(net[,1],net[,2]))

    X=apply(as.matrix(name),1,function(I){ matr[which(nam==I),] })
    a=unique(cla)[1]
    b=unique(cla)[2]
    N1=sum(cla==a)
    N2=sum(cla==b)
    A=X[which(cla==a),]
    B=X[which(cla==b),]
    r=sapply(c(1:dim(X)[2]),function(I){t.test(A[,I],B[,I])[[3]]})
    x=set1
    y=set2
    if(class(set1)=="numeric"){
        rr1=which(r<quantile(r,set1))
        x=name[rr1]
    }
    if(class(set2)=="numeric"){
        rr2=which(r<quantile(r,set2))
        y=name[rr2]
    }

    G=graph.data.frame(as.data.frame(net), directed=F)
    G=set.vertex.attribute(G, "exp", V(G), c(1:dim(X)[2]))
    x.ig=apply(as.matrix(x),1,function(i){which(i==name)})
    y.ig=apply(as.matrix(y),1,function(i){which(i==name)})
    a.ig=apply(as.matrix(x.ig),1,function(i){
        get.all.shortest.paths(G,i,y.ig)
    })
    a=lapply(a.ig,function(i){
        x=sapply(i[[1]],function(j){
            tail(j,1L)
        })
        y=lapply(as.list(y.ig),function(k){
            if(any(k==x)){return(which(k==x))}
            else         {return(NA)}
        })
        z=lapply(y,function(l){
            i[[1]][l]
        })
    })

    Z=lapply(a,function(i){
        lapply(i,function(j){
            x=lapply(j,function(k){
                if(length(k)<2){return(NA);next}
                e=get.vertex.attribute(G,"exp",k)
                y=HotellingsT2(A[,e],B[,e])
                if(cor==FALSE){return(y[[1]][1]);next}else{
                    n=length(e)
                    r1=y[[1]][1]*(N1+N2-2)*n/(N1+N2-1-n)
                    B1=sapply(1:(n-1),function(I){summary(lm(A[,e[I+1]]~A[,e[I]]))[[4]][2,1:2]})
                    BB1=sapply(1:(n-1),function(I){summary(lm(B[,e[I+1]]~B[,e[I]]))[[4]][2,1:2]})
                    r2=sapply(1:(n-1),function(I){(B1[1,I]-BB1[1,I])^2/(B1[2,I]^2+BB1[2,I]^2)})
                    r2=sum(r2)
                    return(sum(r1,r2))
                }
            })

            if(cor==FALSE){xx=which.min(x)}else{xx=which.max(x)}
            if(length(xx)==1){return(j[[xx]])}
            if(length(xx)==0){return(NA)}
        })
    })



    S=lapply(Z,function(I){
        lapply(I,function(J){
            if(is.na(J[1])){return(NA);next}
            e=get.vertex.attribute(G,"exp",J)
            y=HotellingsT2(A[,e],B[,e])
            if(cor==FALSE){return(y[[1]][1]);next}else{
                n=length(e)
                r1=y[[1]][1]*(N1+N2-2)*n/(N1+N2-1-n)
                B1=sapply(1:(n-1),function(I){summary(lm(A[,e[I+1]]~A[,e[I]]))[[4]][2,1:2]})
                BB1=sapply(1:(n-1),function(I){summary(lm(B[,e[I+1]]~B[,e[I]]))[[4]][2,1:2]})
                r2=sapply(1:(n-1),function(I){(B1[1,I]-BB1[1,I])^2/(B1[2,I]^2+BB1[2,I]^2)})
                r2=sum(r2)
                return(sum(r1,r2))
            }
        })
    })



    G=graph.data.frame(as.data.frame(net), directed=F)
    G=set.vertex.attribute(G, "exp", V(G), c(1:dim(X)[2]))
    y.ig=apply(as.matrix(y),1,function(i){which(i==name)})
    r=apply(as.matrix(c(1:perm)),1,function(I){
        sample(1:length(name),length(y),F)
    })
    R=cbind(y.ig,r)
    a.ig=apply(R,1,function(I){
        get.all.shortest.paths(G,I[1],I[-1])
    })

    a=apply(as.matrix(1:length(a.ig)),1,function(I){
        x=sapply(a.ig[[I]][[1]],function(J){
            tail(J,1L)
        })
        y=sapply(R[I,-1],function(k){
            if(any(k==x)){return(which(k==x))}
            else         {return(NA)}
        })
        z=sapply(y,function(l){
            a.ig[[I]][[1]][l]
        })
    })
    if(length(y.ig)==1){a=list(a)}

    Z.r=lapply(a,function(i){
        lapply(i,function(j){
            x=lapply(j,function(k){
                if(length(k)<2){return(NA);next}
                e=get.vertex.attribute(G,"exp",k)
                y=HotellingsT2(A[,e],B[,e])
                if(cor==FALSE){return(y[[1]][1]);next}else{
                    n=length(e)
                    r1=y[[1]][1]*(N1+N2-2)*n/(N1+N2-1-n)
                    B1=sapply(1:(n-1),function(I){summary(lm(A[,e[I+1]]~A[,e[I]]))[[4]][2,1:2]})
                    BB1=sapply(1:(n-1),function(I){summary(lm(B[,e[I+1]]~B[,e[I]]))[[4]][2,1:2]})
                    r2=sapply(1:(n-1),function(I){(B1[1,I]-BB1[1,I])^2/(B1[2,I]^2+BB1[2,I]^2)})
                    r2=sum(r2)
                    return(sum(r1,r2))
                }
            })

            if(cor==FALSE){xx=which.min(x)}else{xx=which.max(x)}
            if(length(xx)==1){return(j[[xx]])}
            if(length(xx)==0){return(NA)}
        })
    })



    S.r=lapply(Z.r,function(I){
        lapply(I,function(J){
            if(is.na(J[1])){return(NA);next}
            e=get.vertex.attribute(G,"exp",J)
            y=HotellingsT2(A[,e],B[,e])
            if(cor==FALSE){return(y[[1]][1]);next}else{
                n=length(e)
                r1=y[[1]][1]*(N1+N2-2)*n/(N1+N2-1-n)
                B1=sapply(1:(n-1),function(I){summary(lm(A[,e[I+1]]~A[,e[I]]))[[4]][2,1:2]})
                BB1=sapply(1:(n-1),function(I){summary(lm(B[,e[I+1]]~B[,e[I]]))[[4]][2,1:2]})
                r2=sapply(1:(n-1),function(I){(B1[1,I]-BB1[1,I])^2/(B1[2,I]^2+BB1[2,I]^2)})
                r2=sum(r2)
                return(sum(r1,r2))
            }
        })
    })

    p=apply(as.matrix(c(1:length(x))),1,function(I){
        apply(as.matrix(c(1:length(y))),1,function(J){
            if(is.na(S[[I]][[J]])){(return(NA))}
            else{return(length(which(S.r[[J]]<S[[I]][[J]]))/sum(!is.na(S.r[[J]])))}
        })
    })

    rr=which(p<pvalue,arr.ind = TRUE)
    if(is.null(dim(p))){rr=cbind(1,rr)}
    a=sapply(Z,function(I){
        sapply(I,function(J){
            length(J)
        })
    })
    if(is.null(dim(a))){a=t(as.matrix(a))}

    aa=sapply(1:length(Z),function(I){
        sapply(1:length(Z[[1]]),function(J){
            1-pchisq(S[[I]][[J]],a[J,I])
        })
    })
    aa=matrix(p.adjust(aa,"BH"),length(Z[[1]]),length(Z))
    rr2=which(aa<0.05,arr.ind = TRUE)
    rr1=apply(rr,1,function(I){sort(c(x[I[2]],y[I[1]]))})
    rr2=apply(rr2,1,function(I){sort(c(x[I[2]],y[I[1]]))})
    if(length(rr1)<2|length(rr2)==0){stop("no significant result")}

    rr3=apply(rr1,2,function(I){
        apply(rr2,2,function(J){
            sum(I==J)==2
        })
    })
    rr4=apply(rr3,2,function(I){any(I)})
    rr5=rr[which(rr4),]
    if(is.null(dim(rr5))){rr5=t(as.matrix(rr5))}

    RR=apply(rr5,1,function(I){list(Z[[I[2]]][[I[1]]])})
    RRR=lapply(RR,function(I){
        x=name[I[[1]]]
        y=cbind(x[-length(x)],x[-1])
    })
    R=RRR[[1]]
    if(length(RRR)!=1){
        R=RRR[[1]]
        for(i in 2:length(RRR)){R=rbind(R,RRR[[i]])}
    R
    print(R)
    }
}
