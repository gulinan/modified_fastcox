library(fastcox)#
library(glmnet)#
data(FHT)#
x=FHT$x#
y=FHT$y#
d=FHT$status#
wy = cbind(time=y,status=d)
setwd('/Users/emeryyi/Desktop')
nfolds=10#
N=nrow(x)#
###Fit the model once to get dimensions etc of output#
cocktail.object=cocktail(x=x,y=y,d=d)#
# cocktail.object=glmnet(x=x,y=wy,family="cox")#
lambda=cocktail.object$lambda#
nz=sapply(predict(cocktail.object,type="nonzero"),length)#
foldid=sample(rep(seq(nfolds),length=N))#
if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=10 recommended")#
 outlist=as.list(seq(nfolds))#
###Now fit the nfold models and store them#
for(i in seq(nfolds)){#
  which=foldid==i#
  if(is.matrix(y))y_sub=y[!which,]else y_sub=y[!which]#
if(is.matrix(d))d_sub=d[!which,]else d_sub=d[!which]#
  outlist[[i]]=cocktail(x[!which,,drop=FALSE],y_sub,d_sub)#
}
nfolds=max(foldid)#
  nobs=as.integer(length(y))#
  weights=rep(1.0,nobs)#
  offset=rep(0.0,nobs)#
  cvraw=matrix(NA,nfolds,length(lambda))#
  for(i in seq(nfolds)){#
    which=foldid==i#
    fitobj=outlist[[i]]#
    coefmat=predict(fitobj,type="coeff")#
    plfull=survpath.deviance(x=x,y=y,d=d,beta=coefmat)#
    plminusk=survpath.deviance(x=x[!which,],y=y[!which],d=d[!which],beta=coefmat)#
    coxy = cbind(time=y,status=d)#
	plfull1=coxnet.deviance(x=x,y=y,d=d,beta=coefmat)#
    plminusk1=coxnet.deviance(x=x[!which,],coxy=y[!which,],beta=coefmat)#
    print(plfull-plfull1)#
	print(plminusk-plminusk1)#
	cvraw[i,seq(along=plfull)]=plfull-plminusk#
  }
nfolds=max(foldid)#
  nobs=as.integer(length(y))#
  weights=rep(1.0,nobs)#
  offset=rep(0.0,nobs)#
  cvraw=matrix(NA,nfolds,length(lambda))#
  for(i in seq(nfolds)){#
    which=foldid==i#
    fitobj=outlist[[i]]#
    coefmat=predict(fitobj,type="coeff")#
    plfull=survpath.deviance(x=x,y=y,d=d,beta=coefmat)#
    plminusk=survpath.deviance(x=x[!which,],y=y[!which],d=d[!which],beta=coefmat)#
    coxy = cbind(time=y,status=d)#
	plfull1=coxnet.deviance(x=x,y=y,beta=coefmat)#
    plminusk1=coxnet.deviance(x=x[!which,],coxy=y[!which,],beta=coefmat)#
    print(plfull-plfull1)#
	print(plminusk-plminusk1)#
	cvraw[i,seq(along=plfull)]=plfull-plminusk#
  }
nfolds=max(foldid)#
  nobs=as.integer(length(y))#
  weights=rep(1.0,nobs)#
  offset=rep(0.0,nobs)#
  cvraw=matrix(NA,nfolds,length(lambda))#
  for(i in seq(nfolds)){#
    which=foldid==i#
    fitobj=outlist[[i]]#
    coefmat=predict(fitobj,type="coeff")#
    plfull=survpath.deviance(x=x,y=y,d=d,beta=coefmat)#
    plminusk=survpath.deviance(x=x[!which,],y=y[!which],d=d[!which],beta=coefmat)#
    coxy = cbind(time=y,status=d)#
	plfull1=coxnet.deviance(x=x,y=coxy,beta=coefmat)#
    plminusk1=coxnet.deviance(x=x[!which,],coxy=y[!which,],beta=coefmat)#
    print(plfull-plfull1)#
	print(plminusk-plminusk1)#
	cvraw[i,seq(along=plfull)]=plfull-plminusk#
  }
for(i in seq(nfolds)){#
    which=foldid==i#
    fitobj=outlist[[i]]#
    coefmat=predict(fitobj,type="coeff")#
    plfull=survpath.deviance(x=x,y=y,d=d,beta=coefmat)#
    plminusk=survpath.deviance(x=x[!which,],y=y[!which],d=d[!which],beta=coefmat)#
    coxy = cbind(time=y,status=d)#
	plfull1=coxnet.deviance(x=x,y=coxy,beta=coefmat)#
    plminusk1=coxnet.deviance(x=x[!which,],y=coxy[!which,],beta=coefmat)#
    print(plfull-plfull1)#
	print(plminusk-plminusk1)#
	cvraw[i,seq(along=plfull)]=plfull-plminusk#
  }
status=d
N=nfolds - apply(is.na(cvraw),2,sum)
N
nfolds
is.na(cvraw),2
dim(cvraw)
tapply(weights*status,foldid,sum)
weights=as.vector(tapply(weights*status,foldid,sum))
status
weights=as.vector(tapply(status,foldid,sum))
cvraw=cvraw/weights
cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
cvm
weights
cvraw
setwd('/Users/emeryyi/Dropbox/Research/googleproject/fastcox/R')
cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
