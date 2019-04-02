library(glmnet)
library(ncvreg)
library(ggplot2)
dat1=read.csv("GSE14814_all_clin.csv")
dat2=read.csv("GSE14814_all_ex.csv",header=TRUE)
rownames(dat2)=dat2[,1]
x=t(dat2[,-1])
y=dat1[,c(-2,-3)]

#a=which(y[,2]=="ACT)

y1=y[,8]-mean(y[,8])
n=length(y1)
x1=scale(x)*sqrt(n/(n-1))
fit1=glmnet(x1,y1,alpha=1)

####solution path####
fit2=ncvreg(x1,y1,lambda=fit1$lambda,penalty='lasso')
fit2$beta=fit2$beta[-1,]
s1=c()
for(i in 1:22215){
  s1=c(s1,sum(fit1$beta[i,]))
}
ip1=which(s1!=0)

p1=rep(rownames(fit2$beta[ip1,]),each=length(fit1$lambda))
lambda1=rep(fit2$lambda,length(ip1))
beta1=c(t(fit2$beta[ip1,]))
df1=data.frame(p1,lambda1,beta1)

ggplot(data=df1, aes(x=lambda1, y=beta1, group=p1)) +
  ggtitle("Solution Path")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(color=p1),size=0.7)+
  xlab("Lambda")+
  ylab("Coefficient")+
  theme(legend.position="none")
###################

cv=cv.glmnet(x1,y1,lambda=fit1$lambda,nfold=n,alpha=1)

##AIC BIC Leave-One-Out##
AB=function(fit,method){
  tLL=fit$nulldev-deviance(fit)
  k=fit$df
  n=fit$nobs
  if(method=="AIC"){
    mm=-tLL+2*k
  }else{
    mm=log(n)*k - tLL
  }
  return(mm)
}

aic=AB(fit1,"AIC")
bic=AB(fit1,"BIC")
k1=which(aic==min(aic))
k2=which(bic==min(bic))
p1=(which(coef(fit1,s=fit1$lambda[k1])!=0)-1)[-1]
p2=(which(coef(fit1,s=fit1$lambda[k2])!=0)-1)[-1]
p3=(which(coef(cv,s=cv$lambda.min)!=0)-1)[-1]

ggplot(data=data.frame(1,fit1$lambda,aic),aes(x=fit1.lambda,y=aic))+
  geom_line(color="limegreen",size=1)+
  geom_point(data =data.frame(1,fit1$lambda[k1],aic[k1]),
             mapping = aes(x = fit1.lambda.k1., y = aic.k1.))+
  geom_text(data =data.frame(1,fit1$lambda[k1],aic[k1]),
            aes(x = fit1.lambda.k1., y = aic.k1.,label=paste("Optimal lambda = ",round(fit1$lambda[k1],5))),hjust=-0.3,vjust=0)+
  ggtitle("AIC")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Lambda")+
  ylab("AIC")

ggplot(data=data.frame(1,fit1$lambda,bic),aes(x=fit1.lambda,y=bic))+
  geom_line(color="limegreen",size=1)+
  geom_point(data =data.frame(1,fit1$lambda[k2],bic[k2]),
             mapping = aes(x = fit1.lambda.k2., y = bic.k2.))+
  geom_text(data =data.frame(1,fit1$lambda[k2],bic[k2]),
            aes(x = fit1.lambda.k2., y = bic.k2.,label=paste("Optimal lambda = ",round(fit1$lambda[k2],5))),hjust=-0.3,vjust=0)+
  ggtitle("BIC")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Lambda")+
  ylab("BIC")

ggplot(data=data.frame(1,cv$lambda,cv$cvm),aes(x=cv.lambda,y=cv.cvm))+
  geom_line(color="limegreen",size=1)+
  geom_point(data =data.frame(1,cv$lambda.min,min(cv$cvm)),
             mapping = aes(x = cv.lambda.min, y = min.cv.cvm.))+
  geom_text(data =data.frame(1,cv$lambda.min,min(cv$cvm)),
            aes(x = cv.lambda.min, y = min.cv.cvm.,label=paste("Optimal lambda = ",round(cv$lambda.min,5))),hjust=-0.3,vjust=0)+
  ggtitle("Leave-One-Out")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Lambda")+
  ylab("CV Error")

######################tree#####################
library(party)
tree=c()
for(i in 100){
  train=sample(1:n,100)
  time=ifelse(y$OS.time>=6,"Over","Less")
  gene=data.frame(time,y$Age,as.numeric(y$Sex),as.numeric(y$Stage),as.numeric(y$Post.Surgical.Treatment),as.numeric(y$Histology.type),x[,p2])
  ct=ctree(time~ ., data = gene,subset=train)
  plot(ct, main = "æ¢ä»¶?Ž¨è«–æ¨¹")
  pred.time= predict(ct,newdata=gene[-train,-1])
  table(time[-train],pred.time)
  tree=c(tree,mean(time[-train]==pred.time))
}
mean(tree)



##randomForest##
library(randomForest)
accuracy=c()
for(i in 1:100){
  train=sample(1:n,100)
  set.seed(NULL)
  rf.gene=randomForest(time~.,data=gene,subset=train,importance=TRUE,ntree=500)
  #plot(rf.gene)
  pred.gene=predict(rf.gene,newdata=gene[-train,-1])
  gene.test=gene[-train,"time"]
  #varImpPlot(rf.gene)
  #table(gene.test,pred.gene)
  accuracy=c(accuracy,mean(gene.test==pred.gene))
}

mean(accuracy)

###########################

##########ACT OBS##########
a=which(y[,2]=="ACT")
n=length(a)
y1=y[a,8]-mean(y[a,8])
x1=scale(x[a,])*sqrt(n/(n-1))
fit1=glmnet(x1,y1,alpha=1)

####solution path####
fit2=ncvreg(x1,y1,lambda=fit1$lambda,penalty='lasso')
fit2$beta=fit2$beta[-1,]
s1=c()
for(i in 1:22215){
  s1=c(s1,sum(fit1$beta[i,]))
}
ip1=which(s1!=0)

p1=rep(rownames(fit2$beta[ip1,]),each=length(fit1$lambda))
lambda1=rep(fit2$lambda,length(ip1))
beta1=c(t(fit2$beta[ip1,]))
df1=data.frame(p1,lambda1,beta1)

ggplot(data=df1, aes(x=lambda1, y=beta1, group=p1)) +
  ggtitle("Solution Path (ACT)")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(color=p1),size=0.7)+
  xlab("Lambda")+
  ylab("Coefficient")+
  theme(legend.position="none")

###################

cv=cv.glmnet(x1,y1,lambda=fit1$lambda,nfold=n,alpha=1)

##AIC BIC Leave-One-Out##
AB=function(fit,method){
  tLL=fit$nulldev-deviance(fit)
  k=fit$df
  n=fit$nobs
  if(method=="AIC"){
    mm=-tLL+2*k
  }else{
    mm=log(n)*k - tLL
  }
  return(mm)
}

aic=AB(fit1,"AIC")
bic=AB(fit1,"BIC")
k1=which(aic==min(aic))
k2=which(bic==min(bic))
p1=(which(coef(fit1,s=fit1$lambda[k1])!=0)-1)[-1]
p2=(which(coef(fit1,s=fit1$lambda[k2])!=0)-1)[-1]
p3=(which(coef(cv,s=cv$lambda.min)!=0)-1)[-1]

ggplot(data=data.frame(1,fit1$lambda,aic),aes(x=fit1.lambda,y=aic))+
  geom_line(color="limegreen",size=1)+
  geom_point(data =data.frame(1,fit1$lambda[k1],aic[k1]),
             mapping = aes(x = fit1.lambda.k1., y = aic.k1.))+
  geom_text(data =data.frame(1,fit1$lambda[k1],aic[k1]),
            aes(x = fit1.lambda.k1., y = aic.k1.,label=paste("Optimal lambda = ",round(fit1$lambda[k1],5))),hjust=-0.3,vjust=0)+
  ggtitle("AIC")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Lambda")+
  ylab("AIC")

ggplot(data=data.frame(1,fit1$lambda,bic),aes(x=fit1.lambda,y=bic))+
  geom_line(color="limegreen",size=1)+
  geom_point(data =data.frame(1,fit1$lambda[k2],bic[k2]),
             mapping = aes(x = fit1.lambda.k2., y = bic.k2.))+
  geom_text(data =data.frame(1,fit1$lambda[k2],bic[k2]),
            aes(x = fit1.lambda.k2., y = bic.k2.,label=paste("Optimal lambda = ",round(fit1$lambda[k2],5))),hjust=-0.3,vjust=0)+
  ggtitle("BIC")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Lambda")+
  ylab("BIC")

ggplot(data=data.frame(1,cv$lambda,cv$cvm),aes(x=cv.lambda,y=cv.cvm))+
  geom_line(color="limegreen",size=1)+
  geom_point(data =data.frame(1,cv$lambda.min,min(cv$cvm)),
             mapping = aes(x = cv.lambda.min, y = min.cv.cvm.))+
  geom_text(data =data.frame(1,cv$lambda.min,min(cv$cvm)),
            aes(x = cv.lambda.min, y = min.cv.cvm.,label=paste("Optimal lambda = ",round(cv$lambda.min,5))),hjust=-0.3,vjust=0)+
  ggtitle("Leave-One-Out")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Lambda")+
  ylab("CV Error")

act=colnames(x[,unique(c(p1,p2,p3))])

ss=c()
for(i in 1:length(act)){
  ss=c(ss,sum(act[i]==obs))
}

##############################
