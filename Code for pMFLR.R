library(fdANOVA)
library(R.matlab)
library(fda)
library(MASS)  
library(fdapace)
library(fda.usc)
library(rlist)
library(utils)
library(dplR)
library(rgl)
library(eegkit)
library(fdatest)
library(spectral)
library(pROC)



####import the dataset from 1 to 19
  data1<- readMat(individual[[ind]][1])
  dis_label <- readMat(individual[[ind]][2])

#diatance label
dis_label <- dis_label[["dist.label"]]

label <- matrix(nrow = 0, ncol = 1)
for (i in 1:48){
  label <- rbind(label,matrix(c(rep(dis_label[i,1],2828)),ncol = 1))
}
colnames(label) <- c('dis_label')

#distance data
d1 <- matrix(nrow = 64, ncol = 0)
for (i in 1:48){
  d1 <- cbind(d1, data1$DATA[1:64,1:2828,i])
}
d2 <- t(d1)



#create time

#time <- as.data.frame(seq(0,by=1/500,length.out=dim(d2)[1]))
time <- as.data.frame(rep(seq(1/500,by=1/500,length.out=2828),48))
colnames(time) = c('time')

#channel
channel <- c('Fp1','Fz','F3','F7','FT9','FC5','FC1','C3',
             'T7','TP9','CP5','CP1','Pz','P3','P7','O1',
             'Oz','O2','P4','P8','TP10','CP6','CP2','Cz',
             'C4','T8','FT10','FC6','FC2','F4','F8','Fp2',
             'AF7','AF3','AFz','F1','F5','FT7','FC3','C1',
             'C5','TP7','CP3','P1','P5','PO7','PO3','POz',
             'PO4','PO8','P6','P2','CPz','CP4','TP8','C6',
             'C2','FC4','FT8','F6','AF8','AF4','F2','Iz')
channel <- t(as.matrix(channel))
colnames(d2) <- channel

#merge
data <- cbind(time, label, d2)


#extrace Long distance
Long <- subset(data, dis_label == 1)
Short<- subset(data, dis_label == 0)
Comb<-rbind(Long,Short)

epoch <- 2828 
aa1 <- dim(Comb)[1]/epoch #48
aa <- seq(1,by=epoch,length.out=aa1) #action 
bb <- seq(3, by=1, length.out = 64) #chanel

##Generate name matrix###
nameLong <- matrix(paste0(rep(paste0(rep(c('Long'), 64),  rep(c("_"), 64), channel), 
                              24), "_",  rep(c(1:24), each = 64)), nrow = 24, byrow = TRUE)
nameShort <- matrix(paste0(rep(paste0(rep(c('Short'), 64),  rep(c("_"), 64), channel), 
                               24), "_",  rep(c(1:24), each = 64)), nrow = 24, byrow = TRUE)
nameComb<-rbind(nameLong,nameShort)

######morlet transformation########
#set.seed(2021)
library(WaveletComp)
w.plot.power <- function(data, row, col){
    i <- aa[row] # row=1,2,...,48, aa[row] denotes the order of row in data
    j <- bb[col] # col=1,2,...,64 order of electrode
    out.wave <- WaveletTransform(as.vector(data[i:(i + 2827), j]),dt=1/500,dj=1/20) 
    sqrt(t(out.wave$Power))
}

#seperate gamma beta alpha theta and delta for 64 channels
# 12 features,   i=1,...48, j=1,...64;
s <- function(obs,ele){
  
  #row denotes the order of row in data48, col denotes order of electrode64
  data <-  w.plot.power(data = Comb, row = obs, col = ele)
  set.seed(1)
  comb1 <- cbind((data[,sample(c(117:160), 3)]),(data[,sample(c(98:116), 3)]),(data[,sample(c(86:97), 3)]), 
                 (data[,sample(c(20:85), 3)])) #(delta,theta,alpha,beta)
  
  cols<- paste0(nameComb[obs, ele], "_", rep(c("delta","theta","alpha","beta"),each=3))
  
  feature <- comb1
  
  colnames(feature) <- cols
  
  #rownames(feature) <- deparse(substitute(data)) #using the name of dataframe as the row name
  
  feature #2828**12
  
}

#generate dataframe for 1 person:48 matrices 2828*768

#column name of each observation
name_rhy <- rep(rep(c("delta","theta","alpha","beta"),each=3),64)
name_ele <- rep(paste0(channel),each=12)
name_obs <- matrix(paste0(name_ele,"_",name_rhy),nrow=1)

#import the data
set.seed(1)
tf.data<-list()
for (obs in 1:48) {
  tf.data[[obs]]<-matrix(0,nrow = 2828, ncol = 64*12)
  for (ele in 1:64) {
    tf.data[[obs]][,c((12*(ele-1)+1):(12*ele))]<-s(obs,ele)
  }
}

#generate list with 768 elements: 2828*48

newdata<-list()
for (i in 1:768) {
  newdata[[i]]<-matrix(0,nrow = 2828, ncol = 48)
  for (j in 1:48) {
    newdata[[i]][,j]<-tf.data[[j]][,i]
  }
}
names(newdata) <- as.vector(name_obs)  


library(WaveletComp)
ma<-data.frame(a=Long$F7[1:2828],b=Short$F7[1:2828])
par(mfrow=c(1,1))
a<-analyze.wavelet(ma,"a",loess.span = 0,dt = 1/500, dj = 1/30,lowerPeriod = 1/30,upperPeriod = 1, make.pval = TRUE, n.sim = 10)
wt.image(a, color.key = "quantile", n.levels = 20,legend.params = list(lab = "wavelet power levels"),main="Long$F7")
b<-analyze.wavelet(ma,"b",loess.span = 0,dt = 1/500, dj = 1/30,lowerPeriod = 1/30,upperPeriod = 1, make.pval = TRUE, n.sim = 10)
wt.image(b, color.key = "quantile", n.levels = 20,legend.params = list(lab = "wavelet power levels"),main="Short$F7")


############################### feature selection #######################

label<-c(rep(1,24),rep(0,24))
pvalue<-NULL
for (i in 1:768) {
  fanova<-fanova.tests(x =newdata[[i]], label, test = "GPF",parallel = FALSE, nslaves = NULL,
                       params=list(paramFP = list(B.FP = 1000, basis = "b-spline",norder = 3)))
  pvalue[i]<-fanova$GPF$pvalueGPF
}

select.fanova<-which(pvalue<=0.1)
feature.select[[ind]]<-select.fanova


#########pMFLR############

Data<-newdata
label<-c(rep(1,24),rep(0,24))
#Time<-c(1:2828)/2828
Time<-seq(1/500,by=1/500,length.out=2828)
select<-select.fanova
name_select<-name_obs[select]
Select<-Data[select]
for (i in 1:length(Select)) {
  Select[[i]]<-t(Select[[i]])
}
index<-c(1:length(Select))

knots <- seq(0,1,length.out=300)
norder <- 4
nbasis <- length(knots) + norder - 2
T <- Time
Trange = c(0,max(Time))
lambda=1e6
n<-48

#Create b-spline basis
bbasisV <- create.bspline.basis(rangeval=Trange,nbasis=nbasis,norder=norder)
curv.Lfd <- int2Lfd(2) #typical for smoothing a cubic spline (order 4),  penalize the squared acceleration
curv.fdParV <- fdPar(bbasisV, curv.Lfd,lambda) #Define A Functional Parameter Object

#choose smoothing parameter
lambdas = 10^seq(-12,5,by=0.5)
mean.gcv = rep(0,length(lambdas))

#Initialize dataframes
SmoothV <- NULL
SmoothVfd <-NULL
A <- NULL
for(j in 1:length(Select)){
  for(ilam in 1:length(lambdas)){
    # Set lambda
    curv.fdParVi <- curv.fdParV
    curv.fdParVi$lambda <- lambdas[ilam]
    # Smooth
    Smoothi <- smooth.basis(T,t(Select[[j]]),curv.fdParVi)
    # Record average gcv
    mean.gcv[ilam] <- mean(Smoothi$gcv)
  }
  #plot(lambdas,mean.gcv,type='b',log='x')
  best = which.min(mean.gcv)
  lambdabest = lambdas[best]
  curv.fdParV$lambda = lambdabest
  bj = smooth.basis(T,t(Select[[j]]),curv.fdParV)
  SmoothV[[j]] <-bj
  #Create A matrices
  c <- SmoothV[[j]]$fd
  SmoothVfd[[j]] <- c
  d <- SmoothVfd[[j]]$coef
  A[[j]] <- t(d)
}

c
#Create PSI matrix
Psi1 <- inprod(bbasisV,bbasisV)
e<- NULL
APsi <- NULL
#multiply Am by Psim
for(i in 1:length(Select)){
  e <- A[[i]]%*%Psi1
  APsi[[i]] <- e
}

#Perform PCA in APsi matrices
##################################################
eigen_value<-NULL
eigen_vector<-NULL
for(i in 1:length(Select)){
  h <- eigen(cor(APsi[[i]]))
  eigen_value[[i]]<-h$value
  eigen_vector[[i]]<-h$vector
}

#Calculate proportion of variance and choose number of PCs
PropVar <- NULL
d<-NULL
e<-NULL
for ( j in 1:length(Select)){
  for( i in 1:length(Select)){
    d[i] <- eigen_value[[j]][i]/sum(eigen_value[[j]])
  }
  e <- cumsum(d)
  PropVar[[j]] <- e
}
NPCV <- NULL
NPC <- NULL
for(j in 1:length(Select)){
  d <- NULL
  for(i in 1:length(Select)){
    if(PropVar[[j]][i] < 0.9){
      d[i] <- 1
    }
    else{ d[i] <- 0}
  }
  NPCV[[j]]<-d
  NPC[j] <- sum(NPCV[[j]])
}

#Select # PCs
min(NPC)
max(NPC)
table(NPC)
numPC <-max(NPC) # mode 

PC <- NULL
for(i in 1:length(Select)){
  d <- eigen_vector[[i]][,1:numPC]
  PC[[i]]<-APsi[[i]]%*%d
}




#traning and validation
for (m in 1:n) {
  
  set.seed(m)
  samp<-c(sample(1:36,0.9*n, replace=TRUE),sample(37:48,0.1*n, replace=TRUE))
  Y<-rep(0,n)
  for (a in 1:n) {
    Y[a]<-label[samp[a]]
  }
  n1<-0.9*n
  train <- samp[1:n1]
  valid <- samp[-(1:n1)]
  #valid <- NULL
  #for(i in 1:n){
  #  if (i %in% train){}
  #  else{
  #    valid <- rbind(valid,i)
  #  }
  #}
  PCtrain <- NULL
  PCvalid <- NULL
  x<- NULL
  
  a <- NULL
  #Training Design
  for(j in 1:length(Select)){
    a <- NULL
    x <- NULL
    for(i in 1:n1){
      a <- PC[[j]][train[i],]
      x<-rbind(x,a)
    }
    PCtrain[[j]]<-x
  }
  
  #Validation Design
  for(j in 1:length(Select)){
    a <- NULL
    x <- NULL
    for(i in 1:(n-n1)){
      a <- PC[[j]][valid[i],]
      x<-rbind(x,a)
    }
    PCvalid[[j]]<-x
  }
  
  #Training Y
  Ytrain <- NULL
  for(i in 1:n1){
    Ytrain[i] <- Y[train[i]]
  }
  
  #Valid Y
  Yvalid <- NULL
  for(i in 1:(n-n1)){
    Yvalid[i] <- Y[valid[i]]
  }
  
  #grp***
  library(robustHD)
  library(grplasso)
  library(calibrate)
  
  #########Validaton Set as Validation
  
  Int <- ones(n1,1)
  lasX<-NULL
  lasX<-cbind(lasX,Int)
  for(i in 1:length(Select)){
    lasX <- cbind(lasX,PCtrain[[i]])
  }
  index <- NULL
  index <- c(index,NA)
  for(i in 1:length(Select)){
    b <-rep(i,numPC)
    index <- c(index,b)
  }
  
  
  #############grplasso##################
  lambdalas <- lambdamax(lasX, y=Ytrain, index=index, penscale = sqrt, model = LogReg()) * 0.5^(0:10)
  lasfit <- grplasso(lasX, y=Ytrain, index= index, lambda = lambdalas, model= LogReg(), penscale=sqrt, control = grpl.control(update.hess = "lambda", trace=0))
  Coefs <- lasfit$coefficients
  dim(Coefs)
  # Check the model
  # as.integer(rowMeans(larsfit$fitted.values)-Ytrain) Correct
  
  CoSUM <- NULL
  for(i in 1:dim(Coefs)[2]){
    CoSUM[i]<-sum(Coefs[2:(length(Select)*numPC + 1),i])
  }
  CoCo <- NULL
  VarNUM <- NULL
  for(j in 1:dim(Coefs)[2]){
    g <- NULL
    for(i in 1:(numPC*length(Select) + 1)){
      if(Coefs[i,j] ==0){
        g[i] <- 0}
      else{g[i] <- 1}
    }
    CoCo[[j]]<-g
  }
  SumCoCo <- NULL
  for(i in 1:dim(Coefs)[2]){
    SumCoCo[i] <- sum(CoCo[[i]][-1])
  }
  
  
  #View number of functional variables selected by grplasso
  #SumCoCo
  
  #Find Lihat and Pihat and classify yhat
  ModelCoef <- Coefs[,2]
  Alpha <- ModelCoef[1]
  Betastar <- ModelCoef[2:(numPC*length(Select) +1)]
  Betastar <- t(Betastar)
  Betastar <- t(Betastar)
  M<-NULL
  for(i in 1:length(Select)){
    b <- PCvalid[[i]]%*%Betastar[(numPC*(i-1)+1):(numPC*i)]
    M[[i]] <- b
  }
  
  APsiB <- 0*ones((n-n1),1)
  for(i in 1:length(Select)){
    APsiB <- APsiB + M[[i]]
  }
  
  
  #Find Lhat
  Lhatvalid <- NULL
  Lhatvalid = Alpha*ones((n-n1),1) + APsiB
  
  
  #Pistar1
  #find Pi(i)??s
  Pstar <- NULL
  YstarV <- NULL
  SEV<-NULL
  for (i in 1:(n-n1)){
    if(Lhatvalid[i]>700){
      Pstar[[i]]<-1
    }
    else{
      Pstar[[i]] <- (exp(Lhatvalid[i]))/(exp(Lhatvalid[i])+1)
    }
    if(Pstar[[i]]>=0.6){
      YstarV[[i]] <- 1
    }
    else{YstarV[[i]]<- 0}
    SEV[[i]]<-(Pstar[[i]]-Yvalid[[i]])^2
  }
  
  a<-Yvalid
  df<-data.frame(a=a,p=as.vector(Pstar))
  odf<-df[order(df$a),]
  if(mean(a)==1){
    auc[m]<-1
  }else{
    roc_obj<-roc(odf$a,odf$p)
    
    auc[m]<-auc(roc_obj)
  }
  
  TP[m]<-(length(which(YstarV==1 & Yvalid== 1)))/(length(which(YstarV==1 & Yvalid== 1))+length(which(YstarV==0 & Yvalid== 1))) 
  FN[m]<-(length(which(YstarV==0 & Yvalid== 1)))/(length(which(YstarV==1 & Yvalid== 1))+length(which(YstarV==0 & Yvalid== 1))) 
  FP[m]<-(length(which(YstarV==1 & Yvalid== 0)))/(length(which(YstarV==1 & Yvalid== 0))+length(which(YstarV==0 & Yvalid== 0))) 
  TN[m]<-(length(which(YstarV==0 & Yvalid== 0)))/(length(which(YstarV==1 & Yvalid== 0))+length(which(YstarV==0 & Yvalid== 0)))
  Acc[m]<-length(which(YstarV==Yvalid))/30
  SE[m]<-mean(SEV)
  
}
