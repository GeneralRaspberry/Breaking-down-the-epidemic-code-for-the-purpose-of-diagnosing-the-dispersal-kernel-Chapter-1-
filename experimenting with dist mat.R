packagedelivery<-function(fry,leela){
  if(fry == TRUE){
    
    
    require(leela,character.only = TRUE)
  } else{
    x<-grepl(leela,search())
    n<-0
    for (j in x){
      n<-n+1
      if (j == TRUE){
        detach(pos=n, unload=TRUE, character.only = TRUE)
      }
    }
  }
}



packagedelivery(TRUE,"spatstat")
packagedelivery(TRUE, "dplyr")
packagedelivery(TRUE, "ggplot2")
packagedelivery(TRUE, "ggpubr")
packagedelivery(TRUE,"RColorBrewer")
packagedelivery(TRUE, "rdist")

#generate a marks object

radiusCluster<-100
lambdaParent<-.02
lambdaDaughter<-30
randmod<-1
hosts<-1000
dim<-2000

numbparents<-rpois(1,lambdaParent*dim)

xxParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)
yyParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)

numbdaughter<-rpois(numbparents,(lambdaDaughter))
sumdaughter<-sum(numbdaughter)


#theta<-2*pi*runif(sumdaughter)
thetaLandscape<-2*pi*runif(sumdaughter)

rho<-radiusCluster*sqrt(runif(sumdaughter))

# xx0=rho*cos(theta)
# yy0=rho*sin(theta)
xx0=rho*cos(thetaLandscape)
yy0=rho*sin(thetaLandscape)


xx<-rep(xxParent,numbdaughter)
yy<-rep(yyParent,numbdaughter)

xx<-xx+xx0

yy<-yy+yy0
cds<-data.frame(xx,yy)
is_outlier<-function(x){
  x > dim| x < 0
}
cds<-cds[!(is_outlier(cds$xx)|is_outlier(cds$yy)),]
while (nrow(cds)<hosts){
  dif<-hosts-nrow(cds)
  extraparentxx<-sample(xxParent,dif,replace = TRUE)
  extraparentyy<-sample(yyParent,dif,replace = TRUE)
  extrathetaLandscape<-2*pi*runif(dif)
  extrarho<-radiusCluster*sqrt(runif(dif))
  newextracoodsxx<-extrarho*cos(extrathetaLandscape)
  newextracoodsyy<-extrarho*sin(extrathetaLandscape)
  extraxx<-extraparentxx+newextracoodsxx
  extrayy<-extraparentyy+newextracoodsyy
  cdsextra<-data.frame(xx=extraxx,yy=extrayy)
  cds<-rbind(cds,cdsextra)
}
#cds<-rbind(cds,cdsextra)

sampleselect<-sample(1:nrow(cds),hosts,replace=F)
cds<-cds%>%slice(sampleselect)

randfunction<-function(x){
  x<-runif(length(x),0,dim)
}
randselect<-sample(1:nrow(cds),floor(hosts*randmod),replace=F)
cds[randselect,]<-apply(cds[randselect,],1,randfunction)

landscape<-ppp(x=cds$xx,y=cds$yy,window=owin(xrange=c(0,dim),yrange=c(0,dim)))

marks(landscape) <- sample(c(TRUE, rep(FALSE, hosts-1)))

###################################changing the status of a mark#############################################
##################################using pairdist##############################################################

#dist.mat <- exp(-pairdist(ppp)^b / theta^b)
#diag(dist.mat) <- NA

landscape$marks[16]<-TRUE
landscape$marks[502]<-FALSE
dist.mat<-pairdist(landscape)

dl<-data.frame(landscape)

which(landscape$marks)
dist.mat.refined<-dist.mat[landscape$marks,]
dl<-cbind(dl,dist.mat.refined)


#################################plotting with ggplot########################################################
################using colour brewer#########################################################################

myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
ggplot(dl)+geom_point(aes(x,y,colour=dist.mat.refined))+coord_equal()+theme_minimal()+
  scale_color_gradientn(colors = myPalette(1000))

#################################creating a dispersal kernel################################################
dispersalgraphgenerator<-function (x,theta,beta,normtype){
  dist.mat<-pairdist(landscape)
  dl<-data.frame(landscape)
  dist.mat.refined<-dist.mat[landscape$marks,]
  dl<-cbind(dl,dist.mat.refined)
  
dist.mat.kernel<-exp(-dist.mat/theta)*beta
dist.mat.kernel.refined<-dist.mat.kernel[landscape$marks,]
dl<-cbind(dl,dist.mat.kernel.refined)
if (normtype==2){
normkernel<-dist.mat.kernel.refined*normfactor
dl<-cbind(dl,normkernel,normfactor)
} else {
  denominator <- 0
  for(i in 1:length(landscape$marks))
  {
    for(j in 1:length(landscape$marks))
    {
      if(i != j)
      {
        denominator <- denominator + exp(-alpha * pairdist()[i,j])
      }
    }
  }
  normFactor <- nHosts / denominator
}
}
}

theta<-500
beta<-50
alphasqr<-1/(theta*theta)
normfactor<-alphasqr*1/(2*pi)
dl<-dispersalgraphgenerator(landscape,theta,beta)


plot_data_column<-function(data,column){

ggplot(data)+geom_point(aes(x,y,colour=column))+coord_equal()+theme_minimal()+
  scale_color_gradientn(colors=myPalette(1000))
  
}

myplots<-lapply(dl[,4:7], plot_data_column, data=dl)

plot.theta500.beta50<-ggarrange(myplots[[1]],myplots[[2]],myplots[[3]],myplots[[4]],nrow = 2,ncol = 2)

ggsave("theta500beta50.png",plot.theta500.beta50,width=50,height = 50, units= "cm")



 ################################recognising that the dist.mat function is simply an index call#############


landscape$marks[16]<-TRUE
landscape$marks[725]<-TRUE
dist.mat<-pairdist(landscape)
dist.mat.refined<-data.frame(dist.mat[landscape$marks,!landscape$marks])

##############################checking pdist functionality################################################
xtest<-c(3,4,56,6,4,46,4,4,6,4,5,64,4,5)
ytest<-c(3,4,56,6,4,46,4,4,6,4,5,64,4,5)
dftest<-data.frame(xtest,ytest)
distcheck<-pdist(dftest)
