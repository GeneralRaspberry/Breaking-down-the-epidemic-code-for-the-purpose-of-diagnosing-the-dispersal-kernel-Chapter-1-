library('spatstat')


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

#####################to bake an apple pie you must first create the universe###############################
#####################creating an object to look at new.infected############################################
###new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]###

df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
s=sum(!landscape$marks) #output is number of marks not infected

which(!landscape$marks)[s]#output is the host ID for the last non-infected host,
                          #because sum(!landscape$marks) gives a length of noninfected

s=length(!landscape$marks)#length landscape$marks reports ALL mark status, effectively ignoring the !

#!landscape$marks actually reports the status of the mark, for example if there is a mark then status reports false!

s=length(which(landscape$marks))
