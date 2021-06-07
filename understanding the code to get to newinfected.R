##########################you can't have new.infected without events############################################




events <- infection(ppp$marks, dist=dist.mat)


k.norm <- beta * area.host * (b/(2*pi*theta^2*gamma(2/b)))


infection <- function(infected, dist){
  inf <-  matrix(k.norm * dist[infected,!infected],
                 nrow=sum(infected), byrow=FALSE)
  inf[is.na(inf)] <- 0
  inf
}

new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]

infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))

##########################Loading necessary packages#####################################################
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
#####################creating an object to look at df.big##################################################

df.big <- data.frame(time=0, who=which(landscape$marks), t(landscape$marks))

#################################infected##################################################################

infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))

#object times not found


