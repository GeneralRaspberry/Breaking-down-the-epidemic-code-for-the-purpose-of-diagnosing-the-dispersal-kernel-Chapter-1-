rm(list=ls())
##RUnning the experiment using the above tau leap function
##---------------------------------------------------------
library('spatstat')
library('ggplot2')
library('dplyr')
library('reshape2')
library('tidyr')
library('ggpubr')
library('patchwork')
library("spatialEco")

## tau-leap Gillespie algorithm function
tauLeapG <- function(beta, # transmission rate
                     theta, # dispersal scale
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputing the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stoping condition 1: incidence lvl
                     t.end=Inf, # stoping condition 2: time after first simulated time step
                     area.host=10, # surface area occupied by one host
                     delta.t=1, # time step
                     ppp, # point pattern as a ppp object, optinally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  
  ## if the point pattern has no marks, generate some randomly that fits q0
  if (is.null(marks(ppp))){
    inf.start <- max(1, round(ppp$n * q0))
    marks(ppp) <- sample(c(rep(FALSE, ppp$n-inf.start), rep(TRUE, inf.start)))
  }
  
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- exp(-pairdist(ppp)^b / theta^b)
    diag(dist.mat) <- NA
  }
  
  ## function that compute infection event probability, based on the dispersal kernel
  k.norm <- beta * area.host * (b/(2*pi*theta^2*gamma(2/b))) # constant part of the exponential power kernel
  infection <- function(infected, dist){
    inf <-  matrix(k.norm * dist[infected,!infected],
                   nrow=sum(infected), byrow=FALSE)
    inf[is.na(inf)] <- 0
    inf
  }
  
  ## starting time
  time <- 0
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
  
  ## computation loop
  while (any(!ppp$marks) & time <= t.end & mean(ppp$marks) < q.end){
    ## infection event probaility
    events <- infection(ppp$marks, dist=dist.mat)
    ## random proisson draws
    new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]
    ## change marks of newly infected
    ppp$marks[new.infected] <- TRUE
    ## increment time
    time <- time + delta.t
    ## if some infection, increment the big dataframe
    if (length(new.infected) > 0){
      df.big <- rbind(df.big, data.frame(time=time, who=new.infected, t(ppp$marks)))
    }
    ## print a dot per new infection
    # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
  }
  
  ## make compact, time only, version of the big dataframe
  times.i <- unique(df.big[,1])
  times.d <- times.i + sigma
  times <- sort(unique(c(times.i, times.d)))
  infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
  detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
  df.small <- data.frame(time=times, infected=infected, detectable=detectable)
  
  ## out put the simplified time series, and the big one
  list(df.small[df.small$time <= max(df.big$time),], df.big) 
} 




## meta parameters
delta.t <- 40 # time step (ALEX-THIS IS BIGGER THAN THE EXPERIMENT BELOW BECAUSE IT IS TAKING SO MUCH LONGER!)
iterations <- 1000 # how many epidemic to simulate
hosts <- 1000 # number of hosts
dim <- 2000 # dimension of the landscape

## epidemic parameters
sigma <- 0.1 #this is the assymptomatic period, doesn't change yet
##but maybe the assymptomatic period is affecting the time till detection?
beta <- 5 ##The data I sent you, which is called data in R is the 1000 realisations of these parameters
theta <- 10
b <- .4


##Concatenating a list of metric values
##-----------------------------------------


## design a function that will be called
sim_par <- function(i=NULL){
  
  set.seed(seed=NULL)
  
  radiusCluster<-100
  lambdaParent<-.02
  lambdaDaughter<-30
  randmod<-0.1
  
  
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
  
  
  
  data <- data.frame(x=landscape$x, y=landscape$y, id=1:hosts)
  
  kk<-Kest(landscape)
  plot(kk)
  kk_iso<-kk$iso
  kk_pois<-kk$theo
  
  kk_div_na<-kk_iso/kk_pois
  kk_div_0<-replace_na(kk_div_na,0)
  kk_mean<-round(mean(kk_div_0),3)
  
  set.seed(seed=NULL)
  marks(landscape) <- sample(c(TRUE, rep(FALSE, hosts-1)))
  
  output <- tauLeapG(beta = beta, theta = theta, b = b,
                     sigma = sigma, delta.t = delta.t,
                     ppp = landscape)
  temp <- output[[2]][,1:2][order(output[[2]][,2]),]
  
  data.frame(time=temp$time, who=temp$who, x=landscape$x[temp$who], y=landscape$y[temp$who],randmod=randmod, metricrel=kk_mean,sim=i) ## what it exports will be concatenated in a list
}

library("parallel")

## create a cluster with the set number of cores, say nmax-1
cl <- makeCluster(mc <- getOption("cl.cores", 3))
## call the library loading function in them
clusterCall(cl, function() library("spatstat"))
clusterCall(cl,function() library("ggplot2"))
clusterCall(cl,function() library("dplyr"))
clusterCall(cl,function() library("tidyr"))
## export all to the nodes, that's dirty, so run this with a clean environement otherwise your memory will be flooded
clusterExport(cl=cl, varlist=ls())
## call the function in a parallel lapply
par_results <- parLapply(1:1000, fun=sim_par, cl=cl) ## test with 10 first, but then replace 10 by 1000
## stop the cluster
stopCluster(cl)
## call cbind on your list of lines to find the matrix you expect
data <- do.call("rbind", par_results)



head(data)
data<-data.frame(data)
times <- sort(unique(data$time))
data_logistic <- function(i=NULL){
data %>% group_by(sim) %>%
  do(data.frame(time=times, infected=sapply(times, function(x) sum(.$time <= x))))
}
## make a logistic df from this data
cl <- makeCluster(mc <- getOption("cl.cores", 3))
clusterCall(cl,function() library("dplyr"))
clusterCall(cl,function() library("tidyr"))
clusterExport(cl=cl, varlist=c("data","times"),envir = environment())
par_data_logistic<-parLapply(1,fun=data_logistic,cl=cl)
stopCluster(cl)
data_log<-data.frame(par_data_logistic)


## prepare a logistic function of r to fit
temp <- filter(data_log, infected < 250)
temp$simdigit<-as.numeric(temp$sim)

r_calculate<-function(i=NULL){

  logis <- function(t, r, K=1, s=0, q0){
  pmin(
    K*q0*exp(r*(t+s)) / (K + q0*(exp(r*(t+s)) - 1)),
    K) # numerical errors can happen for high r and sigma
}


eval <- function(r, df){
  sum((logis(r=r, t=df$time, K=1000, q0=1) - df$infected)^2) ## sum of square errors between predictions and observations
}
temp.a<-subset(temp,simdigit<=200)
temp.b.<-subset(temp,simdigit>=201 & simdigit<=400)
temp.c.<-subset(temp,simdigit>=401 & simdigit<=600)
temp.d.<-subset(temp,simdigit>=601 & simdigit<=800)
temp.e.<-subset(temp,simdigit>=801 & simdigit<=1000)
# sapply(unique(temp$sim), 
#               function(i) optimize(f = eval, interval = c(0, 0.5), df=filter(temp, sim==i))$minimum)
ra <- sapply(unique(temp.a$sim), 
             function(i) optimize(f = eval, interval = c(0, 0.04), df=filter(temp.a, sim==i))$minimum)
rb <- sapply(unique(temp.b.$sim), 
             function(i) optimize(f = eval, interval = c(0, 0.04), df=filter(temp.b., sim==i))$minimum)
rc <- sapply(unique(temp.c.$sim), 
             function(i) optimize(f = eval, interval = c(0, 0.04), df=filter(temp.c., sim==i))$minimum)
rd <- sapply(unique(temp.d.$sim), 
             function(i) optimize(f = eval, interval = c(0, 0.04), df=filter(temp.d., sim==i))$minimum)
re <- sapply(unique(temp.e.$sim), 
             function(i) optimize(f = eval, interval = c(0, 0.04), df=filter(temp.e., sim==i))$minimum)
r<-c(ra,rb,rc,rd,re)
}
#another cluster
cl <- makeCluster(mc <- getOption("cl.cores", 3))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("temp"),envir = environment())
par_r<-parLapply(1,fun=r_calculate,cl=cl)
stopCluster(cl)

logis <- function(t, r, K=1, s=0, q0){
  pmin(
    K*q0*exp(r*(t+s)) / (K + q0*(exp(r*(t+s)) - 1)),
    K) # numerical errors can happen for high r and sigma
}
r_vec<-unlist(par_r)
mean_r<-mean(r_vec)
pred_data <- data.frame(time=times, infected=logis(r=mean_r, t=times, K=1000, q0=1))
ggplot(temp) + geom_point(aes(x=time, y=infected), size=.2) +
  geom_line(data=filter(pred_data, infected<250), aes(x=time, y=infected), colour="red", size=2)+
  ggtitle("Epidemic growth curve for 1000 simulations")
length(unique(r_vec))

#mdata<-melt(data, id.vars = c('x', 'y', 'id'), value.name = 'time', variable.name = 'sim')
#df$detected<-NULL
## if you want to have multiple sampling events then you can add another loop
df<-data
  res <- NULL
 # res<-data.frame(res)
#sur_par_func <- function(i=NULL){
for (sim in unique(df$sim)){ ## looping over simulations
  print(paste("For ",sim))
  n <- 20 ## we sample 20 hosts
  stti<-sample(1:10,1)
  infmax<-max(df$time)
  samp.time <- seq(from = stti, to = infmax, by = 10) ## at sampling time 1 and 4
  ## for ease we make a temporary dataframe of the current simulation
  temp <- df[df$sim==sim,]
  ## we sample n=20 hosts
  ## loop over the sampling times
for (t in samp.time){##here you see that the hosts sampled is after the loop detecting
    ##time that hosts are detected as infected, thus changing every time you update the sampling
    test <- sample(temp$who, n,replace=FALSE)
    print(paste("Id",test))
    
    ## get those hosts infection time
    inf.time <- temp[temp$who %in% test, "time"]
    #print(paste("at time",inf.time))
    ## if an infection time is anterior to the sampling time, we have a detection event
    m <- sum(inf.time <= t) ## we sum to know how many host are seen as infected
    ## we can also measure the true incidence at sampling time
    q <- mean(temp$time<=t)
    ## and increment the result table in the loop
    res <- rbind(res, data.frame(sim, m, q, t=t, n))
    if (m>=1){
     break
    }
  }
}
#}
#cl <- makeCluster(mc <- getOption("cl.cores", 3))
#clusterCall(cl,function() library("dplyr"))
#clusterCall(cl,function() library("tidyr"))
#clusterExport(cl=cl, varlist=c("df","res"),envir = environment())
#par_sur<-parLapply(1,fun=sur_par_func,cl=cl)
#stopCluster(cl)
#sur_data<-data.frame(par_sur)

avg.sur<-res%>% distinct(res$m,res$sim,.keep_all=TRUE)
avg.sur2<-subset(avg.sur, !m==0)
sum.q<-mean(avg.sur2$q)
anq<-(mean_r*10/20)

ggplot(avg.sur2,aes(x=q))+geom_histogram(color="black",fill="lightblue",aes(y=..count../sum(..count..)),binwidth=0.002)+
  geom_vline(aes(xintercept=mean(avg.sur2$q)),color="blue",linetype="dashed")+
  geom_vline(aes(xintercept=(anq)),color="blue")+
  ylab("Probability")+
  xlim(NA,1)+
  ylim(NA,0.2)


absdif<-abs(anq-sum.q)
reldif<-absdif/sum.q #change this to absdif/anq!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

metricrelmean<-mean(df$metricrel)
randmod<-mean(df$randmod)
tableofvalues<-data.frame(c(randmod=randmod,metricrel=metricrelmean,prediction=anq,simulation=sum.q,absdif=absdif,reldif=reldif))
alarm()
tableofvaluesRf0.1theta10beta5<-tableofvalues
save(data=tableofvaluesRf0.1theta10beta5,file="tableofvaluesRf0.1theta10beta5.Rda")
