library(readxl)
splinetru <- read_excel("D:/Download/data_uas.xlsx")
splinetru

library(ggplot2)

#boxplot
stacked_df <- stack(splinetru[, -1])

bp <- boxplot(stacked_df$values ~ stacked_df$ind,
        col = rainbow(ncol(splinetru)), xaxt = "n")
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.3, bp$names, srt = 45, xpd = TRUE,adj = 0.965)

sp1 <- ggplot(data = splinetru, aes(x = Paud , y = IKPS)) +
  geom_point(color = 'red') +
  stat_smooth()
sp2 <- ggplot(data = splinetru, aes(x = Imunisasi , y = IKPS)) +
  geom_point(color = 'red') +
  stat_smooth()
sp3 <- ggplot(data = splinetru, aes(x = Sanitasi , y = IKPS)) +
  geom_point(color = 'red') +
  stat_smooth()
sp4 <- ggplot(data = splinetru, aes(x = Jaminan_kes , y = IKPS)) +
  geom_point(color = 'red') +
  stat_smooth()
sp5 <- ggplot(data = splinetru, aes(x = Konsumsi, y = IKPS)) +
  geom_point(color = 'red') +
  stat_smooth()

library(gridExtra)
grid.arrange(sp1, sp2, sp3, sp4,sp5, ncol=2)


#============TITIK KNOT OPTIMAL=================
truncated<-function(x,orde,k)
  # fungsi membuat matriks x truncated untuk x=1 variabel 
  # k dibuat menjadi vektor yang elemennya adalah knot knot yang akan digunakan
{
  library(pracma)
  orde=1
  x<-as.vector(x)
  k<-as.vector(k)
  d<-length(k)
  trun<-matrix(0,nrow=length(x),ncol=orde+d+1)
  xtrun<-matrix(0,nrow=length(x),ncol=orde+d+1)
  for (j in 1:(orde+1)){
    xtrun[,j]<-x^(j-1)}
  for (r in 1:d){
    for (i in 1:length(x)){
      trun[i,r+orde+1]<-ifelse(x[i]>=k[r],(x[i]-k[r])^(orde),0)
    }
  }
  for (j in (orde+2):(orde+1+d)){
    xtrun[,j]=trun[,j]
  }
  H=xtrun%*%pinv(t(xtrun)%*%xtrun)%*%t(xtrun)
  hasil=list(xtrun=xtrun,H=H)
  return(hasil)
}

mtruncated<-function(x,knot){
  # fungsi membuat matriks x truncated untuk x=p-variabel
  #knot=list dari k untuk masing masing variabel
  library(pracma)
  x=as.matrix(x)
  p=ncol(x)
  xp=vector("list",p)
  for (r in 1:p){
    xp[[r]]=truncated(x[,r],orde,knot[[r]])$xtrun[,-1]
  }
  xt=do.call("cbind",xp)
  xtruncated=cbind(rep(1,nrow(x)),xt)
  H=xtruncated%*%pinv(t(xtruncated)%*%xtruncated)%*%t(xtruncated)
  hasil=list(xtruncated=xtruncated,H=H)
  return(hasil)    
}

gcv1knot<-function(x,y,b){
  #fungsi mencari knot optimum untuk 1 knot untuk masing masing variabel
  #b adalah banyak titik yang dicoba sebagai knot
  x=as.matrix(x)
  y=as.matrix(y)
  gcv<-as.vector(0)
  knot=matrix(0,nrow=b,ncol=ncol(x))
  knott=matrix(0,nrow=b-2,ncol=ncol(x))
  knots=vector("list",b)
  for (j in 1:ncol(x)){
    knot[,j]=seq(min(x[,j]),max(x[,j]),length.out = b)
    knott[,j]=knot[c(-1,-b),j]
  }
  kn=vector("list",ncol(x))
  for (i in 1:nrow(knott)){
    knots[[i]]=knott[i,]
    kn=as.list(knots[[i]])
    gcv[i]=((1/nrow(x))*t(y)%*%(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kn)$H)%*%y)/(((1/nrow(x))*sum(diag(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kn)$H)))^2)
  }
  optimum=min(gcv)
  indeks.knotoptimum=which.min(gcv)
  knot.optimum=knott[indeks.knotoptimum,]
  hasil<-list(knot=knott,gcv=gcv,gcv.minimum=optimum,knot.optimum=knot.optimum)
  return(hasil)
}

gcv2knot<-function(x,y,b){
  #b adalah banyak titik yang dicoba sebagai knot
  x=as.matrix(x)
  y=as.matrix(y)
  gcv<-as.vector(0)
  knot=matrix(0,nrow=b,ncol=ncol(x))
  knott=matrix(0,nrow=b-2,ncol=ncol(x))
  kn=vector("list",ncol(x))
  a=vector("list",ncol(x))
  for (j in 1:ncol(x)){
    knot[,j]=seq(min(x[,j]),max(x[,j]),length.out = b)
    knott[,j]=knot[c(-1,-b),j]
    a[[j]]=t(combn(knott[,j],2))
  }
  for (i in 1:nrow(a[[1]])){
    kn=lapply(a,"[",i,1:2)
    gcv[i]=((1/nrow(x))*t(y)%*%(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kn)$H)%*%y)/(((1/nrow(x))*sum(diag(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kn)$H)))^2)
  }
  optimum=min(gcv)
  indeks.knotoptimum=which.min(gcv)
  knot.optimum=lapply(a,"[",indeks.knotoptimum,1:2)
  hasil<-list(gcv=gcv,gcv.minimum=optimum,knot.optimum=knot.optimum)
  return(hasil)
}

gcv3knot<-function(x,y,b){
  #b adalah banyak titik yang dicoba sebagai knot
  x=as.matrix(x)
  y=as.matrix(y)
  gcv<-as.vector(0)
  knot=matrix(0,nrow=b,ncol=ncol(x))
  knott=matrix(0,nrow=b-2,ncol=ncol(x))
  kn=vector("list",ncol(x))
  a=vector("list",ncol(x))
  for (j in 1:ncol(x)){
    knot[,j]=seq(min(x[,j]),max(x[,j]),length.out = b)
    knott[,j]=knot[c(-1,-b),j]
    a[[j]]=t(combn(knott[,j],3))
  }
  for (i in 1:nrow(a[[1]])){
    kn=lapply(a,"[",i,1:3)
    gcv[i]=((1/nrow(x))*t(y)%*%(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kn)$H)%*%y)/(((1/nrow(x))*sum(diag(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kn)$H)))^2)
  }
  optimum=min(gcv)
  indeks.knotoptimum=which.min(gcv)
  knot.optimum=lapply(a,"[",indeks.knotoptimum,1:3)
  hasil<-list(gcv=gcv,gcv.minimum=optimum,knot.optimum=knot.optimum)
  return(hasil)
}

kombinasi.knot=function(x,y,b){
  x=as.matrix(x)
  y=as.matrix(y)
  gcv<-as.vector(0)
  k1=gcv1knot(x,y,b)
  knot1=k1$knot.optimum
  k2=gcv2knot(x,y,b)
  knot2=k2$knot.optimum
  k3=gcv3knot(x,y,b)
  knot3=k3$knot.optimum
  knott=vector("list",ncol(x))
  kno=vector("list",ncol(x))
  for (i in 1:ncol(x)){
    knott[[i]]=list(knot1[[i]],knot2[[i]],knot3[[i]])
  }
  kom=expand.grid(knott)
  for (i in 1:nrow(kom)){
    kno=lapply(kom,"[[",i)
    gcv[i]=((1/nrow(x))*t(y)%*%(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kno)$H)%*%y)/(((1/nrow(x))*sum(diag(diag(1,nrow=nrow(x),ncol=nrow(x))-mtruncated(x,kno)$H)))^2)
  }
  optimum=min(gcv)
  indeks.knotoptimum=which.min(gcv)
  knot.optimum=lapply(kom,"[[",indeks.knotoptimum)
  hasil<-list(gcv=gcv,gcv.1knot=k1$gcv.minimum,gcv.2knot=k2$gcv.minimum,gcv.3knot=k3$gcv.minimum,satuknot.optimum=k1$knot.optimum,duaknot.optimum=k2$knot.optimum,tigaknot.optimum=k3$knot.optimum,gcv.minimum=optimum,knot.optimum=knot.optimum)
  return(hasil)
}

hitung<-function(x,y,knot,taraf.alpha){
  # fungsi penaksiran beta dan pengujian hipotesis pada suatu titik knot yang ditentukan dan suatu taraf alpha
  # knot adalah suatu list dari vektor knot maksimal 3 knot untuk masing masing variabel
  library(pracma)
  library(car)
  x=as.matrix(x)
  y=as.matrix(y)
  I=diag(1,nrow=nrow(x),ncol=nrow(x))
  J=matrix(1,nrow=nrow(x),ncol=nrow(x))
  matriksX=mtruncated(x,knot)$xtruncated
  Ak=mtruncated(x,knot)$H
  beta=pinv(t(matriksX)%*%matriksX)%*%t(matriksX)%*%y
  SSE=t(y)%*%(I-Ak)%*%y
  SST=t(y)%*%(I-(1/nrow(x))*J)%*%y
  SSR=SST-SSE
  R.square=as.numeric((SSR/SST)*100)
  dbr=ncol(matriksX)-1
  dbe=nrow(x)-ncol(matriksX)
  dbt=nrow(x)-1
  MSE=SSE/dbe
  MSR=SSR/dbr
  Se=sqrt(MSE*diag(pinv(t(matriksX)%*%matriksX)))
  Fhitung=MSR/MSE
  p.valueF=pf(Fhitung,dbr,dbe,lower.tail=FALSE )
  keputusan=NULL
  keputusan.1=NULL
  if (p.valueF >= taraf.alpha)
  {keputusan.1=c("gagal tolak Hnol")}
  else {keputusan.1=c("tolak Hnol")}
  for (i in 1:length(Se)){
    thitung=beta/Se
    p.value=2*pt(abs(thitung),df=nrow(x)-1,lower=FALSE) 
    if (p.value[i]>=taraf.alpha) {
      keputusan[i]=c("gagal tolak Hnol")
    } else{
      keputusan[i]=c("tolak Hnol")
    }
  }
  ytopi=Ak%*%y
  b<-data.frame(Sumber=c("regresi","error","total"),Sum.of.Square=c(round(SSR,4),round(SSE,4),round(SST,4)),db=c(dbr,dbe,dbt),mean.of.square=c(round(MSR,4),round(MSE,4),"-"),Fhitung=c(round(Fhitung,4),"-","-"),p.value=c(p.valueF,"-","-"),keputusan=c(keputusan.1,"-","-"))
  inferensi=data.frame(beta=beta,standar.error=Se,t.hitung=thitung,p.value=p.value,keputusan=keputusan)
  hasil<-list(ytopi=ytopi,beta=beta,tabel.anova=b,inferensi=inferensi,R.square=R.square)
  return(hasil)
}

spline.truncated.linier.multivariabel=function(x,y,b,taraf.alpha){
  #fungsi utama spline linier multivariabel 
  #fungsi ini memanggil fungsi fungsi yang lain dalam membentuk matriks truncated, penaksiran parameter, dan mencari gcv minimum 
  #b adalah banyaknya titik yang dicoba sebagai knot untuk masing masing variabel dalam trial error pencarian knot optimum
  gcv=kombinasi.knot(x,y,b)
  gcvsatuknot=gcv$gcv.1knot
  satuknot=gcv$satuknot.optimum
  gcvduaknot=gcv$gcv.2knot
  duaknot=gcv$duaknot.optimum
  gcvtigaknot=gcv$gcv.3knot
  tigaknot=gcv$tigaknot.optimum
  knot.optimum=gcv$knot.optimum
  gcv.minimum=gcv$gcv.minimum
  beta=hitung(x,y,knot.optimum,taraf.alpha)$beta
  tabel.anova=hitung(x,y,knot.optimum,taraf.alpha)$tabel.anova
  rsquare=hitung(x,y,knot.optimum,taraf.alpha)$R.square
  inferensi=hitung(x,y,knot.optimum,taraf.alpha)$inferensi
  ytopi=hitung(x,y,knot.optimum,taraf.alpha)$ytopi
  residual=y-ytopi
  kenormalan=ks.test(residual,"pnorm",mean=mean(residual),sd=sd(residual))
  if(kenormalan$p.value >= taraf.alpha){
    keputusan.kenormalan=c("residual berdistribusi normal")}
  else { keputusan.kenormalan=c("residual tidak berdistribusi normal")}
  gleyjser=hitung(x,abs(residual),knot.optimum,taraf.alpha)$tabel.anova
  uji.kenormalan.KS=data.frame(D=as.numeric(kenormalan$statistic),p.value=kenormalan$p.value,keputusan=keputusan.kenormalan)
  durbin=durbinWatsonTest(as.vector(residual))
  matrikstruncated=mtruncated(x,knot.optimum)$xtruncated
  p.value.dw=durbinWatsonTest(lm(y~matrikstruncated[,-1]))$p
  if(p.value.dw >= taraf.alpha){
    keputusan.dw=c("gagal tolak Hnol")}
  else { keputusan.dw=c("tolak Hnol")}
  cat("=================================================","\n")
  cat("satu knot optimum untuk masing-masing variabel","\n")
  cat("=================================================","\n")
  cat("","\n")
  cat("GCV minimum = ",gcvsatuknot,"\n")
  cat("","\n")
  for (i in 1:length(satuknot)){
    cat("knot optimum variabel ke-",i,"\n")
    cat(satuknot[[i]],"\n")
    cat("","\n")
  }
  cat("=================================================","\n")
  cat("dua knot optimum untuk masing-masing variabel","\n")
  cat("=================================================","\n")
  cat("","\n")
  cat("GCV minimum = ",gcvduaknot,"\n")
  cat("","\n")
  for (i in 1:length(duaknot)){
    cat("knot optimum variabel ke-",i,"\n")
    cat(duaknot[[i]],"\n")
    cat("","\n")
  }
  cat("=================================================","\n")
  cat("tiga knot optimum untuk masing-masing variabel","\n")
  cat("=================================================","\n")
  cat("","\n")
  cat("GCV minimum = ",gcvtigaknot,"\n")
  cat("","\n")
  for (i in 1:length(tigaknot)){
    cat("knot optimum variabel ke-",i,"\n")
    cat(tigaknot[[i]],"\n")
    cat("","\n")
  }
  cat("==============================================","\n")
  cat("hasil kombinasi knot optimum dan GCV minimum","\n")
  cat("==============================================","\n") 
  cat("GCV minimum = ",gcv.minimum,"\n")
  cat("","\n")
  cat("knot optimum","\n")
  for (i in 1:length(knot.optimum)){
    cat("knot optimum variabel ke-",i,"\n")
    cat(knot.optimum[[i]],"\n")
    cat("","\n")
  }
  cat("","\n")
  cat("=====================================================================================================================","\n")
  cat("tabel anova","\n")
  cat("=====================================================================================================================","\n")
  cat("Sumber         db          SS            MS               Fhit               P-value         keputusan","\n")
  cat("=====================================================================================================================","\n") 
  cat("Regresi      ",tabel.anova[1,3],"     ",round(tabel.anova[1,2],3),"     ",as.character(tabel.anova[1,4]),"       ",as.character(tabel.anova[1,5]),"  ",as.character(tabel.anova[1,6]),"        ",as.character(tabel.anova[1,7]),"\n")
  cat("Error        ",tabel.anova[2,3],"     ",round(tabel.anova[2,2],3),"     ",as.character(tabel.anova[2,4]),"\n")
  cat("Total        ",tabel.anova[3,3],"     ",round(tabel.anova[3,2],3),"\n")
  cat("=====================================================================================================================","\n")
  cat("","\n")
  cat("","\n")
  cat("=====================================================================================================================","\n") 
  cat("Uji kolmogorov-smirnov untuk kenormalan residual ","\n")
  cat("=====================================================================================================================","\n") 
  cat("D              P-value              keputusan","\n")
  cat("=====================================================================================================================","\n") 
  cat(as.numeric(kenormalan$statistic),"         ",kenormalan$p.value,"        ",keputusan.kenormalan,"\n")
  cat("=====================================================================================================================","\n") 
  cat("","\n")
  cat("","\n")
  cat("","\n")
  cat("=====================================================================================================================","\n") 
  cat("uji Gleyjser","\n")
  cat("=====================================================================================================================","\n")  
  cat("Sumber         db          SS             MS             Fhit              P-value               keputusan","\n")
  cat("=====================================================================================================================","\n") 
  cat("Regresi      ",gleyjser[1,3],"      ",round(gleyjser[1,2],3),"      ",as.character(gleyjser[1,4]),"      ",as.character(gleyjser[1,5]),"      ",as.character(gleyjser[1,6]),"      ",as.character(gleyjser[1,7]),"\n")
  cat("Error        ",gleyjser[2,3],"      ",round(gleyjser[2,2],3),"      ",as.character(gleyjser[2,4]),"\n")
  cat("Total        ",gleyjser[3,3],"      ",round(gleyjser[3,2],3),"\n")
  cat("=====================================================================================================================","\n") 
  cat("","\n")
  cat("","\n")
  cat("=====================================================================================================================","\n") 
  cat("uji Durbin Watson ","\n")
  cat("=====================================================================================================================","\n") 
  cat("DW             P-value           keputusan","\n")
  cat("=====================================================================================================================","\n") 
  cat(as.numeric(durbin),"         ",p.value.dw,"        ",keputusan.dw,"\n")
  cat("=====================================================================================================================","\n") 
  hasil=list(beta=beta,ytopi=ytopi)
}


# Extracting Y and X from Data
Y <- splinetru$IKPS
X <- splinetru[, -(1:2)]  # Exclude the Y column

model=spline.truncated.linier.multivariabel(X,Y,b=34,taraf.alpha = 0.05) 

#====Komb Titik Knot====
spline_x1<- bs(splinetru$Paud, knots= c(42.01, 43.70667, 45.40333), degree = 1)
spline_x2 <- bs(splinetru$Imunisasi, knots= 65.28182 , degree = 1)
spline_x3 <- bs(splinetru$Sanitasi, knots= c(71.21, 72.92667, 74.64333), degree = 1)
spline_x4 <- bs(splinetru$Jaminan_kes, knots= c(18.94727, 41.17455), degree = 1)
spline_x5 <- bs(splinetru$Konsumsi, knots= c(6.013939, 21.29121), degree = 1)

model2_komb <- gam(Y ~ spline_x1 + spline_x2 + spline_x3 + spline_x4 + spline_x5, data=splinetru)
summary(model2_komb)
