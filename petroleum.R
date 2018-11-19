library(readxl)
mydata = read_xlsx("PetroleumConsumedHouseHold.xlsx", col_names = TRUE)
petroleum = ts(mydata[,3])
plot.ts(petroleum, ylab = 'Petroleum Consumed', main = 'Plot of Petroleum Consumed')
petroleumT1 = log(petroleum)
petroleumT2 = sqrt(petroleum)
petroleumT3 = petroleum^{-1/2}
par(mfrow = c(2,2))
plot.ts(petroleum, ylab = 'petroleum', main = 'Plot of petroleum')
plot.ts(petroleumT1, ylab = 'ln(petroleum)', main = 'Plot of ln(petroleum)')
plot.ts(petroleumT2, ylab = 'sqrt(petroleum)', main = 'Plot of sqrt(petroleum)')
plot.ts(petroleumT3, ylab = '1/sqrt(petroleum)', main = 'Plot of 1/sqrt(petroleum)')
par(mfrow = c(2,2))
hist(petroleum)
hist(petroleumT1)
hist(petroleumT2)
hist(petroleumT3)
time = 1:395
Model0 = lm(petroleum~time)
Model1 = lm(petroleumT1~time)
Model2 = lm(petroleumT2~time)
Model3 = lm(petroleumT3~time)
par(mfrow = c(2,2))
plot(Model0$residuals)
plot(Model1$residuals)
plot(Model2$residuals)
plot(Model3$residuals)

t = 1:395
model02 = lm(petroleumT1~poly(t,2))
summary(model02)
par(mfrow = c(2,2))
plot(model02$fitted.values)
plot(model02$residuals)
acf(model02$residuals)
pacf(model02$residuals)
AIC(model02)

model03 = lm(petroleumT1 ~poly(t,3))
summary(model03)
par(mfrow = c(2,2))
plot(model03$fitted.values)
plot(model03$residuals)
acf(model03$residuals)
pacf(model03$residuals)
AIC(model03)

model04 = lm(petroleumT1~poly(t,4))
summary(model04)
par(mfrow = c(2,2))
plot(model04$fitted.values)
plot(model04$residuals)
acf(model04$residuals)
pacf(model04$residuals)
AIC(model04)

model05 = lm(petroleumT1~poly(t,5))
summary(model05)
par(mfrow = c(2,2))
plot(model05$fitted.values)
plot(model05$residuals)
acf(model05$residuals)
pacf(model05$residuals)
AIC(model05)
model06 = lm(petroleumT1~poly(t,6))
summary(model06)
par(mfrow = c(2,2))
plot(model06$fitted.values)
plot(model06$residuals)
acf(model06$residuals)
pacf(model06$residuals)
AIC(model06)
model07 = lm(petroleumT1~poly(t,7))
summary(model07)
par(mfrow = c(2,2))
plot(model07$fitted.values)
plot(model07$residuals)
acf(model07$residuals)
pacf(model07$residuals)
AIC(model07)
model08 = lm(petroleumT1~poly(t,8))
summary(model08)
par(mfrow = c(2,2))
plot(model08$fitted.values)
plot(model08$residuals)
acf(model08$residuals)
pacf(model08$residuals)
AIC(model08)

trndseas=function(y,seas=1,lam=1,degtrnd=0){
  
  m=length(lam)
  n=length(y)
  
  # Part of design matrix for estimating trend
  if(degtrnd>0) {
    tm=seq(1/n,1,by=1/n) #use normalized time
    x1=poly(tm,degree=degtrnd,raw=TRUE) #generate the matrix of x
    x1=cbind(rep(1,n),x1)
  } else {
    x1=as.matrix(rep(1,n),ncol=1) #x matrix with only intercept
  }
  
  # Part of design matrix for estimating seasonality
  x2=NULL
  if(seas>1){
    sn=rep(1:seas,length.out=n)
    x2=factor(sn,levels=unique(sn),ordered=TRUE)
    x2=model.matrix(~x2-1) #matrix without the intercept
    m2=ncol(x2)
    m21=m2-1
    x2=x2[,1:m21]-matrix(rep(x2[,m2],m21),ncol=m21,nrow=nrow(x2),byrow=FALSE) #include the negative ones
  }
  
  x=cbind(x1,x2)  # design matrix
  
  xx=t(x)%*%x
  rsq=rep(1,m) #empty vector of length m for computation of r squared of each transformation fit
  m1=ncol(x1)     #degtrnd+1
  m11=m1+1
  mx=ncol(x)      # degtrnd+1+seas-1
  
  for(i in 1:m) { #m is the length of lambda, do Boxcox transformation
    if (lam[i]==0) {
      yt=log(y)
    } else {
      yt=y^lam[i]
    }
    xy=t(x)%*%yt
    coef=solve(xx,xy)
    fit=x%*%coef #this is the trend and seasonal fit
    res=yt-fit
    ssto=(n-1)*var(yt)
    sse=t(res)%*%res
    rsq[i]=1-sse/ssto
  }
  
  ii=which.max(rsq) 
  lamopt=lam[ii]   #choose lambda optimal according to r squared
  if (lamopt==0) {
    yt=log(y)
  } else {
    yt=y^lamopt
  } #optimal transformation done ; yt
  xy=t(x)%*%yt
  coef=solve(xx,xy)
  fit=x%*%coef
  trnd=x1%*%coef[1:m1] #extract the trend part
  season=NULL
  if(seas>1){
    season=c(coef[m11:mx],-sum(coef[m11:mx]))
    #season = x2%*%coef[m11:mx]
  }
  res=yt-fit
  
  result=list(coef=coef,fit=fit,trend=trnd,res=res,season=season,rsq=rsq,lamopt=lamopt)
  return(result)
}
lam = seq(-1,1,by=0.05)
ff = trndseas(petroleumT1,seas = 12,lam = 1,degtrnd = 5)
rsq = ff$rsq
rsq
attributes(ff)
m.fit = ff$trend
ff$season

n = length(petroleumT1)
s.fit = rep(ff$season,length.out=n)
smooth.fit = ff$fit
par(mfrow=c(2,2))
plot.ts(petroleumT1)
plot.ts(m.fit, main='Estimated Trend')
plot.ts(s.fit,main='Estimated Seasonal Component')
plot.ts(petroleumT1,main='Estimated Smooth Part')
points(smooth.fit,type='l',col='red')
par(mfrow=c(1,1))
x = petroleumT1-m.fit-s.fit
acf(x)
pacf(x)
hist(x)
qqnorm(x)
qqline(x)

acf(petroleumT1)
pacf(petroleumT1)
par(mfrow=c(2,2))
pgrm = spec.pgram(petroleumT1, plot = FALSE, log="no")
plot(pgrm$freq,pgrm$spec,type='l',xlab='Frequency',ylab='')
abline(v=1/12, col="red",lty="dotted")
pgrm = spec.pgram(petroleumT1,spans=3, log = "no")
pgrm = spec.pgram(petroleumT1,spans=7, log = "no")
pgrm = spec.pgram(petroleumT1,spans=11, log = "no")

library(forecast)
fit = auto.arima(x, D=0, max.P =8, max.Q = 8, start.P =0, start.Q  = 8)
fit
plot.ts(fit$residuals)
acf(fit$residuals)
pacf(fit$residuals)
Box.test(fit$residuals,lag=10,'Ljung-Box')

par(mfrow=c(1,1))
fit$coef
library(astsa)
coef.ar = fit$coef[1:5]
coef.ma = fit$coef[6:8]
sigma2 = fit$sigma2
fitspec = arma.spec(ar=coef.ar, ma=coef.ma, var.noise = sigma2, log='no')
par(new=TRUE)
plot(fitspec$freq,fitspec$spec,type='l',xlab='Frequency',ylab='')
plot(pgrm$freq,pgrm$spec,type='l', add=TRUE)

acf(fit$residuals)
pacf(fit$residuals)
par(mfrow=c(2,2))
spec.pgram(fit$residuals,plot=TRUE, log='no')
spec.pgram(fit$residuals,spans=3, log = "no")
spec.pgram(fit$residuals,spans=7, log = "no")
spec.pgram(fit$residuals,spans=11, log = "no")

summary(fit)


t2 = 1:384
Yt = petroleumT1[1:384]
modelT = lm(Yt~poly(t2,5))
par(mfrow = c(2,2))
plot(modelT$fitted.values)
plot(modelT$residuals)
acf(modelT$residuals)
pacf(modelT$residuals)
AIC(modelT)
lam = seq(-1,1,by=0.05)
ff1 = trndseas(Yt,seas = 12,lam = 1,degtrnd = 5)
rsq1 = ff1$rsq
rsq1
attributes(ff1)
m.fit1 = ff1$trend
season1 = ff1$season
season1
n1 = length(Yt)
s.fit1 = rep(season1,length.out=n1)
smooth.fit1 = ff1$fit
par(mfrow=c(2,2))
plot.ts(Yt)
plot.ts(m.fit1, main='Estimated Trend on data except for year 2016')
plot.ts(s.fit1,main='Estimated Seasonal Component on data except for year 2016')
plot.ts(Yt,main='Estimated Smooth Part on data except year 2016')
points(smooth.fit1,type='l',col='red')
par(mfrow=c(1,1))
plot.ts(ff1$res,type = 'l', main = "Estimated rough on data except year 2016")
par(mfrow=c(2,2))
x1 = Yt-m.fit1-s.fit1
acf(x1)
pacf(x1)
hist(x1)
qqnorm(x1)
qqline(x1)
fitT = arima(x1,order=c(5,0,3))
par(mfrow=c(2,2))
plot.ts(fitT$residuals)
acf(fitT$residuals)
pacf(fitT$residuals)
Box.test(fitT$residuals,lag=10,'Ljung-Box')

library(Hmisc)
trend2016 = approxExtrap(m.fit[1:395],m.fit1[1:384],xout=m.fit[385:395], method="linear")[1]
trend2016
petroleumT = petroleumT1[385:395]
ff2 = trndseas(petroleumT1[385:395],seas = 12,lam = 1,degtrnd = 5)
ff2$season
h = 7
deg = 2
coef = ff1$coef[1:(deg+1)]
time1 = (n1+(1:h))/n1
predmat = matrix(rep(time1,deg)^rep(1:deg,each=h),nrow=h,byrow=FALSE)
predmat = cbind(rep(1,h),predmat)
predmat
m.fc = predmat %*% coef
s.fc = rep(ff1$season,length.out=n1+h)
s.fc = s.fc[-(1:n1)]
s.fc
fcast = predict(fitTAR6,n.ahead=h)
x.fc = fcast$pred
x.fc
y.fc = m.fc + s.fc + x.fc
y.fc
plot.ts(Yt,xlim=c(0,n1+h))
points(x=n1+1:h, y=y.fc, col='purple',type='b',pch=19)