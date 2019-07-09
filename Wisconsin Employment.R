#Credits
#"Modelling seasonal data with GAMs" by Gavin Simpson
#https://www.fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/

library(tseries)
library(ggplot2)
library(forecast)
library(astsa)
library(mgcv)
library(reshape2)

#Load Data
data <- read.csv(file="/Users/Pepito/Desktop/TS/wisconsin-employment-time-series.csv")
str(data)

#SARIMA modelling

wisc <- ts(data[, "Data"], start=c(1961, 1), end=c(1975, 10), frequency=12)
plot(wisc, ylab="Employees (in thousands)", main="Wisconsin Monthly Employment 1961-1975", col="red")
stl_wisc <- stl(wisc, s.window = "periodic")
plot(stl_wisc, main="STL decomposition of time series")
#Linear upward trend and periodicity of 12 months
#Variance seems to be increasing over time
#These observations suggest non-stationarity
BoxCox.lambda(wisc, method="guerrero") #Suggested Box-Cox transformation parameter is lambda = 0.02
#Log transform the data in an effort to fix heteroskedasticity
log_wisc <- log(wisc)
plot(log_wisc, ylab="Log Employees (in thousands)", main="Log transformed time series", col='red')

#Approach 1
#Remove linear and seasonal trends by differencing at lag 1 and 12
#That is we consider (1-B)(1-B^12)log(Y_t) = W_t
wisc_statio <- diff(diff(log_wisc, lag=12))
plot(wisc_statio, ylab="", main="Log transformed, detrended, deseasonalized time series")
#It appears we have obtained a stationary time series
kpss.test(wisc_statio)
#p-value of KPSS test greater than 0.1
forecast::Acf(wisc_statio, main="W_t ACF", lag=36)
forecast::Pacf(wisc_statio, main="W_t PACF", lag=36)
#ACF suggests q=0, Q=1 while PACF suggests p=0, P=0/1
#Best fit with respect to AIC has p=P=q=0, Q=1, d=D=1
fit_1 <- arima(log_wisc, order=c(0,1,0), seasonal=list(order=c(0,1,1), period=12))
fit_1
#Log likelihood ratio test with more complex model
fit_lr <- arima(log_wisc, order=c(0,1,0), seasonal=list(order=c(0,1,2), period=12))
fit_lr
teststat<- 2*(as.numeric(logLik(fit_lr))-as.numeric(logLik(fit_1)))
pchisq(teststat, df=1, lower.tail=FALSE) #p value for test
#p-value indicates that we can't reject the null
#Diagnostic plots
tsdiag(fit_1) #Get standardized residuals
forecast:Acf(fit_1$resid, main="ACF of residuals")
forecast:Pacf(fit_1$resid, main="PACF of residuals")
cpgram(fit_1$resid, main="Cumulative periodogram of residuals")
#Get normal qq plot of std residuals and ljung box statistic
tmp <- sarima(log_wisc, p = 0, d = 1, q = 0, P=0, D=1, Q=1, S=12)
#Look for dependence in the residuals
forecast:Acf(fit_1$resid^2, lag=24, main='ACF of squared residuals')
forecast:Pacf(fit_1$resid^2, lag=24, main='PACF of squared residuals')
#Forecast with 85 and 95% confidence bands
fore_fit_1 <- forecast(fit_1, h=24)
plot(fore_fit_1, include=48)
#Get back to original scale
fore_fit_1$mean <- exp(fore_fit_1$mean)
fore_fit_1$lower <- exp(fore_fit_1$lower)
fore_fit_1$upper <- exp(fore_fit_1$upper)
fore_fit_1$x <- wisc
plot(fore_fit_1, include = 48)

#Approach 2 (suggested in practical 3)
#Remove deterministic seasonal component from time series
stl_log_wisc <- stl(log_wisc, s.window="periodic", t.window=99, robust = TRUE)
log_wisc_sadj <- seasadj(stl_log_wisc)
#Differentiate once to obtain seemingly stationary time series
plot(log_wisc_sadj, main="Deseasonalized time series")
log_wisc_sadj_diff <- diff(log_wisc_sadj)
plot(log_wisc_sadj_diff, main="First difference of deseasonalized time series")
kpss.test(log_wisc_sadj_diff)
forecast:Acf(log_wisc_sadj_diff, lag=48, main="U_t ACF")
forecast:Pacf(log_wisc_sadj_diff, lag=36, main="U_t PACF")
#Acf and Pacf suggest a SAR(1) component (maybe look for AR/MA ones)
#Best fit has SAR(1) and SMA(1) components
fit_2 <- arima(log_wisc_sadj, order=c(0,1,0), seasonal=list(order=c(1,0,1), period=12))
fit_2
#Log likelihood ratio test with simpler model
fit_lr <- arima(log_wisc_sadj, order=c(0,1,0), seasonal=list(order=c(1,0,0), period=12))
fit_lr
teststat<- 2*(as.numeric(logLik(fit_2))-as.numeric(logLik(fit_lr)))
pchisq(teststat, df=1, lower.tail=FALSE) #p value for test
#p-value indicates that we should reject the null at level alpha=0.05
#Diagnostic plots
tsdiag(fit_2)
forecast:Acf(fit_2$resid)
forecast:Pacf(fit_2$resid)
cpgram(fit_2$resid)
#Get normal qq plot of std residuals and ljung box statistic
tmp <- sarima(log_wisc_sadj, p = 0, d = 1, q = 0, P=1, D=0, Q=1, S=12)
#Look for dependence in the residuals
forecast:Acf(fit_2$resid^2, main='ACF of squared residuals')
forecast:Pacf(fit_2$resid^2, main='PACF of squared residuals')
#Forecast with 85 and 95% confidence bands
fore_fit_2 <- forecast(fit_2, h=24)
fore_fit_2$mean <- fore_fit_2$mean + seasonal(stl_log_wisc)[11:34]
fore_fit_2$lower <- fore_fit_2$lower + seasonal(stl_log_wisc)[11:34]
fore_fit_2$upper <- fore_fit_2$upper + seasonal(stl_log_wisc)[11:34]
fore_fit_2$x <- log_wisc
plot(fore_fit_2, include = 48)
#Get back to original scale
fore_fit_2$mean <- exp(fore_fit_2$mean)
fore_fit_2$lower <- exp(fore_fit_2$lower)
fore_fit_2$upper <- exp(fore_fit_2$upper)
fore_fit_2$x <- wisc
plot(fore_fit_2, include = 48)

#One can perform model selection based on the AIC
#The following loop was used to select models for each approach
#Here we vary the SARMA components while keeping the ARMA ones fixed
AIC.vals <- matrix(0, 2, 2)
for (p in 0:1) {
  for (q in 0:1) {
    f <- arima(log_wisc, order=c(0,1,0), seasonal=list(order=c(p,1,q), period=12))
    AIC.vals[p + 1, q + 1] <- f$aic
  }
}

#Generalized Additive Models

#Create Dataframe from log_wisc time series
Month <-  factor(cycle(wisc), levels = 1:12, labels = month.abb)
wisc_f <- tapply(log_wisc, list(year = floor(time(wisc)), month = Month), c)
wisc_f <- rbind(wisc_f, matrix(0, 2, 12))
rownames(wisc_f)[16] <- "1976"
rownames(wisc_f)[17] <- "1977"
rn <- as.numeric(rownames(wisc_f))
Years <- rn[1]:rn[length(rn)]
wisc_f <- melt(wisc_f,  id=c("month"))
#wisc_f <- wisc_f[c("month","value")]
wisc_f <- wisc_f[c("Var2","value")]
names(wisc_f) <- c("Month","Value")

#The following is based on "Modelling seasonal data with GAMs" by Gavin Simpson
#https://www.fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/

#Add in Year and nMonth for numeric month and a proper Date class
wisc_f <- transform(wisc_f, Year = (Year <- rep(Years, times = 12)),
                    nMonth = rep(1:12, each = length(Years)),
                    Date = as.Date(paste(Year, Month, "15", sep = "-"), format = "%Y-%b-%d"))
#Sort into temporal order
wisc_f <- wisc_f[with(wisc_f, order(Date)), ]
#Add in a Time variable and get rid of last two missing values
wisc_f <- transform(wisc_f, Time = as.numeric(Date) / 1000)
preds <- wisc_f[(nrow(wisc_f)-23):nrow(wisc_f)-2,]
preds <- preds[c("nMonth","Time")]
wisc_f <- wisc_f[1:(nrow(wisc_f)-26),]

head(wisc_f)

#Fitting GAM
fit_3 <- gamm(Value ~ s(nMonth, bs = "cc", k=12) + s(Time), data = wisc_f)
summary(fit_3$gam)

#Smooth terms
layout(matrix(1:2, ncol = 2))
plot(fit_3$gam, scale = 0)
layout(1)

forecast:Acf(resid(fit_3$lme), lag.max = 36, main = "ACF of residuals")
forecast:Pacf(resid(fit_3$lme), lag.max = 24, main = "PACF of residuals")

#A SARIMA model is fitted to the residuals
gam_res <- resid(fit_3$lme)
gam_res <- ts(gam_res, start=c(1961, 1), end=c(1975, 10), frequency=12)
gam_res_statio <- diff(diff(gam_res,lag=12))
plot(gam_res_statio, main="Stationarized time series of residuals", ylab="")
kpss.test(gam_res_statio)
forecast:Acf(gam_res_statio,  main="ACF of stationarized residuals")
forecast:Pacf(gam_res_statio, main="PACF of stationarized residuals")
fit_4 <- arima(gam_res, order=c(0,1,0), seasonal=list(order=c(0,1,1), period=12))
fit_4
#Diagnostics
tsdiag(fit_4)
forecast:Acf(fit_4$resid, lag=36)
forecast:Pacf(fit_4$resid, lag=36)
cpgram(fit_4$resid)
#Get normal qq plot of std residuals and ljung box statistic
tmp <- sarima(gam_res, p = 0, d = 1, q = 0, P=0, D=1, Q=1, S=12)
#Look for dependence in the residuals
forecast:Acf(fit_4$resid^2, main='ACF of squared residuals')
forecast:Pacf(fit_4$resid^2, main='PACF of squared residuals')

#Predict next two years of smooth terms
add <- predict(fit_3$gam, newdata=preds)
add <- as.vector(add)

#Forecast with (most likely wrong) 85 and 95% confidence bands
fore_fit_4 <- forecast(fit_4, h=24)
fore_fit_4$mean <- fore_fit_4$mean + add
fore_fit_4$lower <- fore_fit_4$lower + add
fore_fit_4$upper <- fore_fit_4$upper + add
fore_fit_4$x <- log_wisc
plot(fore_fit_4, include = 48)
#Get back to original scale
fore_fit_4$mean <- exp(fore_fit_4$mean)
fore_fit_4$lower <- exp(fore_fit_4$lower)
fore_fit_4$upper <- exp(fore_fit_4$upper)
fore_fit_4$x <- wisc
plot(fore_fit_4, include = 48)


