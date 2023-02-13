##RScript for Running AP project 

library(dlnm)
library(lme4)
library(MLmetrics)

##athlete data 
data = read.csv("AthleteDataDeid.csv")
data$wind = sqrt(data$ugrd10m^2 + data$vgrd10m^2)
data$spec_hum = data$spfh2m*1000
data$spec_press = data$pressfc/1000
data$celsius

#two pollutant AQI data 
aqi_daily = read.csv("TwoPollutantAQIFinalData.csv")
aqi_daily1 <- as.matrix(aqi_daily[,rep(1:60,each=1)])
colnames(aqi_daily1) <- paste("lag", 1:60, sep="")

#two pollutant summed AQI data 
aqi_sum_daily = read.csv("SummedAQIFinalData.csv")
aqi_sum_daily1 <- as.matrix(aqi_sum_daily[,rep(1:60,each=1)])
colnames(aqi_sum_daily1) <- paste("lag", 1:60, sep="")

#PM2.5 concentration data 
pm_daily = read.csv("PM25FinalData.csv")
pm_daily1 <- as.matrix(pm_daily[,rep(1:60,each=1)])
colnames(pm_daily1) <- paste("lag", 1:60, sep="")

#Ozone concentration data 
oz_daily = read.csv("OzoneFinalData.csv")
oz_daily1 <- as.matrix(oz_daily[,rep(1:60,each=1)])
colnames(oz_daily1) <- paste("lag", 1:60, sep="")


##Run the DLNM code for Two-pollutant AQI
quantile(aqi_daily1, c(0.20, 0.80))
quantile_20 = 36.1
quantile_80 = 55.1 #closest to values used for prediction to 54.3 (80th percentile)
start = 31.1
end = 67.3

lag_df = 5
var_df = 5
lag_selection = 21

aqi_daily1 <- as.matrix(aqi_daily[,rep(0:(lag_selection),each=1)])
colnames(aqi_daily1) <- paste("lag", 0:(lag_selection - 1), sep="")
cbnest1 = crossbasis(aqi_daily1, lag = c(0,(lag_selection - 1)), argvar = list(fun="ns",df= var_df),
                     arglag = list(fun = "ns", df = lag_df))
lm1 = lme4::lmer(data$adj_S ~ cbnest1 + data$PR + data$real_days_year + data$days_since_event +
                   data$farenheit +data$spec_hum + data$wind  + (1|data$Meet) + (1|data$University), data=data)
pred1.aqi <- crosspred(cbnest1, model=lm1, cen=quantile_20,at = start:end,  cumul = TRUE)

exp = quantile_80

###Generate the plots
plot(pred1.aqi, var = exp,  xlab = 'Lag day', ylab = 'Time (seconds)', col="red", ylim = c(-3,10),
     main = "Two-pollutant threshold AQI")
plot(pred1.aqi,"overall",xlab="Two-pollutant threshold AQI", ylab = 'Time (seconds)', 
     col="red",
     main="Overall cumulative association", ylim = c(-10, 35))
abline(v=50,col="black",lwd=2, lty = 2)


##Calculate the cumulative exposure and associated 95% confidence intervals 
check = (exp - start) + 1
a = pred1.aqi$matfit
sum(a[check,])
pred1.aqi$cumfit[check,lag_selection]
pred1.aqi$cumlow[check,lag_selection]
pred1.aqi$cumhigh[check,lag_selection]

##Now we can do the same thing for summed PM2.5 

quantile(aqi_sum_daily1, c(0.20, 0.80))
quantile_20 = 58.5
quantile_80 = 93.5
start =51.5
end = 106.7

lag_df = 5
var_df = 5
lag_selection = 21

aqi_sum_daily1 <- as.matrix(aqi_sum_daily[,rep(0:lag_selection,each=1)])
colnames(aqi_sum_daily1) <- paste("lag", 0:(lag_selection - 1), sep="")
cbnest1 = crossbasis(aqi_sum_daily1, lag = c(0,(lag_selection - 1)), argvar = list(fun="ns",df= var_df), arglag = list(fun = "ns", df = lag_df))
lm1 = lme4::lmer(data$adj_S ~ cbnest1 
                 + data$PR + data$real_days_year + data$days_since_event + data$farenheit + data$wind + data$spec_hum + (1|data$Meet) + (1|data$University), data=data)
pred1.aqi <- crosspred(cbnest1, model=lm1, cen=quantile_20, bylag = 1, at =start:end, cumul = TRUE)

exp = quantile_80

###Generate the plots
plot(pred1.aqi, var = exp, xlab = 'Lag day', ylab = 'Time (seconds)', col="red", ylim = c(-3,10),
     main = "Summed two-pollutant AQI")
plot(pred1.aqi,"overall",xlab="Summed two-pollutant threshold AQI", ylab = 'Time (seconds)', 
     col="red",
     main="Overall cumulative association", ylim = c(-10, 35))

##Calculate the cumulative exposure and associated 95% confidence intervals 
check = (exp - start) + 1
pred1.aqi$cumfit[check,lag_selection]
pred1.aqi$cumlow[check, lag_selection]
pred1.aqi$cumhigh[check,lag_selection]


##Code for PM2.5 
quantile(pm_daily1, c(0.20, 0.80))
quantile_20 = 5
quantile_80 = 10

lag_selection = 21
lag_df = 5
var_df = 5

start = 4
end = 13

pm_daily1 <- as.matrix(pm_daily[,rep(0:lag_selection,each=1)])
colnames(pm_daily1) <- paste("lag", 0:(lag_selection - 1), sep="")
cbnest1 = crossbasis(pm_daily1, lag = c(0,(lag_selection - 1)), argvar = list(fun="ns",df=var_df), arglag = list(fun = "ns", df = lag_df))
lm1 = lme4::lmer(data$adj_S ~ cbnest1 +
                   data$PR + data$real_days_year + data$days_since_event +
                   data$farenheit + data$wind + data$spec_hum + (1|data$Meet) + (1|data$University), data=data)
pred1.pm <- crosspred(cbnest1, model=lm1, cen= quantile_20, bylag = 1, at = start:end, cumul = TRUE)

exp = quantile_80
plot(pred1.pm, var = exp, xlab = 'Lag day', ylab = 'Time (seconds)', col="red",ylim=c(-3,10),
     main = "PM2.5")
plot(pred1.pm,"overall",xlab="PM2.5 (ug/m^3)", ylab = 'Time (seconds)', col="red",
     main="Overall cumulative association", ylim = c(-10, 35))
abline(v=12,col="black",lwd=2, lty = 2)

check = (exp - start) + 1
pred1.pm$cumfit[check,lag_selection]
pred1.pm$cumlow[check,lag_selection]
pred1.pm$cumhigh[check,lag_selection]


##Ozone
quantile(oz_daily1, c(0.20, 0.80))
quantile_20 = 36.9
quantile_80 = 54.9

start = 32.9
end = 59.4

lag_selection = 21
lag_df = 5
var_df = 5

oz_daily1 <- as.matrix(oz_daily[,rep(0:lag_selection,each=1)])
colnames(oz_daily1) <- paste("lag", 0:(lag_selection-1), sep="")
cbnest1 = crossbasis(oz_daily1, lag = c(0,(lag_selection-1)), argvar = list(fun="ns",df=var_df), arglag = list(fun = "ns", df = lag_df))
lm1 = lme4::lmer(data$adj_S ~ cbnest1 +
                   data$PR + data$real_days_year + 
                   data$days_since_event + data$farenheit + data$wind + data$spec_hum + (1|data$Meet) + (1|data$University), data=data)
pred1.oz <- crosspred(cbnest1, model=lm1, cen=quantile_20, bylag = 1, at = start:end, cumul = TRUE)

exp = quantile_80
plot(pred1.oz, var = exp,  xlab = 'Lag day', ylab = 'Time (seconds)', col="red", ylim = c(-3,10),
     main = "Ozone")

pred2.oz <- crosspred(cbnest1, model=lm1, cen=36.9, bylag = 1, at = 32.9:57.5, cumul = TRUE)
plot(pred2.oz, "overall",xlab="Ozone (ppm)", ylab = 'Time (seconds)', col="red",
     main="Overall cumulative association", ylim = c(-10, 35))
abline(v=54,col="black",lwd=2, lty = 2)


check = exp - start + 1 
pred1.oz$cumfit[check, lag_selection]
pred1.oz$cumlow[check, lag_selection]
pred1.oz$cumhigh[check, lag_selection]




