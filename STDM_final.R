install.packages(knitr)
install.packages('RCurl')
install.packages("forecast")
install.packages("ggfortify")
install.packages('nnet')
library(nnet)
library(forecast)
library(ggfortify)
library(sp)
library(sf)
library(ggplot2)
library(raster)
library(dplyr)
library(tmaptools)
library(here)
library(maptools)
library(spdep)
library(rgdal)
library(tmap)
library(gridExtra)
library(gstat)
library(OpenStreetMap)
library(spacetime)
library(knitr)
library(RCurl)
library(reshape)
library(lattice)
source("starima_package.R")

#import the data
rawdata <- read.csv("https://raw.githubusercontent.com/YaleGu/STDM_final/main/hungary_chickenpox.csv")
data <- rawdata[1:(length(rawdata)-1)]
data$Date <-  as.Date(data$Date,"%d/%m/%Y")

#check the distribution
hist(rawdata$Total, 20, main='Frequency histogram of total cases', xlab='cases')
abline(v=mean(rawdata$Total),col='red')

data_matrix<-data.matrix(data[,3:ncol(data)])
hist(data_matrix)
abline(v=mean(data_matrix), col="red")

#data overview
plot(rowMeans(data[,3:ncol(data)]), xlab = "Date", ylab = "Case", type="l", axat="n")

#data overview(county-level)
data2 <- data.frame(t(data))
colnames(data2)[1:ncol(data2)] <- data2[2,]
data2 <- data2[-c(1:2),]
newdata <- melt(data, id.vars = 1:2, measure.vars = 3:ncol(data))
colnames(newdata)[3:4] <- c("county", "cases") 


newdata$Date <-as.Date(newdata$Date,"%Y-%m-%d")
xyplot(cases ~ Date | county, xlab = newdata$Date, type='l',
       layout = c(5, 4),
       data=newdata,
       main = "Cases in Hungary")




#import the map and merge attributes
map <- st_read(here::here("map.shp"))
data2$county <- rownames(data2)
map2 <- merge(map, data2,
      by.x = "NAME_1", by.y = "county")

#plot the map
tm_shape(map2)+ 
  tm_fill("2005-01-03", style="jenks", palette="Purples")+
  tm_borders("white")+
  tm_compass(position=c("left","top"))+
  tm_scale_bar()


#temporal autocorrelation
acf(rawdata$Total, lag.max=80)
pacf(rawdata$Total, lag.max=80)

#spatial autocorrelation
W <- nb2listw(poly2nb(map2))
W2 <- listw2mat(W)


map3 <- as_Spatial(map2)
data_matrix2<-data.matrix(map3@data)

hungary_cases_avg <- rowMeans(data_matrix2[,2:ncol(data_matrix2)])
moran.test(x=hungary_cases_avg, listw=W)
moran.mc(x=hungary_cases_avg, listw=W, nsim=9999)
lm <- localmoran(x=hungary_cases_avg, listw=W)

#spatio-temporal autocorrelation
Wmat <- listw2mat(W)
stacf(data_matrix, Wmat, 50)

stpacf(data_matrix, Wmat, 50)

#merge the clustering results and plot
clustering <- read.csv("https://raw.githubusercontent.com/YaleGu/STDM_final/main/clustering%20result.csv")
map4 <- merge(map, clustering,
              by.x = "NAME_1", by.y = "X")

map4$cluster <- as.factor(map4$cluster,listw=W)

ggplot(map4) + 
  geom_sf(aes(fill=cluster))+
  scale_fill_brewer(palette = "Set3") + 
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


#ARIMA
decom <- stl(rawdata$Total,t.window=5, s.window="periodic")
autoplot(decom)

p1 <- autoplot(acf(rawdata$Total,plot=FALSE))
p2 <- autoplot(pacf(rawdata$Total,plot=FALSE))
p3 <- autoplot(acf(diff(rawdata$Total),plot=FALSE)) 
p4 <- autoplot(pacf(diff(rawdata$Total),plot=FALSE))
grid.arrange(p1,p2, p3, p4)

p1 <- autoplot(rawdata$Total)
p2 <- autoplot(ts(diff(rawdata$Total)))
p3 <- autoplot(ts(diff(rawdata$Total, lag=2)))
grid.arrange(p1,p2,p3)


fit <- auto.arima(rawdata[1:400,"Total"],lambda = "auto")
summary(fit)
checkresiduals(fit)
fit %>% forecast(h=50) %>% autoplot()

NRMSE_fit <- NRMSE(res=fit$residuals, obs=rawdata$Total)
tsdiag(fit)

plot(pre.ar$pred,type='line',col='red')
par(new=TRUE)
plot(rawdata[401:522,"Total"],type='line',col='blue')


pre.ar<-predict(fit, n.ahead=122)
pre.Ar <- Arima(rawdata[401:522,"Total"], model=fit)
matplot(cbind(pre.Ar$fitted, pre.Ar$x), type="l")

#ANN model
a <- as.matrix(rawdata$Total)
b <- as.matrix(a[-1,])
c <- as.matrix(b[401:521,])
case.nnet <- nnet(a[1:400,1], b[1:400,1], decay=5e-6, linout = TRUE, size=6)
case.pred<-predict(case.nnet, c)
matplot(cbind(c, case.pred),ylab="cases", xlab="Time (in weeks)", type="l")
RMSE_ANN <- sqrt(mean((c - case.pred)^2))
