library(psych)
pairs.panels(train[3:4])

plot(train$Global_active_power[1:3000])
library(PerformanceAnalytics)
chart.Correlation(train[3:4])

library(corrplot)
x <- train[3:4]
corrplot(x, type="upper", order="hclust")

######


library("fUnitRoots")
library(tseries)
#urkpssTest(x_series, type = c("tau"), lags = c("long"),use.lag = NULL, doplot = TRUE)
tsstationary = diff(x_series , differences=100)
plot(tsstationary)


library(MASS)
require(fitdistrplus)
set.seed(1)
dat <- rnorm(50,0,1)
descdist(dat, discrete = FALSE)
f1 <- fitdist(dat,"norm")
f2 <- fitdist(dat,"logis")
f3 <- fitdist(dat,"cauchy")
plotdist(dat,"norm",para=list(mean=f1$estimate[1],sd=f1$estimate[2]))