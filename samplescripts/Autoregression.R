
dailyts <- ts(dayeven$flow, start=1, frequency = 365) 
#decompose the dataset
components.ts = decompose(dailyts)
x=1:1417
par(mfrow=c(4,1))
plot(x,dailyts,type="l", ylab="Daily Time series")
plot(x,components.ts$seasonal,type="l",col="blue", ylab="Seasonality")
plot(x,components.ts$trend,type="l", col="green", ylab="Trend")
plot(x,components.ts$random,type="l", col="red", ylab="Random")


plot(components.ts, col =c("blue","green","black","red"))
acf(dailyts  )
pacf(dailyts  )
# go with finding the difference of the data and do analysis

second.comp = decompose(components.ts$random)
plot(second.comp)


hourlyts <- ts(houreven$flow, start=1, frequency = 365) 
#decompose the dataset
first.h.ts = decompose(hourlyts)
plot(first.h.ts)
# go with finding the difference of the data and do analysis

second.h.comp = decompose(first.h.ts$random)
plot(second.h.comp)

third.h.comp = decompose(second.h.comp$random)
plot(third.h.comp)

fourth.h.comp = decompose(third.h.comp$random)
plot(fourth.h.comp)
#fit a curve to data to model the residual and remove the trend


#try unit root test: This test is used to find out that first difference or regression which should be used on the trending data to make it stationary. In Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test, small p-values suggest differencing is required.
library("fUnitRoots")
urkpssTest(dailyts, type = c("tau"), lags = c("long"),use.lag = NULL, doplot = TRUE)
tsstationary = diff(dailyts , differences=100)
plot(tsstationary)



urkpssTest(hourlyts, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)
h.tsstationary = diff(hourlyts , differences=1000)
plot(h.tsstationary)

