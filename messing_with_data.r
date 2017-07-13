library(readr)



train <- read_csv("C:/Users/Evan Chisholm/Desktop/CMPT318/train.txt", 
                  col_types = cols(Date = col_date(format = "%d/%m/%Y"), 
                                   Time = col_time(format = "%H:%M:%S")
                    ))
View(train)

x <- train$Global_active_power[0:10000]
mean(x,na.rm=TRUE)
plot(x)
#acf(x, lag.max = NULL, type = "covariance", plot= TRUE, na.action=na.pass)
