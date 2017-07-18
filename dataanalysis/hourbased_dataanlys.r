

library(readr)


train <- read_csv("/home/heather/Code/cmpt-318/train.txt", #switch filename 
                  col_types = cols(Date = col_date(format = "%d/%m/%Y"), 
                                   Time = col_time(format = "%H:%M:%S")
                    ))

active_powers = numeric(25) 
reactive_powers = numeric(25) 
voltages=numeric(25)
intensities=numeric(25)
hours=1:25
for (i in hours){
  x<-subset(train,i*0<Time & Time<i*3600) #grab all entries from a particular hour
  
  active_power <- x$Global_active_power[0:10000]
  active_powers[i]=mean(active_power,na.rm=TRUE) #load the means into a vector
  
  reactive_power <- x$Global_reactive_power[0:10000]
  reactive_powers[i]=mean(reactive_power,na.rm=TRUE)
  
  voltage <- x$Voltage[0:10000]
  voltages[i]=mean(voltage,na.rm=TRUE)
  
  intensity <- x$Global_intensity[0:10000]
  intensities[i]=mean(intensity,na.rm=TRUE)
  
}

#plot vectors
plot(active_powers)
plot(reactive_powers)
plot(voltages)
plot(intensities)

cor(active_powers,reactive_powers)
cov(active_powers,reactive_powers)

cor(active_powers,voltages)
cov(active_powers,voltages)

cor(active_powers,intensities)
cov(active_powers,intensities)

cor(reactive_powers,voltages)
cov(reactive_powers,voltages)

cor(reactive_powers,intensities)
cov(reactive_powers,intensities)

cor(voltages,intensities)
cov(voltages,intensities)



#active powers and intensities strongly related

