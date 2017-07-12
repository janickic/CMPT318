###
#train data
traintest1 <- list (x= data.frame(uni_train$power), N=length(uni_train$Date))
traintest1$x <- data.matrix(traintest1$x, rownames.force = NA)
traintest1$x <- scale(traintest1$x)
class(traintest1) <- "hsmm.data"
# end
uni_test<- data.frame(as.numeric(test1$Global_active_power))
uni_test <- setNames(uni_test,c("power"))
uni_test$power[is.na(uni_test$power)] <- mean(na.omit(uni_test$power))
xtest1 <- list (x= data.frame(uni_test$power), N=length(uni_test$power))
xtest1$x <- data.matrix(xtest1$x, rownames.force = NA)
xtest1$x <- scale(xtest1$x)
class(xtest1) <- "hsmm.data"
# create M-step



library(mhsmm)
J <- 2
init0 <- rep(1/J, J)
B0 <- list(mu = c(0.5,3), sigma = c(2.5, 1.5))

P <- matrix(c(0, 1/2,
              1/2,0
                     ), nrow = J)

M <- length(uni_train$power)
d0 <- cbind(dunif(1:M, 1000, 8000), dunif(1:M, 1000, 8000))


startval <- hsmmspec(init0, P, parms.emis = B0, list(d = d0, type = "nonparametric"), dens.emis = dnorm.hsmm)


nnn <- traintest1$x[,1]
UniHSMM<- hsmmfit(nnn, startval, mstep = mstep.norm, M=1000)
#summary(hmv)
UniHSMMtrain <- predict.hsmm(UniHSMM,traintest1$x)
UniHSMMtest <- predict.hsmm(UniHSMM,xtest1$x)
plot(UniHSMMtrain)
#plot(UniHSMMtest)
UniHSMMtrain$s
UniHSMMtest$s

