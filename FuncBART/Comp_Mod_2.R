
### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

library(tidyverse)
library(fda)

# 2 classes, 100 obervation per class, 200 obs total


# fUnction that creates exponentially distributed partition of unerlying T space

randomInterval_normal <- function(numInts, int_sd){
  tempInts= c(rep(0,numInts),1)
  for(i in 2:(numInts)){
    tempInts[i] = abs(rnorm(1, mean = (i-1)/numInts, sd=int_sd))
  }
  tempInts = sort(tempInts/max(tempInts))
  return(tempInts)
}



#Gets the AverageValue
scalarExtract <- function(intervals, coeffs, basis ){
  scalar_data = NULL
  for (i in 1:(length(intervals)-1)) {
    xvals = seq(by =0.01,from=intervals[i],to=intervals[i+1])
    yvals = t(coeffs)%*%t(eval.basis(xvals, basis))
    scalar_data = cbind(scalar_data,rowMeans(yvals) )
  }
  return(scalar_data)
}


#gets scalar extraction for M trees, each extraction is different interval
m_MV_ScalarExtraction <- function(m, numInts, fp, theta=1, coeffs, basis){
  m_scalar_data = array(dim = c(ncol(coeffs),numInts*fp,m))

  rInts = NULL
  ngrps = nrow(coeffs)/fp

  for (j in 1:m) {
    temp = randomInterval_normal(numInts,theta)
    print(temp)
    temp_scalar = NULL
    #For loop for each functional predictor
    for (k in 0:(fp-1)) {
      coef2 = coeffs[(1+ngrps*k):(ngrps*(k+1)), ]
      temp2= scalarExtract(temp, coef2, basis)
      temp_scalar=cbind(temp_scalar, temp2)
    }

    m_scalar_data[,,j]=temp_scalar
    rInts = rbind(rInts, temp)


  }
  #Saves the intervals and the scalar data
  ret=list()
  ret[[1]] = m_scalar_data
  ret[[2]] = rInts
  return(ret)
}


plotFunctionInterval <- function(data, numInts, intervals, theta=1, obs=190, ylab = "X(t)", xlab = "t" ){
  bb = create.bspline.basis(rangeval = c(0,0.2), norder=4)
  discrete_data = NULL

  plot(data, ylab = ylab, xlab = xlab)

  ints2= intervals

  for (i in 1:numInts) {
    xvals = seq(by =0.01,from=ints2[i],to=ints2[i+1])
    yvals = t(data$coefs)%*%t(eval.basis(xvals, data$basis))

    points(x= xvals, y=yvals[obs,], col="darkgreen", pch=16)
    meany = rep(mean(yvals[obs,]), length(xvals))
    discrete_data = cbind(discrete_data,rowMeans(yvals) )

    lines(x= xvals, y = meany, col="green", lwd=3)
  }

}

######################################################################

numInts=6

##########################################################
###########################################################
#testing on new functional dataset
load("~/Desktop/Functional Data/mfds-master/data/ArabicDigits.rda")



ArabicDigits$MFCC1[,1]
ArabicDigits$class[ArabicDigits$class==0]

grid = seq(.1, from=0, to=1)
numPoinst=length(grid)+2
bb <- create.bspline.basis(rangeval = c(0,1), breaks = grid, nbasis = numPoinst)

#zeros
arabic_zeros_C1 = ArabicDigits$MFCC1[,ArabicDigits$class==0]

f0 = Data2fd(arabic_zeros_C1, basisobj = bb)
plot(f0)

#ones
arabic_one_C1 = ArabicDigits$MFCC1[,ArabicDigits$class==1]

f1 = Data2fd(arabic_one_C1, basisobj = bb)
plot(f1)

#Channal 1

digits_c1=cbind(arabic_zeros_C1, arabic_one_C1 )
f_c1= Data2fd(digits_c1, basisobj = bb)

plot(f_c1)


#zeros
arabic_zeros_C2 = ArabicDigits$MFCC2[,ArabicDigits$class==0]

f0 = Data2fd(arabic_zeros_C2, basisobj = bb)
plot(f0)

#ones
arabic_one_C2 = ArabicDigits$MFCC2[,ArabicDigits$class==1]

f1 = Data2fd(arabic_one_C2, basisobj = bb)
plot(f1)

#Channal 2

digits_c2=cbind(arabic_zeros_C2, arabic_one_C2 )
f_c2= Data2fd(digits_c2, basisobj = bb)

plot(f_c2)



#zeros
arabic_zeros_C3 = ArabicDigits$MFCC3[,ArabicDigits$class==0]

f0 = Data2fd(arabic_zeros_C3, basisobj = bb)
plot(f0)

#ones
arabic_one_C3 = ArabicDigits$MFCC3[,ArabicDigits$class==1]

f1 = Data2fd(arabic_one_C3, basisobj = bb)
plot(f1)

#Channal 3

digits_c3=cbind(arabic_zeros_C3, arabic_one_C3 )
f_c3= Data2fd(digits_c3, basisobj = bb)

plot(f_c3)
f_c3$coefs

################################
#lines(mean.fd(f0), col="lightgrey", lwd = 5)
#lines(mean.fd(f1), col="lightgrey", lwd = 5)

# create a for loop to loop over all 13 channels of digits

sD2coefs= NULL

for (i in 1:13) {

  #zeros
  arabic_zeros = ArabicDigits[[i]][,ArabicDigits$class==0]
  f0 = Data2fd(arabic_zeros, basisobj = bb)
  #ones
  arabic_one = ArabicDigits[[i]][,ArabicDigits$class==1]
  f1 = Data2fd(arabic_one, basisobj = bb)


  #Channal 1
  digits=cbind(arabic_zeros, arabic_one)
  f_c1= Data2fd(digits, basisobj = bb)

  sD2coefs= rbind(sD2coefs, f_c1$coefs)

}



#####################################################################
######################################################################
#Compare other models using intervals
m=200
num_predictors=13
num_intervals = 6

s_data_test = m_MV_ScalarExtraction(m, num_intervals, num_predictors, 1/num_intervals^2, sD2coefs, f_c1$basis)


X=s_data_test[[1]][,,1]
s_data_test[[2]]

#Logistic regression

cc = c(ArabicDigits$class[ArabicDigits$class==0], ArabicDigits$class[ArabicDigits$class==1])

#create Train and Test sets
df= data.frame(interval=X, class = cc)
df = as_tibble(df)
df$class <- as.integer(df$class)-1



#use 50% of dataset as training set and 30% as test set
set.seed(1)
sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.05,0.95))
train  <- df[sample, ]
test   <- df[!sample, ]

#############################################################################
#############################################################################
############################################################################


#Logistic Regresion
library(glmnet)
m1=glm(class~., data = train, family = "binomial")
summary(m1)
m1$fitted.values

sd(m1$residuals)

p1 = predict(m1, test)
p1 = 1/(1+exp(-(p1)))


# Testing on BART
library(BART)
library(precrec)

m2 = pbart(x.train = as.matrix(train[,-(num_predictors*num_intervals+1)]), y.train = train$class,
           x.test = as.matrix(test[,-(num_predictors*num_intervals+1)]))

m2$prob.test.mean

m2$varcount.mean
barplot(m2$varcount.mean)
m2$varprob.mean



# TEsting using DART
m3 = pbart(x.train = as.matrix(train[,-(num_predictors*num_intervals+1)]), y.train = train$class,
           x.test = as.matrix(test[,-(num_predictors*num_intervals+1)]), sparse = TRUE)
m3$prob.test.mean

m3$varcount.mean
barplot(m3$varcount.mean)
m3$varprob.mean

#By function
matDART = matrix(m3$varprob.mean)
dim(matDART) = c(num_intervals, num_predictors)
matDART
barplot(colSums(matDART), main="functional predictor probability")
barplot(colSums(t(matDART)), main="Interval probability")



m3$proc.time


#Testing using Radnom Forest
library(randomForest)

m4 <- randomForest(class ~ ., data=train, importance=TRUE,
                   proximity=TRUE)


p4 = predict(m4, test)
###############################################################
###############################################################
#Adding functionalBart




###################################
###################################
#Calculate sigma value from linear regression
#Sigma_ols


#X = s_data_test[[1]]


X_train  <- s_data_test[[1]][sample, , ]

X_test   <- s_data_test[[1]][!sample, ,]

y_train <- df$class[sample]
y_test <-  df$class[!sample]

test$class-y_test

model_sigma_ols = lm(df$class~s_data_test[[1]][,,1])

summary(model_sigma_ols)
s_ols=sd(m1$residuals)


#find the degrees of freedom such that
nu=3
nu/qchisq(0.9,df=nu)

lambda = 0.9*pchisq(s_ols,df=nu)/nu

# inverse chi-square distribtuions
#lambda*nu/pchisq(s_ols,df=nu)

source("~/Desktop/testPack2/bartModelMatrix.R")
newX = list();
xInfo_list = list();
newX_test= list();
#Code to make it fit in C++ row dominant
for(i in 1:m){

  temp = bartModelMatrix(X_train[,,i], numcut=100, usequants=FALSE,
                         cont=FALSE, xinfo=matrix(0.0,0,0), rm.const=TRUE)
  xInfo_list[[i]] = temp$xinfo
  newX[[i]] = t(temp$X)

  newX_test[[i]] = t(X_test[,,i])

}


#PARAMETERS ################################################


funcP=num_predictors
intvls=num_intervals
p =num_predictors*num_intervals
obs_train=sum(sample)
obs_test=sum(!sample)
nc = rep(100, p)
ndpost=1000L
nskip=100L
keepevery=1L
nkeeptrain=ndpost
binaryOffset=qnorm(mean(y_train))
iter=1000
burn=500

#Priors
alpha=0.95
beta=2
tau=0.5


dart=TRUE
#DART_Functional_Predictors
dart_fp = TRUE
theta_fp = 0
omega_fp = 1
a_fp = 0.5
b_fp = 1
rho_fp = num_predictors
aug_fp = FALSE

#DART_Intervals
dart_int = TRUE
theta_int = 0
omega_int = 1
a_int = 1
b_int = 1
rho_int = num_intervals
aug_int = FALSE

lambda=4^2



Rcpp::sourceCpp("src/listTest3.cpp")


r1 = functionalBART_res(newX, #input matrix of scalar extraction data
                        y_train,    #y response vector
                        p,    #Number or intervals
                        obs_train,  #Observations
                        funcP,
                        intvls,
                        newX_test,
                        obs_test,

                        nu,
                        lambda,

                        xInfo_list,   #xInfoList
                        nc,            #Number of cut points
                        nkeeptrain,
                        binaryOffset=binaryOffset,
                        iter,
                        burn,
                        alpha,
                        beta,
                        tau,
                        #Functional Predictors
                        dart,

                        dart_fp,
                        theta_fp,
                        omega_fp,
                        a_fp,
                        b_fp,
                        rho_fp,
                        aug_fp,

                        dart_int,
                        theta_int,
                        omega_int,
                        a_int,
                        b_int,
                        rho_int,
                        aug_int
)


plot(r1$sigma_values, type ="l", main="Sigma value over iteration")
abline(v=burn, col="red")

r1$sigma_values[iter+burn]
mean(r1$sigma_values[burn:iter+burn])
#dev.off()

#probabilities

plot(pnorm(colSums(r1$yhats)), col="blue",
     main="As Probabilities")
#plots
plot(colSums(r1$yhats0), ylim=c(-2,2))
points(colSums(r1$yhats), col="green",
       main= "As values")




barplot(colSums(r1$vp_funcP)/sum(colSums(r1$vp_funcP)),
        main="functional predictor probability")



barplot(colSums(r1$vp_intvls)/sum(colSums(r1$vp_intvls)),
        main = "interval probability")



#Test values

plot(colSums(r1$yhats_test), col="red", ylim=c(-2,2))

y_out=0
ytest_out=0

prob_predictor = 0
prob_interval=0
for (i in 1:m) {
  y_out=y_out+colSums(r1$yhatsCube[,,i])/iter
  ytest_out=ytest_out+colSums(r1$yhatsCube_test[,,i])/iter

  prob_predictor= prob_predictor + colSums(r1$vp_funcP_Cube[,,i])/iter
  #print( colSums(r1$vp_funcP_Cube[,,i])/iter )
  prob_interval= prob_interval + colSums(r1$vp_intvls_Cube[,,i])/iter
}

#predictions
plot(pnorm(y_out), col="blue",
     main="Prediction As Probabilities")
plot(pnorm(ytest_out), col="red",
     main="FuncBART Test Values", ylab = "Probability")


###############################################################
###############################################################
# Functional Logistic Regression






###############################################################
###############################################################
#Compare all methods

CV1 = cbind(p1, m2$prob.test.mean,
            m3$prob.test.mean, p4, pnorm(ytest_out), test$class )

msmdat2 <- mmdata(CV1[,1:5],
                  c(CV1[,6]), modnames = c("LogReg", "BART", "DART", "RF", "funcBART" ))

mscurves <- evalmod(msmdat2)
par(mfrow=c(1,1))
autoplot(mscurves)
aa <- auc(mscurves)
aa$aucs

# Results
resultsDF = data.frame(matrix(aa$aucs, ncol=5))
names(resultsDF) = c("LogReg", "BART", "DART", "RF" , "funcBART")
resultsDF


par(mfrow=c(2,2))
barplot(prob_predictor/sum(prob_predictor), main="predictors FunctBART",
        names.arg  = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8","X9", "X10", "X11", "X12", "X13" ) , space=0)

barplot(prob_interval/sum(prob_interval), main="intervals FuncBART",
        names.arg  = c("1", "2", "3", "4", "5", "6") , space=0)

#To compare to DART
#By function

barplot(colSums(matDART), main="predictors DART",
        names.arg  = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8","X9", "X10", "X11", "X12", "X13") , space=0)
barplot(colSums(t(matDART)), main="Intervals DART",
        names.arg  = c("1", "2", "3", "4", "5", "6") , space=0)

par(mfrow=c(1,1))



