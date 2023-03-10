# For the purpose of using Multivariate syntheic Functional Data

### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))
library(FADPclust)
library(fda)
library(tidyverse)

data("simData2") # 2 classes, 100 obervation per class, 200 obs total


class_1 = simData2[[1]]$coefs[,1:100]
#class_1$basis = simData2[[1]]$basis

class_2 = simData2[[1]]$coefs[,101:200]


#True signal of class!




#Generate a 1000 for each class

truth_1 = c(3,3,2,0,3,3,3,0,0,1)
#plot(truth_1)

truth_2 = c(1,-3,2,0,3,3,3,0,0,1)
#plot(truth_2)

c1_s_data = NULL
c2_s_data = NULL
for(i in 1:1000){
  c1_s_data = cbind(c1_s_data, truth_1 + rnorm(10,sd=1/2) )
  c2_s_data = cbind(c2_s_data, truth_2 + rnorm(10,sd=1/2) )
}
sim_data_1 = cbind(c1_s_data, c2_s_data)

####################
# Simulated Data
bb = create.bspline.basis(rangeval = c(0,1), norder=4)

par(mfrow=c(2,2))
plot(Data2fd(c1_s_data, basis=bb ))
plot(Data2fd(c2_s_data, basis=bb ))


truth_1 = c(1,3,2,0,3,3,3,0,0,-1)
#plot(truth_1)

truth_2 = c(1,3,2,0,3,-3,3,0,0,1)
#plot(truth_2)

c1_s_data = NULL
c2_s_data = NULL
for(i in 1:1000){
  c1_s_data = cbind(c1_s_data, truth_1 + rnorm(10,sd=.5) )
  c2_s_data = cbind(c2_s_data, truth_2 + rnorm(10,sd=.5) )
}
sim_data_2 = cbind(c1_s_data, c2_s_data)
############
# See simulated data
plot(Data2fd(c1_s_data, basis=bb ))
plot(Data2fd(c2_s_data, basis=bb ))


par(mfrow=c(1,1))



################ Build Model

m=200
num_predictors=8
num_intervals = 6


sD2coefs= rbind(sim_data_1, sim_data_2)

for(i in 1:(num_predictors-2)){
  noise1 = rnorm(prod(dim(sim_data_1)), sd=10)
  dim(noise1) = dim(sim_data_1)

  sD2coefs= rbind(sD2coefs, noise1)
}


s_data_test = m_MV_ScalarExtraction(m, num_intervals, num_predictors,
                                    1/num_intervals^2, sD2coefs, simData2[[3]]$basis)


colors = rainbow(num_intervals)
plot(rep(1,num_intervals+1) ~s_data_test[[2]][1,], cex=0,  ylim=c(0,21),
     ylab="Tree Number", xlab="Intervals on (t) grid",
     main="Intervals for each Tree for scalar extraction")
for(i in 1:20){
  for(j in 1:num_intervals){
    points(rep(i,2) ~s_data_test[[2]][i,j:(j+1)] , pch=1, lwd=2, cex=1.5)
    lines(rep(i,2) ~s_data_test[[2]][i,j:(j+1)], lwd=5.5, col=colors[j])
  }

}






X=s_data_test[[1]][,,1]



#Logistic regression

cc = c(rep(0,1000), rep(1,1000))

#create Train and Test sets
df= data.frame(interval=X, class = cc)
df = as_tibble(df)



#use 50% of dataset as training set and 30% as test set
set.seed(1)
sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.1,0.9))
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

m3$prob.test
m3$varcount.mean
barplot(m3$varcount.mean)
m3$varprob.mean
barplot(m3$varprob.mean)

m3$proc.time



#By function
matDART = matrix(m3$varprob.mean)
dim(matDART) = c(num_intervals, num_predictors)
matDART
barplot(colSums(matDART), main="functional predictor probability")
barplot(colSums(t(matDART)), main="Interval probability")


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
s_data_test[[1]][sample, ,]

X_train  <- s_data_test[[1]][sample, , ]
X_test   <- s_data_test[[1]][!sample, ,]

y_train <- cc[sample]
y_test <-  cc[!sample]

test$class-y_test

model_sigma_ols = lm(cc~s_data_test[[1]][,,1])

summary(model_sigma_ols)
s_ols=sd(model_sigma_ols$residuals)


#find the degrees of freedom such that
nu=3
nu/qchisq(0.9,df=nu)

lambda = 0.9*pchisq(s_ols,df=nu)/nu
1/lambda
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
burn=1000

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
a_int = 0.5
b_int = 1
rho_int = num_intervals
aug_int = FALSE

lambda=4^3



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



#r1$vp_funcP_Cube[,2,1]

#plot(r1$vp_funcP_Cube[,7,4])

par(mfrow=c(1,1))

#plot(r1$sigma_values, type ="l", main="Sigma value over iteration")
#abline(v=burn, col="red")
#r1$sigma_values[iter+burn]

y_out=0
ytest_out=0

prob_predictor = 0
prob_interval=0
for (i in 1:m) {
  y_out=y_out+colMeans(r1$yhatsCube[,,i])
  ytest_out=ytest_out+colMeans(r1$yhatsCube_test[,,i])

  prob_predictor= prob_predictor + colMeans(r1$vp_funcP_Cube[,,i])
  #print( colSums(r1$vp_funcP_Cube[,,i])/iter )
  prob_interval= prob_interval + colMeans(r1$vp_intvls_Cube[,,i])
}



# See how many times the tree chose that interval

#Convert to 2D matrix

#predictions
plot(pnorm(y_out), col="blue",
     main="Prediction As Probabilities")
plot(pnorm(ytest_out), col="red",
     main="Test valuesAs Probabilities")



#dev.off()

#probabilities



#barplot(colMeans(r1$vp_funcP) ,
  #      main="functional predictor probability")


#barplot(colMeans(r1$vp_intvls) ,
#        main = "interval probability")



#Test values

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
        names.arg  = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8") , space=0)

barplot(prob_interval/sum(prob_interval), main="intervals FuncBART",
        names.arg  = as.character(1:num_intervals) , space=0)

#To compare to DART
#By function

barplot(colSums(matDART), main="predictors DART",
        names.arg  = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8") , space=0)
barplot(colSums(t(matDART)), main="Intervals DART",
        names.arg  = as.character(1:num_intervals) , space=0)

par(mfrow=c(1,1))


#####################################################################
#####################################################################
# Building importance time grid


t_grid = seq(from=0, to=1, by=1/100)
intervals= s_data_test[[2]]
prob_grid = rep(1, length(t_grid))

prob_model_intervals = colSums(r1$vp_intvls)/sum(colSums(r1$vp_intvls))

get_grid_probabilites <- function(t_grid, interval, probs){

  prob_grid = rep(1, length(t_grid))
  for(i in 1:length(t_grid)){

    for(j in 1:(length(interval)-1 ) ){

      if(t_grid[i]>= interval[j] & t_grid[i]<= interval[j+1]){
        prob_grid[i] = probs[j]
      }
    }
  }
  return(prob_grid)
}

#get_grid_probabilites(t_grid,intervals[3,], prob_model_intervals )


get_grid_probabilites_M <- function(t_grid, interval, probs){
  ave_prob_grid = NULL
  for(i in 1:nrow(interval)){
    ave_prob_grid =rbind( ave_prob_grid, get_grid_probabilites(t_grid,interval[i,], probs))
  }
  return(ave_prob_grid)
}


ave_probs =colMeans(get_grid_probabilites_M(t_grid,intervals, prob_model_intervals ))

plot(ave_probs, ylim=c(0,max(ave_probs)))
lines(ave_probs)
abline(h=1/num_intervals)



rowMeans(colMeans(r1$probcount))
barplot(rowMeans(colMeans(r1$probcount)))
barplot(rowMeans(colMeans(r1$varcount)))



round(runif(200, min=1, max=10))
