# For the purpose of using Multivariate syntheic Functional Data

### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))
library(FADPclust)
library(fda)
library(tidyverse)


#True signal of class!




#Generate a 1000 for each class

truth_1 = c(3,3,2,0,3,3,3,0,0,1)
#plot(truth_1)

truth_2 = c(3,-3,2,0,3,3,3,0,0,1)
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


truth_1 = c(1,3,2,0,3,3,3,0,0,1)
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
############################################



# fUnction that creates exponentially distributed partition of unerlying T space
randomInterval <- function(numInts, theta){
  ints = rexp(numInts,theta)
  ints = ints/(sum(ints))
  ints2 = cumsum(ints)
  ints2 = c(0, ints2)

  return(ints2)
}

#random interval creation that jitters breakpoints of intervals





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
    #print(temp)
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



#####################################################################
######################################################################
#Compare other models using intervals
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

par(mfrow=c(1,1))
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



barplot(m2$varcount.mean)




# TEsting using DART
m3 = pbart(x.train = as.matrix(train[,-(num_predictors*num_intervals+1)]), y.train = train$class,
           x.test = as.matrix(test[,-(num_predictors*num_intervals+1)]), sparse = TRUE)



barplot(m3$varcount.mean)





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
iter=500
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
a_int = 0.5
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







par(mfrow=c(1,1))

plot(r1$sigma_values, type ="l", main="Sigma value over iteration")
abline(v=burn, col="red")
r1$sigma_values[iter+burn]


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



# See how many times the tree chose that interval

#Convert to 2D matrix

#predictions
plot(pnorm(y_out), col="blue",
     main="Prediction As Probabilities")
plot(pnorm(ytest_out), col="red",
     main="Test valuesAs Probabilities")



#dev.off()

#probabilities



barplot(colSums(r1$vp_funcP)/sum(colSums(r1$vp_funcP)),
        main="functional predictor probability")


barplot(colSums(r1$vp_intvls)/sum(colSums(r1$vp_intvls)),
        main = "interval probability")



#Test values

###############################################################
###############################################################
# Differing Intervals


#gets scalar extraction for M trees, each extraction is different interval
m_intervals_ScalarExtraction <- function(m, numInts, fp, theta=1, coeffs, basis){
  m_scalar_data = array(dim = c(ncol(coeffs),numInts*fp,m))

  rInts = NULL
  ngrps = nrow(coeffs)/fp

  for (j in 1:m) {
    temp = randomInterval_normal(numInts,1/numInts^2)
    #print(temp)
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



#random intervals per tree
m=200
num_interval = round(runif(m, min=4, max=10))

sD2coefs



intervals_list =list()
scalar_list = list()

for (i in 1:m) {

  ri_1 = randomInterval_normal(num_interval[i],0.1)
  temp = m_intervals_ScalarExtraction(1, num_interval[i], num_predictors, theta=.1, sD2coefs, simData2[[3]]$basis)
  scalar_list[[i]] = temp[[1]]
  intervals_list[[i]]= temp[[2]]
}



colors = rainbow(10)
seq(from=0, to=1, by=1/9)
plot(rep(1,10+1) ~ seq(from=0, to=1, by=1/10), cex=0,  ylim=c(0,20),
     ylab="Tree Number", xlab="Intervals on (t) grid",
     main="Intervals for each Tree for scalar extraction")
for(i in 1:20){
  for(j in 1:num_interval[i]){
    points(rep(i,2) ~ intervals_list[[i]][j:(j+1)] , pch=1, lwd=2, cex=1.5)
    lines(rep(i,2) ~intervals_list[[i]][j:(j+1)] , lwd=5.5, col=colors[j])
  }
}
#######################################################################

X_train = list()
X_test = list()

for (i in 1:m) {
  X_train[[i]] =  scalar_list[[i]][sample, , ]
  X_test[[i]] =  scalar_list[[i]][!sample,  ,]
}



source("~/Desktop/testPack2/bartModelMatrix.R")
newX = list();
xInfo_list = list();
newX_test= list();
#Code to make it fit in C++ row dominant
for(i in 1:m){

  temp = bartModelMatrix(as.matrix(X_train[[i]]), numcut=100, usequants=FALSE,
                         cont=FALSE, xinfo=matrix(0.0,0,0), rm.const=TRUE)
  xInfo_list[[i]] = temp$xinfo
  newX[[i]] = t(temp$X)

  newX_test[[i]] = t(as.matrix(X_test[[i]]))

}

#################################################################################
#################################################################################



funcP=num_predictors
max_intvls=10

interval_size=num_interval

p =num_predictors*max_intvls
obs_train=sum(sample)
obs_test=sum(!sample)
nc = rep(100, p)
ndpost=1000L
nskip=100L
keepevery=1L
nkeeptrain=ndpost
binaryOffset=qnorm(mean(y_train))
iter=500
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
a_int = 0.5
b_int = 1
rho_int = max_intvls
aug_int = FALSE

lambda=4^2

nu=3


Rcpp::sourceCpp("src/fbart_intervals.cpp")


r2 = functionalBART_intervals(newX, #input matrix of scalar extraction data
                              y_train,    #y response vector
                              p,    #Number or intervals
                              obs_train,  #Observations
                              funcP,
                              max_intvls,
                              interval_size,
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


y_out2 = rowMeans(colMeans(r2$yhatsCube))
ytest_out2 = rowMeans(colMeans(r2$yhatsCube_test))
prob_predictor2 = rowMeans(colMeans(r2$vp_funcP_Cube))


colMeans(r2$interval_prob[[20]])


par(mfrow=c(2,2))
barplot(colMeans(r2$interval_prob[[1]]),
        names.arg  = as.character(1:ncol(r2$interval_prob[[1]]) ) , space=0,
        main = "Tree 1 interval probabilities"
)

barplot(colMeans(r2$interval_prob[[2]]),
        names.arg  = as.character(1:ncol(r2$interval_prob[[2]]) ) , space=0,
        main = "Tree 2 interval probabilities"
)
barplot(colMeans(r2$interval_prob[[3]]),
        names.arg  = as.character(1:ncol(r2$interval_prob[[3]]) ) , space=0,
        main = "Tree 3 interval probabilities"
)

barplot(colMeans(r2$interval_prob[[20]]),
        names.arg  = as.character(1:ncol(r2$interval_prob[[20]]) ) , space=0,
        main = "Tree 20 interval probabilities"
)
par(mfrow=c(1,1))


#predictions
plot(pnorm(y_out2), col="blue",
     main="Prediction As Probabilities (INT)")
plot(pnorm(ytest_out2), col="red",
     main="Test valuesAs Probabilities (INT)")

barplot(prob_predictor2,
        main="functional predictor probability (INT)",
        names.arg  = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8") , space=0)

###### Get the probability grid
###################################################################
# reconstructing importance grid

t_grid = seq(from=0, to=1, by=1/100)
prob_grid = rep(1, length(t_grid))

intervals_list


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


get_grid_probabilites(t_grid,intervals_list[[1]], colMeans(r2$interval_prob[[1]]) )

probs_grid=rep(0,101)
for (i in 1:m) {

  probs_grid= probs_grid+get_grid_probabilites(t_grid,intervals_list[[i]], colMeans(r2$interval_prob[[i]]) )

}
probs_grid/m


plot(probs_grid/m ~ t_grid, ylim=c(0,max(probs_grid/m)),
     main = "Importance Score at (t)", ylab="Importance (Averaged Prob)")
lines(probs_grid/m ~t_grid)








###############################################################
###############################################################
#Compare all methods

CV1 = cbind(p1, m2$prob.test.mean,
            m3$prob.test.mean, p4, pnorm(ytest_out), pnorm(ytest_out2), test$class )

msmdat2 <- mmdata(CV1[,1:6],
                  c(CV1[,7]), modnames = c("LogReg", "BART", "DART", "RF", "funcBART", "funcBART_INT" ))

mscurves <- evalmod(msmdat2)
par(mfrow=c(1,1))
autoplot(mscurves)
aa <- auc(mscurves)
aa$aucs

# Results
resultsDF = data.frame(matrix(aa$aucs, ncol=6))
names(resultsDF) = c("LogReg", "BART", "DART", "RF" , "funcBART", "funcBART_INT" )
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

#get_grid_probabilites(t_grid, intervals[1,], colMeans(r1$vp_intvls_Cube[,,1]) )
#plot(get_grid_probabilites(t_grid, intervals[1,], colMeans(r1$vp_intvls_Cube[,,1]) ))




get_grid_probabilites_M <- function(t_grid, interval, probs){
  ave_prob_grid = NULL
  for(i in 1:nrow(interval)){
    ave_prob_grid =rbind( ave_prob_grid, get_grid_probabilites(t_grid,interval[i,], colMeans(probs[,,i])))
  }
  return(ave_prob_grid)
}

get_grid_probabilites_M(t_grid,intervals, r1$vp_intvls_Cube )


ave_probs = colMeans(get_grid_probabilites_M(t_grid,intervals, r1$vp_intvls_Cube ))

plot(ave_probs ~ t_grid, ylim=c(0,max(ave_probs)),
     main = "Importance Score at (t)", ylab="Importance (Averaged Prob)")
lines(ave_probs~ t_grid)
abline(h=1/num_intervals)


