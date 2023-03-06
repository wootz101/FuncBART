# Creating differen tnumbe rof intervals per tree


### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))
library(FADPclust)
library(fda)
library(tidyverse)

data("simData2") # 2 classes, 100 obervation per class, 200 obs total


# fUnction that creates exponentially distributed partition of unerlying T space
randomInterval <- function(numInts, theta){
  ints = rexp(numInts,theta)
  ints = ints/(sum(ints))
  ints2 = cumsum(ints)
  ints2 = c(0, ints2)

  return(ints2)
}

#random interval creation that jitters breakpoints of intervals

numInts=6



randomInterval_normal <- function(numInts, int_sd){
  tempInts= c(rep(0,numInts),1)
  for(i in 2:(numInts)){
    tempInts[i] = abs(rnorm(1, mean = (i-1)/numInts, sd=int_sd))
  }
  tempInts = sort(tempInts/max(tempInts))
  if(tempInts[1]!=0){
    tempInts[1]=0
  }
  return(tempInts)
}

randomInterval_normal(numInts, .1)



#Gets the AverageValue
scalarExtract <- function(intervals, coeffs, basis ){
  scalar_data = NULL
  for (i in 1:(length(intervals)-1)) {
    xvals = seq(by =0.01,from=intervals[i],to=intervals[i+1])
    yvals = t(coeffs)%*%t(eval.basis(xvals, basis))
    scalar_data = cbind(scalar_data, rowMeans(yvals) )
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

numInts=6
ints2= randomInterval_normal(numInts, .1)

#See the DATA with intervals
par(mfrow=c(1,3))
plotFunctionInterval(simData2[[1]], 6, ints2, 1, 10, ylab ="X_1(t)")
plotFunctionInterval(simData2[[2]], 6, ints2, 1, 10, ylab ="X_2(t)")
plotFunctionInterval(simData2[[3]], 6, ints2, 1, 10, ylab ="X_3(t)")

# Rest to default
dev.off()
plotFunctionInterval(simData2[[1]], 6, ints2, 1, 10, ylab ="X_1(t)")
plotFunctionInterval(simData2[[2]], 6, ints2, 1, 10, ylab ="X_2(t)")
plotFunctionInterval(simData2[[3]], 6, ints2, 1, 10, ylab ="X_3(t)")

#####################################################################
#random intervals per tree
m=20
num_interval = round(runif(m, min=3, max=10))

num_predictors=3

noise1 = rnorm(prod(dim(simData2[[3]]$coefs)), sd=10)
dim(noise1) = dim(simData2[[3]]$coefs)

noise2 = rnorm(prod(dim(simData2[[3]]$coefs)), sd=10)
dim(noise2) = dim(simData2[[3]]$coefs)

#sD2coefs= rbind(simData2[[1]]$coefs, simData2[[2]]$coefs, simData2[[3]]$coefs)
sD2coefs= rbind(simData2[[3]]$coefs, noise1, noise2)




intervals_list =list()
scalar_list = list()

for (i in 1:m) {

  ri_1 = randomInterval_normal(num_interval[i],0.1)
  temp = m_MV_ScalarExtraction(1, num_interval[i], 3, theta=.1, sD2coefs, simData2[[3]]$basis)
  scalar_list[[i]] = temp[[1]]
  intervals_list[[i]]= temp[[2]]

}



colors = rainbow(15)
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

#################################################################################
#################################################################################




source("~/Desktop/testPack2/bartModelMatrix.R")
newX = list();
xInfo_list = list();
newX_test= list();
#Code to make it fit in C++ row dominant
for(i in 1:m){

  temp = bartModelMatrix(as.matrix(scalar_list[[i]]), numcut=100, usequants=FALSE,
                         cont=FALSE, xinfo=matrix(0.0,0,0), rm.const=TRUE)
  xInfo_list[[i]] = temp$xinfo
  newX[[i]] = t(temp$X)

  newX_test[[i]] = t(as.matrix(scalar_list[[i]]))

}
cc = c(rep(0,100), rep(1,100))
y_train <- cc

#PARAMETERS ################################################


funcP=num_predictors
intvls=10

interval_size=num_interval

p =num_predictors*10
obs_train=200
obs_test=200
nc = rep(100, p)
ndpost=1000L
nskip=100L
keepevery=1L
nkeeptrain=ndpost
binaryOffset=qnorm(mean(y_train))
iter=500
burn=100

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
rho_int = 3
aug_int = FALSE

lambda=4^2

nu=3


Rcpp::sourceCpp("src/fbart_intervals.cpp")


r1 = functionalBART_intervals(newX, #input matrix of scalar extraction data
                        y_train,    #y response vector
                        p,    #Number or intervals
                        obs_train,  #Observations
                        funcP,
                        10,
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

r1$yhatsCube
y_out=0
ytest_out=0

prob_predictor = 0
prob_interval=0
for (i in 1:m) {
  y_out=y_out+colSums(r1$yhatsCube[,,i])/iter
  ytest_out=ytest_out+colSums(r1$yhatsCube_test[,,i])/iter

  prob_predictor= prob_predictor + colSums(r1$vp_funcP_Cube[,,i])/iter
  #print( colSums(r1$vp_funcP_Cube[,,i])/iter )
  #prob_interval= prob_interval + colSums(r1$vp_intvls_Cube[,,i])/iter
}

plot(y_out)
plot(ytest_out)

r1$interval_prob[[1]]
colMeans(r1$interval_prob[[1]])

barplot(colMeans(r1$interval_prob[[1]]))

barplot(colMeans(r1$interval_prob[[2]]))
barplot(colMeans(r1$interval_prob[[3]]))
barplot(colMeans(r1$interval_prob[[4]]))

r1$vp_funcP_Cube
barplot(prob_predictor)

plot(colMeans(r1$vp_funcP_Cube))

rowMeans(colMeans(r1$vp_funcP_Cube))


###################################################################
###################################################################
# reconstructing importance grid

t_grid = seq(from=0, to=1, by=1/100)
prob_grid = rep(1, length(t_grid))

intervals_list
colMeans(r1$interval_prob[[2]])

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


get_grid_probabilites(t_grid,intervals_list[[1]], colMeans(r1$interval_prob[[1]]) )

probs_grid=rep(0,101)
for (i in 1:m) {

  probs_grid= probs_grid+get_grid_probabilites(t_grid,intervals_list[[i]], colMeans(r1$interval_prob[[i]]) )

}
probs_grid/m






plot(probs_grid/m, ylim=c(0,max(probs_grid/m)))
lines(probs_grid/m)



