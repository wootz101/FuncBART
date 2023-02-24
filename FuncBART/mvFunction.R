# Testing out the multivariate function FBART model
# Load in multivariate example



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
    temp = randomInterval(numInts,theta)
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



numInts=6
ints2= randomInterval(numInts, 1)

#See the DATA with intervals
par(mfrow=c(3,1))
plotFunctionInterval(simData2[[1]], 6, ints2, 1, 10, ylab ="X_1(t)")
plotFunctionInterval(simData2[[2]], 6, ints2, 1, 10, ylab ="X_2(t)")
plotFunctionInterval(simData2[[3]], 6, ints2, 1, 10, ylab ="X_3(t)")

# Rest to default
dev.off()
plotFunctionInterval(simData2[[1]], 6, ints2, 1, 10, ylab ="X_1(t)")
plotFunctionInterval(simData2[[2]], 6, ints2, 1, 10, ylab ="X_2(t)")
plotFunctionInterval(simData2[[3]], 6, ints2, 1, 10, ylab ="X_3(t)")

#####################################################################
#####################################################################

m=10

sD2coefs= rbind(simData2[[1]]$coefs, simData2[[2]]$coefs, simData2[[3]]$coefs)

#The Noisy 3rd function with no value

noise1 = rnorm(prod(dim(simData2[[3]]$coefs)), sd=100)
dim(noise1) = dim(simData2[[3]]$coefs)

noise2=rnorm(prod(dim(simData2[[3]]$coefs)), sd=100)
dim(noise2) = dim(simData2[[3]]$coefs)

revData = matrix(rev(simData2[[1]]$coefs))
dim(revData) = dim(simData2[[3]]$coefs)


sD2coefs= rbind( noise1, simData2[[1]]$coefs, noise2)

#sD2coefs= rbind( noise1, revData, noise2)



s_data_test = m_MV_ScalarExtraction(m, 6, 3, 1, sD2coefs, simData2[[1]]$basis)
#s_data_test[[2]]

s_data_test[[1]][,,1]
###################################
###################################
#Calculate sigma value from linear regression
#Sigma_ols



y=c(rep(0,100), rep(1,100))
X = s_data_test[[1]]

model_sigma_ols = lm(y~s_data_test[[1]][,,1])

summary(model_sigma_ols)
s_ols=sd(model_sigma_ols$residuals)


#find the degrees of freedom such that
nu=10
nu/qchisq(0.9,df=nu)

lambda = 0.9*pchisq(s_ols,df=nu)/nu

# inverse chi-square distribtuions
lambda*nu/pchisq(s_ols,df=nu)

source("~/Desktop/testPack2/bartModelMatrix.R")
newX = list();
xInfo_list = list();
newX_test= list();
#Code to make it fit in C++ row dominant
for(i in 1:m){

  temp = bartModelMatrix(X[,,i], numcut=100, usequants=FALSE,
                         cont=FALSE, xinfo=matrix(0.0,0,0), rm.const=TRUE)
  xInfo_list[[i]] = temp$xinfo
  newX[[i]] = t(temp$X)

  newX_test[[i]] = t(temp$X)+rnorm(200, sd=6)

}


#PARAMETERS ################################################
################################################
################################################

funcP=3
intvls=6
p =18
obs=200
nc = rep(100, p)
ndpost=1000L
nskip=100L
keepevery=1L
nkeeptrain=ndpost
binaryOffset=qnorm(mean(y))
iter=200
burn=200

#Priors
alpha=0.95
beta=2
tau=0.5


#DART
dart = TRUE
theta = 0
omega = 1
a = 0.5
b = 1
rho = p
aug = TRUE


#lambda=lambda



Rcpp::sourceCpp("src/listTest3.cpp")


r1 = functionalBART_res(newX, #input matrix of scalar extraction data
               y,    #y response vector
               p,    #Number or intervals
               obs,  #Observations
               funcP,
               intvls,
               newX_test,
               obs,

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
               dart,
               theta,
               omega,
               a,
               b,
               rho,
               aug
)


r1$sigma_values

plot(r1$sigma_values, type ="l", main="Sigma value over iteration")
abline(v=burn, col="red")
#dev.off()
r1$yhats
colSums(r1$yhats)
#probabilities
pnorm(colSums(r1$yhats))
plot(pnorm(colSums(r1$yhats)), col="blue",
     main="As Probabilities")
#plots
plot(colSums(r1$yhats0), ylim=c(-2,2))
points(colSums(r1$yhats), col="green",
       main= "As values")



r1$vc_funcP
r1$vp_funcP

barplot(colSums(r1$vp_funcP)/sum(colSums(r1$vp_funcP)),
        main="functional predictor probability")


r1$vc_intvls
r1$vp_intvls
barplot(colSums(r1$vp_intvls)/sum(colSums(r1$vp_intvls)),
        main = "interval probability")



#Test values
r1$yhats_test
colSums(r1$yhats_test)
plot(colSums(r1$yhats0), ylim=c(-2,2))
points(colSums(r1$yhats_test), col="red")

plot(pnorm(colSums(r1$yhats_test)), col="blue",
     main="As Probabilities TEST")






r1$probcount




# See how many times the tree chose that interval
r1$varcount
#Convert to 2D matrix
dim(r1$varcount) <- c(iter*m, p)

colSums(r1$varcount)
colSums(r1$varcount)/sum(colSums(r1$varcount))

barplot(colSums(r1$varcount)/sum(colSums(r1$varcount)),
        main="Variable Count percentage",
        names.arg  = c(rep("X_1", 6),rep("X_2", 6),rep("X_3", 6)), space=0)

allV = colSums(r1$varcount)/sum(colSums(r1$varcount))
allV = matrix(allV, nrow=6, ncol=3)
allV = t(allV)
allV
colSums(allV)

barplot( colSums(allV),
        main="Variable Count Only (t)", xlab="interval")


barplot( rowSums(allV),
         main="Variable Count Functional Predictors", xlab="X_i(t)")


# See the average probabilites from m trees
r1$probcount

dim(r1$probcount) <- c(iter*m, p)

colSums(r1$probcount)
colSums(r1$probcount)/sum(colSums(r1$probcount))
barplot(colSums(r1$probcount)/sum(colSums(r1$probcount)),
        main="Probability Count percentage")

allP = colSums(r1$probcount)/sum(colSums(r1$probcount))
allP = matrix(allP, nrow=6, ncol=3)
allP = t(allP)
allP
colSums(allP)

barplot( colSums(allP),
         main="Probability Only (t)", xlab="interval")

barplot( rowSums(allP),
         main="Probability Functional Predictors", xlab="X_i(t)")


















