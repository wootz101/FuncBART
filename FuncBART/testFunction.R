
### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))
#setwd("~/Desktop/Functional Data")



library(FADPclust)
library(fda)
library(tidyverse)

data(simData1)
#simData1

plot(simData1)

bb = create.bspline.basis(rangeval = c(0,0.2), norder=4)
#eval.basis(seq(by =0.01,from=0,to=0.2), simData1$basis)


#create the intervals
numInts = 6
ints = rexp(numInts,300)
ints = ints/(sum(ints))
ints

ints2 = cumsum(ints)
ints2 = c(0, ints2)

plot(simData1)
discrete_data = NULL


# fUnction that creates exponentially distributed partition of unerlying T space
randomInterval <- function(numInts, theta){
  ints = rexp(numInts,theta)
  ints = ints/(sum(ints))
  ints2 = cumsum(ints)
  ints2 = c(0, ints2)

  return(ints2)
}

randomInterval(8,1)

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
mScalarExtraction <- function(m, numInts, theta=1, coeffs, basis){
  m_scalar_data = array(dim = c(dim(coeffs)[2],numInts,m))
  rInts = NULL

  for (j in 1:m) {
    temp = randomInterval(numInts,theta)
    m_scalar_data[,,j]= scalarExtract(randomInterval(numInts,theta),
                                      coeffs, basis)
    rInts = rbind(rInts, temp)


  }
  #Saves the intervals and the scalar data
  ret=list()
  ret[[1]] = m_scalar_data
  ret[[2]] = rInts
  return(ret)
}

m=10

s_data_test = mScalarExtraction(m, 6, 1, simData1$coefs, simData1$basis)




plotFunctionInterval <- function(data, numInts, theta=1){
  bb = create.bspline.basis(rangeval = c(0,0.2), norder=4)

  plot(data, ylab = "X(t)", xlab = "t")

  ints2= randomInterval(numInts, theta)
  for (i in 1:numInts) {
    xvals = seq(by =0.01,from=ints2[i],to=ints2[i+1])
    yvals = t(data$coefs)%*%t(eval.basis(xvals, data$basis))

    points(x= xvals, y=yvals[190,], col="darkgreen", pch=16)
    meany = rep(mean(yvals[190,]), length(xvals))
    discrete_data = cbind(discrete_data,rowMeans(yvals) )

    lines(x= xvals, y = meany, col="green", lwd=3)
  }

}

plotFunctionInterval(simData1, 6, 1)



plot(simData1, ylab = "X(t)", xlab = "t")

ints2= randomInterval(numInts, 1)
for (i in 1:numInts) {
  xvals = seq(by =0.01,from=ints2[i],to=ints2[i+1])
  yvals = t(simData1$coefs)%*%t(eval.basis(xvals, simData1$basis))

  points(x= xvals, y=yvals[190,], col="darkgreen", pch=16)
  meany = rep(mean(yvals[190,]), length(xvals))
  discrete_data = cbind(discrete_data,rowMeans(yvals) )

  lines(x= xvals, y = meany, col="green", lwd=3)
}





###############################################################
###############################################################

#Test dataset for functionalbart


y=c(rep(0,100), rep(1,100))
X = s_data_test[[1]]


source("~/Desktop/testPack2/bartModelMatrix.R")
newX = list();
xInfo_list = list();
newX_test= list();
for(i in 1:m){


  temp = bartModelMatrix(X[,,i], numcut=100, usequants=FALSE,
                         cont=FALSE, xinfo=matrix(0.0,0,0), rm.const=TRUE)
  xInfo_list[[i]] = temp$xinfo
  newX[[i]] = t(temp$X)

  newX_test[[i]] = t(temp$X)+rnorm(200, sd=6)

}


#PARAMETERS

nc = rep(100, 6)
ndpost=1000L
nskip=100L
keepevery=1L
nkeeptrain=ndpost
binaryOffset=qnorm(mean(y))
iter=100

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
rho = 6
aug = TRUE






r1 = listTest3(newX, #input matrix of scalar extraction data
               y,    #y response vector
               #newX_test,
               6,    #Number or intervals
               200,  #Observations
               newX_test,
               200,
               xInfo_list,   #xInfoList
               nc,            #Number of cut points
               nkeeptrain,
               binaryOffset=binaryOffset,
               iter,
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



r1$yhats
colSums(r1$yhats)
#probabilities
pnorm(colSums(r1$yhats))
#plots
plot(colSums(r1$yhats0), ylim=c(-m,m))
points(colSums(r1$yhats), col="green")


# See how many times the tree chose that interval
r1$varcount
colSums(r1$varcount)
colSums(r1$varcount)/sum(colSums(r1$varcount))
barplot(colSums(r1$varcount)/sum(colSums(r1$varcount)),
        main="Variable Count percentage")

# See the average probabilites from m trees
r1$probcount
colSums(r1$probcount)
colSums(r1$probcount)/sum(colSums(r1$probcount))
barplot(colSums(r1$probcount)/sum(colSums(r1$probcount)),
        main="Probability Count percentage")

#Test values
r1$yhats_test
colSums(r1$yhats_test)
plot(colSums(r1$yhats0), ylim=c(-m,m))
points(colSums(r1$yhats_test), col="red")






dim(r1$varcount)
r1$varcount[1,,]



sumVC=r1$varcount[1,,]
for(i in 2:iter){

  sumVC= sumVC + colSums(r1$varcount[i,,])

}

sum(colSums(sumVC))
