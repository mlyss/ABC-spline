library(DEoptim)

#calculate pivot points from a column of some positive values 
toCumulative <- function (numbersColumnAsFrame) {
  baseFrame <- data.frame(x=1, y=numbersColumnAsFrame[,1],row.names = rownames(numbersColumnAsFrame))
  baseFrame <- baseFrame[order(baseFrame$y, decreasing = T),]
  sum<-apply(baseFrame, MARGIN =  2, FUN = sum ) 
  iRowCumulated<-rep(0, ncol(baseFrame))
  for(i in rownames(baseFrame)){
    iRowCumulated<-iRowCumulated+baseFrame[i,]/sum
    baseFrame[i,]<-iRowCumulated
  }
  baseFrame<-rbind(baseFrame, c(0,0))
  baseFrame <- baseFrame[order(baseFrame$x, decreasing = F),]
  
  for(i in 1:nrow(baseFrame)) {
    for(j in 1:ncol(baseFrame)){
      if(baseFrame[i,j]>1){
        baseFrame[i,j]<-1
      }
    }
  }
  
  return(baseFrame)
}


# calculate cumulative curve pivot points based sample of 500 from uni[0;1]
cumData <- toCumulative(data.frame("V1"=runif(500,0,1)))


#list of functions to approximate pivot points
FUN<-list(
  "1_Pareto"= function(x, a)x^a,
  "2_Kakwani_Podder"= function(x, a)(1-((1-x)^a[1])*exp(-a[2]*x)),
  "3_Rasche"= function(x, a)(1-(1-x^a[1])^a[2]),
  "4_Ortega"= function(x, a)(1-((1-x)^a[1])*(1-x^a[2])),
  "5_Aggarwal_Arnold"= function(x, a)(x*(2*a[1]+a[2]-1+(1-a[1])*x)/(a[1]+a[2]*x)),
  "6_Choticapanich"= function(x, a)(exp(a*x)-1)/(exp(a)-1),
  "7_Pearl_Reed"= function(x, a)( exp(1-x^a)-exp(1) )/( 1-exp(1) ),
  "8_Fourier"= function(x, a)sin((pi*x^a)/2),
  "9_Ballou"= function(x, a)((a+1)*x)/(a+x),
  "10_Ballou_Pareto_1"= function(x, a)((a[1]+1)*x^a[2])/(a[1]+x),
  "11_funBallou_Pareto_2"= function(x, a)((a[1]+1)*x)/(a[1]+x^a[2]),
  "12_funBallou_Pareto_3"= function(x, a)((a[1]+1)*x^a[2])/(a[1]+x^a[2]),
  "13_funBallou_Pareto_4"= function(x, a)(((a[1]+1)*x)/(a[1]+x))^a[2],
  "14_Fourier_Ballou"= function(x, a) sin((pi/2)*(((a[1]+1)*x)/(a[1]+x))^a[2]),
  "15_Radical"= function(x, a)(a[1]*x+(1-a[1])*x^2)^a[2]
)


# count error sum of squares for points in dataSet and function fun with parameters vector x
#if number of points = 0 ESS =0, if ESS is counted as NA or NaN, ESS becomes +infinitive
ESS<-function(x, fun, dataSet){
  returnESS<-0
  if(nrow(dataSet)>0){
    for(i in 1:nrow(dataSet)){
      returnESS<-returnESS+(fun(x=dataSet[i,1],a=x)-dataSet[i,2])^2
    }
  }
  if(is.na(returnESS)|is.nan(returnESS)) {
    returnESS<-Inf
  }
  
  return(returnESS)
}


#get starting values of parameters for optimization
optimParam<-function(funID){
  if(funID==1){return(0.53)}
  if(funID==2){return(c(2, 0.05))}
  if(funID==3){return(c(0.9, 2.1))}
  if(funID==4){return(c(1.1, 0.9))}
  if(funID==5){return(c(10, 20))}
  if(funID==6){return(-2)}
  if(funID==7){return(0.7)}
  if(funID==8){return(0.8)}
  if(funID==9){return(0.6)}
  if(funID==10){return(c(0.5, 1.1))}
  if(funID==11){return(c(1.1, 1.8))}
  if(funID==12){return(c(0.2, 1.4))}
  if(funID==13){return(c(0.17, 2))}
  if(funID==14){return(c(6, 0.9))}
  if(funID==15){return(c(2, 0.9))}
}

#spline function used in optimizing
#xABxBC_par1_par2_par3 is vector of optimized values: limits x(AB) and x(BC), parameters of each function
splineESS<-function(xABxBC_par1_par2_par3, fun1_id, fun2_id, fun3_id, cumData1){
  
  xAB<-xABxBC_par1_par2_par3[1]
  xBC<-xABxBC_par1_par2_par3[2]
  
  if(xAB>xBC|xAB<0|xBC<0|xAB>1|xBC>1){
    return(Inf)
    break()
  }
  
  par1_indexes<-3:(2+length(optimParam(fun1_id)))
  par2_indexes<-(max(par1_indexes)+1):(max(par1_indexes)+length(optimParam(fun2_id)))
  par3_indexes<-(max(par2_indexes)+1):(max(par2_indexes)+length(optimParam(fun3_id)))
  
  ESS1<-ESS(xABxBC_par1_par2_par3[par1_indexes], FUN[[fun1_id]], cumData1[cumData1[,1]<xAB,])
  ESS2<-ESS(xABxBC_par1_par2_par3[par2_indexes], FUN[[fun2_id]], cumData1[cumData1[,1]>=xAB & cumData1[,1]<xBC,])
  ESS3<-ESS(xABxBC_par1_par2_par3[par3_indexes], FUN[[fun3_id]], cumData1[cumData1[,1]>=xBC,])
  
  return(ESS1+ESS2+ESS3)
  
}



#the obtained result depends on the method of optimization
#some possible ways of optimization: using optim(), DEoptim()
#the case of (fun1_id=2, fun2_id=4, fun3_id=4) is one of complex ones


# "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH" ????
test1<-optim(par = c(0.3,0.6,optimParam(2),optimParam(4),optimParam(4)), fn = splineESS, method = "L-BFGS-B",
             fun1_id=2, fun2_id=4, fun3_id=4, cumData1=cumData)


# parameter limits have to be defined better
# is obtained result reliable?
test2<-DEoptim(fn=function(x){splineESS(xABxBC_par1_par2_par3 = x,fun1_id = 2, fun2_id = 4, fun3_id = 4, cumData1 = cumData)},
               lower=c(0,0,-99,-99,-99,-99,-99,-99), upper=c(1,1,99,99,99,99,99,99),control = DEoptim.control(itermax=30) )



#a data frame for optimization results using different functions combinations is created and filled
splineRes<-data.frame(matrix(ncol = 7, nrow = 0))
colnames(splineRes)<-c("fun1_id", "fun2_id", "fun3_id", "xAB", "xBC", "sumESS","funs_params")

for (AfunID in c(2,3,4,5,6,9,11,14,15)) {
  for (BfunID in c(2,3,4,5,6,9,11,14,15)){
    for (CfunID in c(2,3,4,5,6,9,11,14,15)){
      current_optim<-optim(par = c(0.3,0.6,optimParam(AfunID),optimParam(BfunID),optimParam(CfunID)),
                           fn = splineESS, method = "L-BFGS-B",
                           fun1_id=AfunID, fun2_id=BfunID, fun3_id=CfunID, cumData1=cumData)
      
      splineRes<-rbind(splineRes, data.frame("fun1_id"=AfunID,
                                             "fun2_id"=BfunID,
                                             "fun3_id"=CfunID,
                                             "xAB"=format(current_optim$par[1], scientific=F),
                                             "xBC"=format(current_optim$par[2], scientific=F),
                                             "sumESS"=current_optim$value,
                                             "funs_params"=paste( format(current_optim$par[-(1:2)], scientific=F), collapse = "; "  )
      ) )
      
      print(paste( c( as.character(Sys.time())," - ",AfunID, BfunID, CfunID ), collapse = "" ))
    }
  }
}



