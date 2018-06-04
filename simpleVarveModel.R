
simpleVarveModel = function(signal,H,shape=1.5,mean = 1,SNR=.25){
  #simple PSM for varves
  #signal = the climate forcing variable, at annual timescales, (e.g., a vector of temperature)
  #H = Hurst coefficient
  #shape=of the gamma distribution
  #mean of the varve series output
  
  #first, gammify the input
  gamSig = gammify(signal,shape=shape,mean = mean)
  
  
  #create gammafied autocorrelated fractal brownian motion series
  gamNoise = gammify(FGN::SimulateFGN(length(signal),H),shape = shape, mean = mean)
  
  #combine the signal with the noise, based on the SNR
  varves = (gamSig*SNR + gamNoise*(1/SNR))/(SNR+1/SNR)
  return(varves)
}

gammify <- function (X,shape = 1.5, mean = 1,jitter=FALSE){ 
  #   Transform each column of data matrix X to a gaussian distribution using the inverse Rosenblatt transform.
  #
  # inspired by gaussianize.r and.R, split.m in normal.m by Van Albada, S.J., Robinson P.A. (2006)
  # Transformation of arbitrary distributions to the normal distribution with application to EEG
  # test-retest reliability. J Neurosci Meth, doi:10.1016/j.jneumeth.2006.11.004
  #
  #  Modified by matlab code written 26/06/2015 by Julien Emile-Geay (USC)
  #translated to R and added jitter option 7/06/2017 by Nick McKay (NAU) 
  
  if(!is.matrix(X)){
    X=as.matrix(X)
  }
  p=NCOL(X)
  n=NROW(X) 
  
  if(jitter){
    #add tiny random numbers to avoid ties
    X=array(rnorm(p*n,mean=0,sd=sd(as.vector(X))/1e6),c(n,p))+X
  }
  
  
  Xn    = matrix(0,n,p);
  for (j in 1:p){
    # Sort the data in ascending order and retain permutation indices
    R=rank(X[,j])
    # The cumulative distribution function
    CDF = R/n - 1/(2*n);
    # Apply the inverse Rosenblatt transformation
    Xn[,j] = qgamma(CDF,shape = shape,rate = shape/mean)  # Xn is now gamma distributed
  }
  
  return(Xn)
}
