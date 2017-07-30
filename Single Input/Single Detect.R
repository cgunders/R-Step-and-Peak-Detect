#necessary libraries
library(changepoint)
library(Peaks)
library.dynam('Peaks', 'Peaks', lib.loc=NULL)

#data class that will be passed back from the function
setClass(Class="pkstp",
         representation(
           cpt="data.frame",
           pks="data.frame",
           bsl="vector"
         )
)

#x and y are vectors containing x and y data of time series trace
#k and msl are changepoint detection parameters, k is floating point, msl is int
#sig and thresh are peak detection parameters, both floating points
detect <- function(x,y,k,msl,sig,thresh) {
  
  #run changepoint detection using PELT with user-specificed 
  #penalty value k and minimum changepoint length msl
  z = cpt.mean(y, penalty='Manual', pen.value=k*log(length(y)),
               method='PELT', minseglen=msl)
  
  #initialize vectors
  a = vector()
  c = vector()
  b = cpts(z)
  
  #calculate initial elements of step height and changepoint location vectors
  a[1] = mean(y[b[1]:b[2]]) - mean(y[1:b[1]])
  c[1] = mean(y[1:b[1]])
  
  #calculate final elements of step height and changepoint location vectors
  a[length(b)] = mean(y[b[length(b)]:length(y)]) - 
    mean(y[b[length(b)-1]:b[length(b)]])
  c[length(b)] = mean(y[b[length(b)-1]:b[length(b)]])
  
  #calculate middle elements of step height and changepoint location vectors
  m = length(b) - 1
  
  for(j in 2:m){
    
    a[j] = mean(y[b[j]:b[j+1]]) - mean(y[b[j-1]:b[j]])
    c[j] = mean(y[b[j-1]:b[j]])
    
  }
  
  #store changepoint data into data frame
  g = data.frame(a,c,b)
  names(g) = c('Step height','Mean','Cpt loc')
  
  #initialize detected baseline vector
  lower_baseline = vector()
  
  n = length(g[,2])
  
  #calculate initial and final elements of new baseline vector
  lower_baseline[1:g[1,3]] = g[1,2]
  lower_baseline[g[length(g[,3]),3]:length(x)] = 
    g[length(g[,2]),2] + g[length(g[,1]),1]
  
  #calculate middle elements of new baseline vector
  for(i in 2:n){
    
    lower_baseline[g[i-1,3]:g[i,3]] = g[i,2]
    
  }
  
  #subtract steps from old time series
  new_y = y - lower_baseline
  
  #search for peak locations in subtracted time series with user-specified 
  #sigma sig and threshold thresh
  info = SpectrumSearch(new_y,sigma=sig,threshold=thresh,background=FALSE,
                        iterations=13,markov=FALSE,window=3)
  pks = sort(info$pos, decreasing = FALSE, na.last = NA)
  
  #put peak locations back into subtracted baseline to get height
  pkh = c(new_y[pks])
  
  #store peak location and height into data frame
  h = data.frame(pks,pkh)
  names(h) = c('Peak loc','Peak height')
  
  #export
  data = new("pkstp")
  data@cpt <- g
  data@pks <- h
  data@bsl <- lower_baseline
  
  return(data)
}
