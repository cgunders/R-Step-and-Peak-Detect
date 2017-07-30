library(abf2)
library(changepoint)
library(Peaks)
library.dynam('Peaks', 'Peaks', lib.loc=NULL)
library(ggplot2)
library(stringr)
require(gridExtra)

#data class that will be store x and y, fit parameters, and fit results

setClass(Class="TraceData",
         representation(
             x="vector",
             y="vector",
             par="vector",
             cpt="data.frame",
             pks="data.frame",
             bsl="vector",
             name="character",
             volt="numeric"
             )
)

#step and peak detection function

#x and y are vectors containing x and y data of time series trace
#k and msl are changepoint detection parameters, k is floating point, msl is int
#sig and thresh are peak detection parameters, both floating points
detect <- function(TraceData) {

    x = TraceData@x
    y = TraceData@y
    k = TraceData@par[1]
    msl = TraceData@par[2]
    sig = TraceData@par[3]
    thresh = TraceData@par[4]
    
    
    #run changepoint detection using PELT with user-specificed 
    #penalty value k and minimum changepoint length msl
    step_locs = cpts(cpt.mean(y, penalty='Manual', pen.value=k*log(length(y)),
                              method='PELT', minseglen=msl))
    
    #initialize vectors
    step_heights = vector()
    step_means = vector()
    
    #number of steps
    n_steps = length(step_locs)
    
    #calculate initial elements of step height and changepoint location vectors
    step_heights[1] = mean(y[step_locs[1]:step_locs[2]]) - mean(y[1:step_locs[1]])
    
    step_means[1] = mean(y[1:step_locs[1]])
    
    #calculate final elements of step height and changepoint location vectors
    step_heights[n_steps] = mean(y[tail(step_locs,1):length(y)]) - 
                    mean(y[step_locs[n_steps-1]:tail(step_locs,1)])
    
    step_means[n_steps] = mean(y[step_locs[n_steps-1]:tail(step_locs,1)])

    #calculate middle elements of step height and changepoint location vectors
    for(j in 2:(n_steps-1)){
        
        step_heights[j] = mean(y[step_locs[j]:step_locs[j+1]]) - mean(y[step_locs[j-1]:step_locs[j]])
        step_means[j] = mean(y[step_locs[j-1]:step_locs[j]])

    }

    #store changepoint data into data frame
    steps = data.frame(step_locs,step_heights,step_means)
    names(steps) = c('Step locations', 'Step heights', 'Step means')
    
    #initialize step baseline vector
    step_baseline = vector()

    #calculate initial and final elements of step baseline vector
    step_baseline[1:step_locs[1]] = step_means[1]
    
    step_baseline[step_locs[n_steps]:length(x)] = 
                    step_means[n_steps] + step_heights[n_steps]

    #calculate middle elements of step baseline vector
    for(i in 2:n_steps){
    
        step_baseline[step_locs[i-1]:step_locs[i]] = step_means[i]

    }
    
    #subtract steps from time series
    new_y = y - step_baseline
    
    #search for peak locations in subtracted time series with user-specified 
    #sigma sig and threshold thresh
    info = SpectrumSearch(new_y,sigma=sig,threshold=thresh,background=FALSE,
                          iterations=13,markov=FALSE,window=3)
    peak_locs = sort(info$pos, decreasing = FALSE, na.last = NA)

    #put peak locations back into subtracted baseline to get peak heights
    peak_heights = c(new_y[peak_locs])
    
    #store peak location and height into data frame
    peaks = data.frame(peak_locs,peak_heights)
    names(peaks) = c('Peak locations', 'Peak heights')
    
    #export
    TraceData@cpt <- steps
    TraceData@pks <- peaks
    TraceData@bsl <- step_baseline
    
    return(TraceData)
    
}

#plot for initial inspection

plot_results <- function(TraceData) {
    
    theme_set(theme_classic())
    
    lay <- rbind(c(1,2), c(3,4))

    plotdata <- data.frame(TraceData@x, TraceData@y)
                       #TraceData[[i]]@cpt[,1], TraceData[[i]]@cpt[,2],
                       #TraceData[[i]]@cpt[,3], TraceData[[i]]@pks[,1],
                       #TraceData[[i]]@pks[,2], TraceData[[i]]@bsl)

    plot1 <- ggplot(plotdata, aes(x=TraceData@x)) +
        geom_line(aes(y=TraceData@y), color="black") +
        labs(x=("t / s"), y=("i / pA"), 
             title= str_c(TraceData@name, 
                          as.character(TraceData@volt), 
                          sep=", ")) +
        theme(text = element_text(size=20))

    plot2 <- ggplot(plotdata, aes(x=TraceData@x)) +
        geom_line(aes(y=TraceData@y), color="black") +
        geom_line(aes(y=TraceData@bsl), color="red", lwd=0.75) +
        labs(x=("t / s"), y=("i / pA"), 
             title= str_c(TraceData@name, 
                          as.character(TraceData@volt), 
                          sep=", ")) +
        theme(text = element_text(size=20))

    plot3 <- ggplot(plotdata, aes(x=TraceData@x)) +
        geom_line(aes(y=TraceData@y), color="black") +
        labs(x=("t / s"), y=("i / pA"), 
             title= str_c(TraceData@name, 
                          as.character(TraceData@volt), 
                          sep=", ")) +
        annotate("text",  x=TraceData@pks[,1]/1000, 
            y=c(TraceData@y[TraceData@pks[,1]]), 
            label=as.character(round(TraceData@pks[,2],1)),
            color="blue", size = 3, hjust = 1, vjust = 0) +
        theme(text = element_text(size=20))

    plot4 <- ggplot(plotdata, aes(x=TraceData@x)) +
        geom_line(aes(y=TraceData@y), color="black") +
        geom_line(aes(y=TraceData@bsl), color="red", lwd=0.75) +
        labs(x=("t / s"), y=("i / pA"), 
             title= str_c(TraceData@name, 
                          as.character(TraceData@volt), 
                          sep=", ")) +
        annotate("text",  x=TraceData@cpt[,1]/1000, 
            y=TraceData@cpt[,3], 
            label=as.character(round(TraceData@cpt[,2],1)),
            color="red", size = 3, hjust = -0.25, vjust = 1.5) +
        theme(text = element_text(size=20))

    grid.arrange(plot1,plot2,plot3,plot4,layout_matrix=lay)
}

#data import from .abf files

ablist = dir('.../AnalysisFolder',
             full.names=TRUE)

#store data into data objects

tracelist = list()

for(i in 1:length(ablist)){

    ab = abfload(filename=ablist[i])
    
    tracelist[[i]] <- new("TraceData")
    
    #store x and y from trace
    tracelist[[i]]@x <- ab$s
    tracelist[[i]]@y <- as.vector(ab$traces)
    
    #store fitting parameters
    tracelist[[i]]@par <- c(150, 200, 5, 15)
    
    #store name of pore used in experiment
    tracelist[[i]]@name <- str_extract(ablist[i], "([:alnum:])+(?=_)")
    
    #store voltage used for experiment
    tracelist[[i]]@volt <- as.integer(str_extract(str_extract(
        ablist[i],"(?<=_).\\d.+(?=_)"), "\\d+"))/(-1000)
    
    #perform and store step and peak detection
    tracelist[[i]] = detect(tracelist[[i]])
    
    #optional plotting
    #plot_results(tracelist[[i]])
    
}

print("Finished")
