# R-Step-and-Peak-Detect

This is a script I created to measure steps and peaks in electrochemical current-time traces I collected during experiments. It uses the Changepoint package to find the changepoints in the trace. The changepoints are subtracted from the trace to generate a flat baseline with peaks that are detected using the Peaks package. 

The Multi File Input script MultiDetect.R is more readable. It sees how many traces I have stored as .abf files in the specified folder and imports them. There is an optional plotting function used during fitting.

Detection thresholds for both step and peak detection are important to avoid over-fitting. The minimum segment length in the changepoint detection can be used to avoid fitting peaks, leaving them for the peak detection algorithm afterward. I evaluate the goodness of the fit by inspection afterward.

InitialCurrentTrace.png shows the current-time trace as initially imported.

DetectedTrace.png shows the current-time trace with changepoint detection overlayed in red and detected peak heights in blue.

DetectedTraceZoomIn.png shows the same, but zoomed in on a segment and with changepoint step heights in red.

SubtracedStepsPeakDetect.png shows the trace after steps are subtracted from the original trace with peaks remaining. This is the trace that step detection is performed on.

Figure.png shows a figure that I made with the useful information (just step heights currently) obtained using this algorithm on my experimental data.

Credit:

Uses abf2 library: https://github.com/cran/abf2

Uses Changepoint library: https://github.com/rkillick/changepoint/

Uses Peaks library: https://cran.r-project.org/web/packages/Peaks/index.html
