# R-Step-and-Peak-Detect

Published at https://doi.org/10.1021/acs.langmuir.7b04090

This is a script I created to measure steps and peaks in electrochemical current-time traces I collected during experiments. It uses the Changepoint package to find the changepoints in the trace. The changepoints are subtracted from the trace to generate a flat baseline with peaks that are detected using the Peaks package. 

The Multi File Input script MultiDetect.R is more readable. It sees how many traces I have stored as .abf files in the specified folder and imports them. There is an optional plotting function used during fitting.

Detection thresholds for both step and peak detection are important to avoid over-fitting. The minimum segment length in the changepoint detection can be used to avoid fitting peaks, leaving them for the peak detection algorithm afterward. I evaluate the goodness of the fit by inspection afterward.

![alt text](https://github.com/cgunders/R-Step-and-Peak-Detect/blob/master/InitialCurrentTrace.png "Initial Current Trace")

The current-time trace as initially imported.

![alt text](https://github.com/cgunders/R-Step-and-Peak-Detect/blob/master/DetectedTrace.png "Detected Current Trace")

The current-time trace with changepoint detection overlayed in red and detected peak heights in blue.

![alt text](https://github.com/cgunders/R-Step-and-Peak-Detect/blob/master/DetectedTraceZoomIn.png "Detected Current Trace Zoom-In")

The same, but zoomed in on a segment and with changepoint step heights in red.

![alt text](https://github.com/cgunders/R-Step-and-Peak-Detect/blob/master/SubtractedStepsPeakDetect.png "Subtracted Current Trace For Peak Detection")

The trace after steps are subtracted from the original trace with peaks remaining. This is the trace that step detection is performed on.

Credit:

Uses abf2 library: https://github.com/cran/abf2

Uses Changepoint library: https://github.com/rkillick/changepoint/

Uses Peaks library: https://cran.r-project.org/web/packages/Peaks/index.html
