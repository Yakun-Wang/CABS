# CABS
Conditional Adaptive Bayesian Spectrum Analysis (CABS)
"Adaptive Bayesian Spectral Analysis of Nonstationary Biomedical Time Series"
by Bruce, Hall, Buysse, and Krafty (2017)

Author: Yakun Wang

Description: Instructions for implementing CABS estimation from "Adaptive Bayesian
Spectral Analysis of Nonstationary Biomedical Time Series"
by Bruce, Hall, Buysse, and Krafty (2017)

Dependencies:
Code was developed using Rstudio (Version 1.1.463), so code may not function properly
on other versions of Rstudio. Several R packages need to be installed which is mentioned in the code.


Quick start guide:
Follow the steps below to simulate data, run the CABS estimation procedure,
and create data visualizations of the method to simulated slowly varying and 
piecewise AR processes detailed in the paper.

Download the demo.R as well as CABS.R file and extract the contents to a folder of your
choosing. In what follows, we will assume the folder is saved with the
following path 'C:\CABSdemo'.

Open the file 'C:\CABSdemo\demo.R'. This script contains
all necessary code to reproduce the data and MCMC sampler results for the
piecewise and slowly varying AR processes.

Follow the instructions in the comments of the demo file to simulate
data from the piecewise and slowly varying AR processes described in the
paper and apply the CABS procedure to obtain spectral estimates and summary
plots of the MCMC sampler fit. 

When you want to call the cabs function in "CABS.R" file, you need to source the code
with command 'source('C:\CABSdemo\CABS.R')'.

Using CABS on other data:
Two inputs are necessary to run the MCMC sampler for the CABS estimation
procedure.

'x' is a matrix whereby each column contains a realization of the time
series for a given subject. The rows are indexed by time
(t=1,...,T) and the columns are indexed by subject (l=1,...,L). For
example, {x}_11 contains the first time series instance for the first
subject.

'u' is a vector containing the corresponding covariate values for each
subject. The columns are indexed by subject, which is the same as for
the time series matrix 'x'. For example, {u}_1 is the covariate value
for the first subject.

And some parameters mentioned in the 'CABS.R' need to be set up for apply the
CABS function to the simulation data x and u. 

Once you have created these inputs based on your data, you can pass
them into the 'CABS' function just as in the demo file for estimation.

MIT License

Copyright (c) 2019 Yakun Wang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
