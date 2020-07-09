We recommend instead for easier reading of documentation:  
https://localepi.github.io/LEMMA/

Forecasts and scenarios for California counties:  
https://github.com/LocalEpi/LEMMA-Forecasts/tree/master/Forecasts  
https://github.com/LocalEpi/LEMMA-Forecasts/tree/master/Scenarios

# LEMMA
LEMMA (Local Epidemic Modeling for Management &amp; Action) is designed to provide regional (e.g. city or county-level) projections of the SARS-CoV-2 (COVID-19) epidemic under various scenarios. Daily projections with uncertainty bounds are made for hospitalizations, ICU use, ventilator use, active cases, and total cases. As detailed below, LEMMA allows for a range of user-specified parameterizations (including on the model structure) and is fit using case series data of COVID-19 hospitalizations.

## Contributors
LEMMA is a collaborative effort between experts in Medicine, Public Health, and Data Science, including but not limited to

- Joshua Schwab - UC Berkeley
- Laura B. Balzer - UMass Amherst
- Elvin Geng - Washington University
- James Peng - UC San Francisco
- Maya L. Petersen - UC Berkeley

We have moved our model fitting from R to Stan. Our Stan implementation is based on the "Santa Cruz County COVID-19 Model" (https://github.com/jpmattern/seir-covid19) by Jann Paul Mattern (UC Santa Cruz) and Mikala Caton (Santa Cruz County Health Services Agency). We are very grateful to Paul and Mikala for generously sharing their code and helping us.

## Installation
1. Install RStudio. (https://rstudio.com/products/rstudio/download/#download)
2. Create a folder to store your LEMMA inputs and outputs. For example, create a folder "MyFolder" within Documents.
3. For this step there are several choices, depending on your local machine.

    **Installing From Source (available for all platforms, but requires a C/C++ compiler)**
    
    If you do not have a toolchain including a C/C++ compiler, see here for installing one:
    
    https://support.rstudio.com/hc/en-us/articles/200486498
```{r}
install.packages("remotes")  #if you do not already have the remotes package
remotes::install_github("LocalEpi/LEMMA")
```
    **Installing From Binary (for MacOS using R 3.6)**
```{r}
install.packages("https://github.com/joshuaschwab/LEMMAstan/blob/master/LEMMA_MacOS_R-3-6-3.tgz?raw=true", repos = NULL, type = "binary")
```
    
    **Installing From Binary (for Windows using R 3.6)**
```{r}
install.packages("https://github.com/joshuaschwab/LEMMAstan/blob/master/LEMMA_Win_R-3-6-3.zip?raw=true", repos = NULL, type = "binary")
```
    
    **Installing From Binary (for Windows using R 4.0)**
```{r}
install.packages("https://github.com/joshuaschwab/LEMMAstan/blob/master/LEMMA_Win_R-4-0.zip?raw=true", repos = NULL, type = "binary")
```

4. Copy and paste these lines into the RStudio console:
    
```{r}
setwd("~/Documents/MyFolder")   # replace "~/Documents/MyFolder" with the path/folder you created
file.copy(system.file("extdata", "template.xlsx", package = "LEMMA", mustWork = TRUE), "example.xlsx", overwrite = TRUE)
```

LEMMA is in early development and is changing rapidly. Please restart RStudio and repeat step 3 once per day.

If you get an error something like 
"lazy-load database 'Library/R/3.6/library/LEMMA/R/LEMMA.rdb' is corrupt"
restart R and try step 3 again.


## Running LEMMA
1. Edit the Excel file ~/Documents/MyFolder/example.xlsx and save under a new name. For example, "MyCity.xlsx"
2. Run the following code.
```{r}
LEMMA::CredibilityIntervalFromExcel("MyCity.xlsx")
```

## Additional Documentation 
(work in progress)

https://localepi.github.io/LEMMA/articles/case_study.html

https://localepi.github.io/LEMMA/articles/faq.html

## Input
The provided spreadsheet provides a template and example of the inputs needed to run LEMMA. These are inputted in four tabs.

### Sheet 1: Parameters with Distributions
Briefly, LEMMA requires parameters related to the epidemic modeling (e.g., basic reproductive number, duration of infectiousness, percent of infected persons who are hospitalized). LEMMA also allows the user to specify the timing and impact of public health interventions, such as school closures and shelter-in-place orders. Interventions may occur before the current date to reflect such public health interventions. Interventions may also occur after the current date and can be used to simulate epidemic if measures are implemented or lifted at a future date. Explanations for specific parameters are provided below. Users can input a mean and standard deviation for each parameter. Each parameter will be drawn from a normal distribution. 

- Basic reproductive number R0 before Intervention1: This should represent initial epidemic growth before any public health interventions were implemented.

### Sheet 2: Interventions

- Intervention Date: Specify the date of the public health intervention. You can set standard deviation to zero to indicate that the intervention definitely starts on a certain day, or set standard deviation to a positive number to indicate there is some uncertainty about when the intervention took/takes effect. Intervention Date can be in the past or in the future. 

- Re Multiplier: Specify the impact of the first intervention in terms of multiplicative reductions in the basic reproductive number. Suppose this is 60%. Then the effective reproductive number after the first intervention would be $Re = 0.6 * R0$, where R0 is the basic reproductive number provided on the Parameters with Distributions sheet. If Re Multiplier for the second intervention is 50%, the effective reproductive number after the second intervention would be $Re = 0.5 * 0.6 * R0$.

- Days to reach new Re: LEMMA assumes the effects of interventions do not happen instantaneously. Therefore, specify the number of days to reach the new effective reproductive number. 

### Sheet 3: Model Inputs
Specify the number of people in the area of interest; this could be a region, city, or hospital catchement area. Specify the starting and final date of projections. 

### Sheet 4: Data
Provide hospital, ICU and/or death case series data. PUI (Persons Under Investigation, or "Probable" cases) can be entered if available. Any entries (either an entire column or specific rows) can be left blank if the data is not available. 
- Hospitalizations: Number of patients with COVID19 hospitalized on a given day, *including* those in ICU.
- ICU: Number of patients in ICU with COVID19 on a given day
- Cumulative Deaths: Total number of persons who have died due to COVID19 by a given day
- Cumulative Hospitalizations: Total number of patients who have been hospitalized with COVID19 by a given day

### Sheet 5: PUI Details
If PUIs are used on the Data sheet, a mean and standard deviation for the fraction of PUIs who are actually COVID19 positive can be entered. If PUIs are not used for a given category on the Data sheet, the mean and standard deviation for that category on the PUI Details sheet will be ignored.

### Sheet 6: Internal
This allows for more nuanced changes including changing the file names (output.filestr) and plotting details (plot.observed.data.long.term, plot.observed.data.short.term).
Control arguments passed to Stan can also be specified (iter, cores, refresh, max_treedepth, adapt_delta). See `?rstan::stan` for details. 
Note that iter is now the number of *posterior* draws (in previous versions of LEMMA, the number of prior iterations was specified.) The total number of posterior draws is 2 * iter because there are 4 chains, each with iter draws, but 50% of these are warmup draws which are not used.

## Output 
Given the above input, LEMMA runs specified number of simulations, where in each simulation the parameters are sampled from the specified prior distributions. The default total number of simulations is 2000, but can be changed with iter on the Internal tab (see Sheet 6: Internal, above). The posterior distribution is provided for the model outputs: number hospitalized, number in ICU, number of deaths, cumulative number of admissions.

The main output is provided in pdf format. Plots include short term and long term projections of umber hospitalized, in the ICU, cumulative deaths and cumulative hospital admissions (these are only shown for categories in which data was entered on the Data sheet). A histogram of Re and a plot of Re over time are shown up to 14 days before the last observed data. It is difficult to estimate Re beyond that date because it takes at least two weeks for changes in Re to be reflected in hospitalizations. Posterior distributions for each parameter are also shown. 

The Excel spreadsheet output provides detailed output for daily number hospitalized, in the ICU, cumulative deaths and cumulative hospital admissions.


### License
 
The MIT License (MIT)

Copyright (c) 2020 LEMMA

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
