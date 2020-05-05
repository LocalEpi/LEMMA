# LEMMA
LEMMA (Local Epidemic Modeling for Management &amp; Action) is designed to provide regional (e.g. city or county-level) projections of the SARS-CoV-2 (COVID-19) epidemic under various scenarios. Daily projections with uncertainty bounds are made for hospitalizations, ICU use, ventilator use, active cases, and total cases. As detailed below, LEMMA allows for a range of user-specified parameterizations (including on the model structure) and is fit using case series data of COVID-19 hospitalizations.

# Contributors:

LEMMA is a collaborative effort between experts in Medicine, Public Health, and Data Science, including but not limited to

- Maya L. Petersen - UC Berkeley
- Joshua Schwab - UC Berkeley
- Laura B. Balzer - UMass Amherst
- Elvin Geng - Washington University
- James Peng - UC San Francisco

# Installation: 

1) Install RStudio. (https://rstudio.com/products/rstudio/download/#download)
2) Create a folder to store your LEMMA inputs and outputs. For example, create a folder "MyFolder" within Documents.
3) Copy and paste these lines into the RStudio console, one at a time:
```{r}
install.packages("remotes")  
remotes::install_github("LocalEpi/LEMMA")
setwd("~/Documents/MyFolder")   # replace "~/Documents/MyFolder" with the path/folder you created, or use Session > Set Working Directory > Choose Directory in RStudio
file.copy(system.file("extdata", "SF-April13.xlsx", package = "LEMMA", mustWork = TRUE), "example.xlsx")
```

** LEMMA is in early development and is changing rapidly. Please restart RStudio and repeat step 3 once per day. **

If you get an error something like 
"lazy-load database 'Library/R/3.6/library/LEMMA/R/LEMMA.rdb' is corrupt"
restart R and try step 3 again.

# Running LEMMA
1) Edit the Excel file ~/Documents/MyFolder/example.xlsx and save under a new name. For example, "MyCity.xlsx"
2) Run the following code.
```{r}
LEMMA::CredibilityIntervalFromExcel("MyCity.xlsx")
```
# Case Studies (package vignette)
To build the vignette:
```{r}
remotes::install_github("LocalEpi/LEMMA", build_vignettes = T)
browseVignettes(package="LEMMA")
```


# Input:

The provided spreadsheet provides a template and example of the inputs needed to run LEMMA. These are inputted in four tabs.

### Tab1: Parameters with Distributions
Briefly, LEMMA requires parameters related to the epidemic modeling (e.g., basic reproductive number, duration of infectiousness, percent of infected persons who are hospitalized). LEMMA also allows the user to specify the timing and impact of public health interventions, such as school closures and shelter-in-place orders. The current implementation allows for 3 interventions; the first two are assumed to occur before the current date, and the third can be used to simulate epidemic if measures are implemented or lifted at a future date. Explanations for specific parameters are provided below.
- Priors (Row 3): Specify the sampling probabilities for each parameter. You can change these as long as they sum to 100% The E-column corresponds to your best guess and should be given the most weight.
- Basic reproductive number R0 before Intervention1: This should represent initial epidemic growth before any public health interventions were implemented.
- Number of Days from Infection to Becoming Infectious (Latent Period): This must be an whole number. Specifying 0 indicates there is no latent period and an SIR (versus SEIR) epidemic model will be run. 
- Date of first intervention: Specify the date of the first public health intervention. 
- Re multiplier: Specify the impact of the first intervention in terms of multiplicative reductions in the basic reproductive number. Supposed this is 60%. Then the effective reproductive number after the first intervention would be Re= 0.6*R0, where R0 is the basic reproductive number provided on Row 4. 
- Days to reach new Re: LEMMA assumes the effects of interventions do not happen instantaneously. Therefore, specify the number of days to reach the new effective reproductive number. 
- Patients in hospital are infectious: This is a true/false indicator if patients are infectious in the hospital.
- Contant rate to hospital (if FALSE, fixed number of days to hospital): This is a true/false indicator if there is a constant rate at which infectious persons are hospitalized. (Note this is different than the "Percent of Infected that are Hospitalized" (Row 8).)

### Tab2: Model Inputs
Specify the number of people in the area of interest; this could be a region, city, or hospital catchement area. Specify the final date of projections. 

### Tab3: Hospitalization Data
Provide hospital case series data. On each date, specify the a lower and upper bound on the number of persons hospitalized with COVID-19. This bounding is intended to account for uncertainty due to persons who are hospitalized but under investigation (i.e., not yet confirmed COVID-19 positive or not). You may enter any range of dates. 

### Tab4: Internal
This allows for more nuanced changes including changing the file names (output.filestr), the number of iterations (main.iterations), and tolerance of projections to the observed hospital case seres (required.in.bounds). 

# Output: 
Given the above input, LEMMA runs specified number of simulations, where in each simulation the parameters are sampled from the specified prior distributions. (The default number of simulations is 10,000, but can be changed with main.iterations on the Internal tab.) The projections that are sufficiently close to the observed hospitalization case series are saved and provide a posterior distribution for the model outputs: number hospitalized, number in ICU, number on ventilators, active cases, and total cases. We recommend finding the parameterizations resulting in at least 1,000 posterior projections and urge  **CAUTION** when interpretting results with few posterior projections. 

The main output is provided in pdf format. 
- Page 1 is a graph of projected hospitalizations by date with uncertainity bounds. The yellow line is the projection under your best guess of parameters (represented by column E of the 'Parameters with Distributions' Tab). The red line is the median of the posterior projections. Uncertainty bands are represented by the shaded regions, representing quantiles of the posterior projections. 25%-75% of the posterior distribution falls within the dark gray band, 15%-85% within the medium gray band, and 5%-95% within the light gray bands. The black dots are the observed hospitalizations (currently assumed to be the average of the lower and upper bounds). 
- Page 2 provides the prior distribution for the models considered. These models are derived by a combination of 3 parameters where 0 denotes false and 1 true. 'hasE' indicates if there is a latent period; 'infect in hosp' indicates if patients are infectious in the hospital, and 'rate to hospital' indicates if there is a constant rate at which infectious persons are hospitalized. 
- Page 3 provides the posterior distribution for the models considered. Importantly, at the top of this graph and all 'posterior' pages, you will see the number of simulations included in the posterior distribution (niter=#). 
- Page 4 provides the prior distribution on the "current" effective reproductive number (currentRe). This is the reproductive number after the first two inteventions have taken place (but not the third which is assumed to be at a future date). 
- Page 5 provides the posterior distribution for the current effective reproductive number (currentRe).
- The following pages provide the prior and posterior distributions for all the specified parameters. 

The resulting spreadsheet provides detailed output for daily number hospitalized, in the ICU, on ventilators, as well as active cases and resolved cases. For each, the best guess and quantiles are provided.

# License:
 
The MIT License (MIT)

Copyright (c) 2020 LEMMA

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
