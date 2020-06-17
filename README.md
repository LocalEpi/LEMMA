# LEMMA
LEMMA (Local Epidemic Modeling for Management &amp; Action) is designed to provide regional (e.g. city or county-level) projections of the SARS-CoV-2 (COVID-19) epidemic under various scenarios. Daily projections with uncertainty bounds are made for hospitalizations, ICU use, ventilator use, active cases, and total cases. As detailed below, LEMMA allows for a range of user-specified parameterizations (including on the model structure) and is fit using case series data of COVID-19 hospitalizations.

## Contributors
LEMMA is a collaborative effort between experts in Medicine, Public Health, and Data Science, including but not limited to

- Joshua Schwab - UC Berkeley
- Laura B. Balzer - UMass Amherst
- Elvin Geng - Washington University
- James Peng - UC San Francisco
- Maya L. Petersen - UC Berkeley

We are in the process of moving our model fitting from R to Stan. Our Stan implementation is based on the "Santa Cruz County COVID-19 Model" (https://github.com/jpmattern/seir-covid19) by Jann Paul Mattern (UC Santa Cruz) and Mikala Caton (Santa Cruz County Health Services Agency). We are very grateful to Paul and Mikala for generously sharing their code and helping us.

## Installation
1. Install RStudio. (https://rstudio.com/products/rstudio/download/#download)
2. Create a folder to store your LEMMA inputs and outputs. For example, create a folder "MyFolder" within Documents.
3. Copy and paste these lines into the RStudio console, one at a time:
```{r}
install.packages("remotes")  
remotes::install_github("LocalEpi/LEMMA")
setwd("~/Documents/MyFolder")   # replace "~/Documents/MyFolder" with the path/folder you created
file.copy(system.file("extdata", "SF-April13.xlsx", package = "LEMMA", mustWork = TRUE), "example.xlsx")
```

If you would like to try our significantly improved version (currently in beta):
```{r}
install.packages("remotes")  
remotes::install_github("LocalEpi/LEMMA@Stan")
setwd("~/Documents/MyFolder")   # replace "~/Documents/MyFolder" with the path/folder you created
file.copy(system.file("extdata", "template.xlsx", package = "LEMMA", mustWork = TRUE), "example.xlsx")
```

** LEMMA is in early development and is changing rapidly. Please restart RStudio and repeat step 3 once per day. **

If you get an error something like 
"lazy-load database 'Library/R/3.6/library/LEMMA/R/LEMMA.rdb' is corrupt"
restart R and try step 3 again.


## Running LEMMA
1. Edit the Excel file ~/Documents/MyFolder/example.xlsx and save under a new name. For example, "MyCity.xlsx"
2. Run the following code.
```{r}
LEMMA::CredibilityIntervalFromExcel("MyCity.xlsx")
```

## Documentation 
(work in progress)

https://localepi.github.io/LEMMA/articles/case-studies.html

https://localepi.github.io/LEMMA/articles/faq.html


## Input
The provided spreadsheet provides a template and example of the inputs needed to run LEMMA. These are inputted in four tabs.

### Sheet 1: Parameters with Distributions
Briefly, LEMMA requires parameters related to the epidemic modeling (e.g., basic reproductive number, duration of infectiousness, percent of infected persons who are hospitalized). LEMMA also allows the user to specify the timing and impact of public health interventions, such as school closures and shelter-in-place orders. The current implementation allows for 3 interventions; the first two are assumed to occur before the current date, and the third can be used to simulate epidemic if measures are implemented or lifted at a future date. Explanations for specific parameters are provided below.

- Priors (Row 3): Specify the sampling probabilities for each parameter. You can change these as long as they sum to 100% The E-column corresponds to your most likely prior and should be given the most weight.

- Basic reproductive number R0 before Intervention1: This should represent initial epidemic growth before any public health interventions were implemented.

- Number of Days from Infection to Becoming Infectious (Latent Period): This must be an whole number. Specifying 0 indicates there is no latent period and an SIR (versus SEIR) epidemic model will be run.

- Date of first intervention: Specify the date of the first public health intervention. 

- Re multiplier: Specify the impact of the first intervention in terms of multiplicative reductions in the basic reproductive number. Suppose this is 60%. Then the effective reproductive number after the first intervention would be Re= 0.6*R0, where R0 is the basic reproductive number provided on Row 4. 

- Days to reach new Re: LEMMA assumes the effects of interventions do not happen instantaneously.
Therefore, specify the number of days to reach the new effective reproductive number. 

- Patients in hospital are infectious: This is a true/false indicator if patients are infectious in the hospital.

- Contant rate to hospital (if FALSE, fixed number of days to hospital): This is a true/false indicator if there is a constant rate at which infectious persons are hospitalized. (Note this is different than the "Percent of Infected that are Hospitalized" (Row 8).)

### Sheet 2: Model Inputs
Specify the number of people in the area of interest; this could be a region, city, or hospital catchement area. Specify the starting and final date of projections. 

### Sheet 3: Hospitalization Data
Provide hospital case series data. On each date, specify the a lower and upper bound on the number of persons hospitalized with COVID-19. This bounding is intended to account for uncertainty due to persons who are hospitalized but under investigation (i.e., not yet confirmed COVID-19 positive or not). You may enter any range of dates. 

### Sheet 4: Internal
This allows for more nuanced changes including changing the file names (output.filestr), the number of iterations (main.iterations), and tolerance of projections to the observed hospital case seres (required.in.bounds). 


## Output 
Given the above input, LEMMA runs specified number of simulations, where in each simulation the parameters are sampled from the specified prior distributions. (The default number of simulations is 10,000, but can be changed with main.iterations on the Internal tab.) The projections that are sufficiently close to the observed hospitalization case series are saved and provide a posterior distribution for the model outputs: number hospitalized, number in ICU, number on ventilators, active cases, and total cases. We recommend finding the parameterizations resulting in at least 1,000 posterior projections and urge  **CAUTION** when interpretting results with few posterior projections. 

The main output is provided in pdf format. 

- Page 1 is a graph of short term hospitalization projections by date with uncertainity bounds. The dark blue line is the median of the posterior projections. Uncertainty bands are represented by the shaded regions, representing quantiles of the posterior projections. 25%-75% of the posterior distribution falls within the dark purple band, 15%-85% within the medium purple band, and 5%-95% within the light purple bands. Green and red Xs denote the lower bound and upper bounds of the hospitalization case series data. (You can change these labels on the Internal Tab.)

- Page 2 provides a graph of long term hospitalization projections by date with uncertainity bounds.

- The remaining pages provide the prior and posterior distributions for all the specified parameters. Importantly, at the top of all 'posterior' pages, you will see the number of simulations included in the posterior distribution (niter=#). Pages 31 and 32 provide the prior and posterior distributions for the models considered. These models are derived by a combination of 3 parameters where 0 denotes false and 1 true. 'hasE' indicates if there is a latent period; 'infect in hosp' indicates if patients are infectious in the hospital, and 'rate to hospital' indicates if there is a constant rate at which infectious persons are hospitalized. 

The Excel spreadsheet output provides detailed output for daily number hospitalized, in the ICU, on ventilators, as well as active cases and resolved cases. For each, the User's Prior Projection and quantiles are provided.


### License
 
The MIT License (MIT)

Copyright (c) 2020 LEMMA

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
