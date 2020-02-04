
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BSSRed (Blinded Sample Size Reestimation with Evnt-Driven Data.)

<!-- badges: start -->

<!-- badges: end -->

The goal of BSSRed is to provide the methods published by T. Friede et
al. (Friede and Schmidli 2010), which main aspect lies on reestimation
procedures

## Installation

You can install the released version of BSSRed from
[GITHUB](https://github.com/ChristophAnten/BSSRed) with:

``` r
library("devtools")
install_git(url="https://github.com/ChristophAnten/BSSRed.git",
            upgrade=FALSE)
library("BSSRed")
```

## Study Exapmle

This is a basic example to evaluate and varify a study plan with event
driven data.

For this example it is assumed, that we have some prior information from
other studies on witch we base our assumptions on the hazard and drop
out rates.

We assume, that in the control group 30% of our patients have an event
after 24 month and 20% of the patients across both group drop out of the
study within 24 month. Assuming a 30% reduction of hazard in the control
group.

The recruitment is assumed to start a bit slowly in the beginning with
11 in the control and 22 in the treatment group, rising to 79 in the
control and 159 in the treatment group.

``` r
library(BSSRed)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!

<div id="refs" class="references">

<div id="ref-friede2010blinded">

Friede, Tim, and Heinz Schmidli. 2010. “Blinded Sample Size Reestimation
with Count Data: Methods and Applications in Multiple Sclerosis.”
*Statistics in Medicine* 29 (10): 1145–56.

</div>

</div>
