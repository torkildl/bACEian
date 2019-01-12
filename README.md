# bACEian

This repo documents my tests with bayesian twin analysis on simulated data. The objective is to see whether it is possible and efficient to estimate twin models in Stan (or in the future, other bayesian estimation programs). 

Look at the script "fit.R" for the current code, and "model.stan" for the Stan model specification for more info.

The work is inspired by Espen Eilertsen's paper in Norsk Epidemiologi a few years back. Here's the URL:
https://www.ntnu.no/ojs/index.php/norepid/article/view/2017 He also put some stuff on the web, including the backbone of the script I am using here ("model.stan")


TODO:
- Include other types of relatives
- Try other efficiency tweaks
