# bACEian

This repo documents some tests with bayesian twin analysis on simulated data. The objective is to see whether it is possible and efficient to estimate twin models in Stan (or in the future, other bayesian estimation programs).

The work is inspired by Espen Eilertsen's paper in Norsk Epidemiologi a few years back. Here's the URL:
https://www.ntnu.no/ojs/index.php/norepid/article/view/2017 He also put some stuff on the web, that I used.

Since then Ole Røgeberg (Frisch centre) helped tremendously by rewriting the model which increased the efficiency of the Stan runs.

What is in the repo now, are model definitions in Stan and a R script (rogeberg.R) that show a couple of experiments with twin data.

TODO:
- Include other types of relatives
- Better diagnostics?

