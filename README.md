# Levinsohn and Petrin's (2003) estimator (__LP__)

This routine obtains fast and efficient estimators of production function coefficients in large samples. I provide the different stages of the __LP__ estimator with more details than I could find elsewhere. Then, I introduce my estimator. 

## The different stages in __LP__

![Stages.](http://www.evens-salies.com/tfp_lp_nlls.png)

## The routine ```tfp_lp_nlls.do```

This repository includes two programs that reproduce Levinsohn and Petrin's (2003) LS estimator in __Stata__ but with a higher speed than in ```levpet``` and ```profest```. To achieve more speed, I implemented a nonlinear least squares (__NLLS__) estimator. The NLLS estimator estimates the capital elasticity and the coefficients of the markov equation simultaneously. Basicaly, step 4-6 are merged into one __Stata__ command. The routine does not rely on the bootstrap to achieve efficiency of estimated coefficients. Bootstrapping is more relevant in small samples. The present routines should be used in large samples of firms and households (millions) as we use e.g. in my laboratory @Sciences Po.

## Result

As far as I know, the speed of productivity estimation is not very considered as a subject per se. My routine, however, goes 3 to 6 times faster than that of Petrin, A., Poi, B. and Levinsohn, J., 2004, Production function estimation in Stata using inputs to control for unobservables, _The Stata Journal_, 4, 113-123, or Rovigatti, G. and Mollisi, V., 2018, Theory and practice of total-factor productivity estimation: the control function approach using Stata, _The Stata Journal_, 18, 618-662.

## The variables that you must have in your data sheet

 I decided not to make a command, for the code is short so it is more interesting I think for you to copy-paste it into your own code. You'll just have to rename your variables (firm id, time, added-value, etc.) to:
 
 - firm id ```siren```,
 - time ```year```,
 - log of added value ```l_v```,
 - capital stock ```l_k```,
 - log of total hours of work ```l_hours```,
 - log of materials ```l_m```.

__Bon voyage__
