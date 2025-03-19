# Levinsohn and Petrin's (2003) estimator (__LP__)

This routine obtains fast and efficient estimators of production function coefficients. The routine simplifies LP a bit so it should be useful in large samples (I'll post a benchmark comparison soon). 

The section below documents the different stages of the __LP__ estimator (in French so far ...). Then, I suggest a <ins> Nonlinear Least Squares approach </ins> to LP. Mail me with _warm cheers_ at [evens.salies@sciencespo.fr](mailto:evens.salies@sciencespo.fr) if my presentation below is useful to your work.

## Estimation stages in __LP__

<img src="http://www.evens-salies.com/tfp_lp_nlls.png" alt="Estimation stages" width="200">

## The routine ```tfp_lp_nlls.do```

This repository includes a program that implement Levinsohn and Petrin's (2003) LS estimator in __Stata__ but with a higher execution speed than in ```levpet``` and ```profest```. To achieve higher speed, I implement a nonlinear least squares (__NLLS__) estimator. The NLLS estimator estimates the capital elasticity and the coefficients of the markov equation simultaneously. Basicaly, steps 4-6 are merged into one __Stata__ command. The routine does not rely on the _bootstrap_ to achieve efficiency of estimated coefficients. Bootstrapping is useful for small samples, whereas this program should be used in large samples of firms and households (millions) as we use e.g. in my research unit [@SciencesPo](https://www.ofce.sciences-po.fr/en/).

As far as I know, the speed of productivity estimation is not very considered as a subject per se. However, when the sample charge is large and unless your computer reaches the speed of a rocket, my routine, goes 6 times faster than the pioneer ```levpet``` command by Petrin, A., Poi, B. and Levinsohn, J., 2004, Production function estimation in Stata using inputs to control for unobservables, _The Stata Journal_, 4, 113-123, or Rovigatti, G. and 3 times faster than ```prodest``` by Mollisi, V., 2018, Theory and practice of total-factor productivity estimation: the control function approach using Stata, _The Stata Journal_, 18, 618-662.

Another advantage is that you don't need to use your own grid search or other optimization algorithm, but instead benefit from the built in algorithms of Stata. Basically, by setting the options in the ```nl``` command. 

## The variables that you must have in your data sheet

I decided not to make a command, for the code is short enough; so, it is more interesting I think for you to copy-paste ```tfp_lp_nlls.do``` into your own code. Before, you'll have to rename your variables for firm id, time, added-value, etc. to:

 - firm id: ```siren```,
 - time: ```year```,
 - log of added value: ```l_v```,
 - capital stock: ```l_k```,
 - log of total hours of work: ```l_hours```,
 - log of materials: ```l_m```.

__Bon voyage__ 
:+1:
