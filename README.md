This repositoty includes algorithms for estimating total factor productivity with or without R&D.

So far, there are two programs that reproduce Levinsohn and Petrin's (2003) LS estimator. To achieve more speed, I implemented a NLLS estimator. Only a few papers do that (__More to say on this later ...__).

The rough idea is that bootstrap is not necessary to achieve efficiency of estimated coefficients in samples of firms and households as we use e.g. in my laboratory @Sciences Po.

As far as I know, the speed of estimators has not been considered as a subject. My routine, however, goes 3 to 6 times faster than that of PETRIN, A., POI, B. and LEVINSOHN, J., 2004, Production function estimation in Stata using inputs to control for unobservables, _The Stata Journal_, 4, 113-123, or ROVIGATTI, G. and MOLLISI, V., 2018, Theory and practice of total-factor productivity estimation: the control function approach using Stata, _The Stata Journal_, 18, 618-662.

My objective in this research is to obtain fast and efficient estimators of production function coefficients when the sample size includes millions of individuals.

You can join me on this if you want.
.

