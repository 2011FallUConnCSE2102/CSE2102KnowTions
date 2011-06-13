ArmSampler implements an adaptive rejection Metropolis sampler (ARMS) that
can sample from virtually any univariate distribution. The method performs
best if a log-concave density is being sampled from, but it also works for
other densities, for which an additional Metropolis step is inserted. The
greater the difference to log-concave shape, the more Metropolis rejections
must be expected.

The central class is org.knowceans.arms.ArmSampler. When sampling a distri-
bution, override the abstract method logpdf(double x, Object params), as 
shown in the examples. The souce tree src-dep represents dependencies from
the knowceans-tools package for convenience but without guarantee that they
are the latest available versions.

This implementation is a port of the original C / Fortran implementation by
Wally Gilks available at
http://www.mrc-bsu.cam.ac.uk/BSUsite/Research/ars.shtml.

Please acknowledge this work by referencing the relevant scientific
literature and program code (Web: http://www.arbylon.net/projects).

References:

Gilks, W. R. (1992) Derivative-free adaptive rejection sampling for Gibbs
sampling. Bayesian Statistics 4, (eds. Bernardo, J., Berger, J., Dawid, A.
P., and Smith, A. F. M.) Oxford University Press.

Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995) Adaptive rejection
Metropolis sampling. Applied Statistics, 44, 455-472.

Gilks, W. R. and Wild, P. (1992) Adaptive rejection sampling for Gibbs
sampling. Applied Statistics 41, pp 337-348.