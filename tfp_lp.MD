### A blend of ideas to put in form later

package st0537 ```prodest```, by Rovigatti and Molisi (2018) with the Chile data file. 

@article{rovigatti2018theory,
  title={Theory and practice of total-factor productivity estimation: The control function approach using Stata},
  author={Rovigatti, Gabriele and Mollisi, Vincenzo},
  journal={The Stata Journal},
  volume={18},
  number={3},
  pages={618--662},
  year={2018},
  publisher={SAGE Publications Sage CA: Los Angeles, CA}
}

A distinction between levpet and prodest is:
-	predict ., omega after levpet adds \eta to productivity \omega.
-	predict ., omega after prodest just predivts \omega. To recover \omega+\eta, use the residual.

À bientôt.
