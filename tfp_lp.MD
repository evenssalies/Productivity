### A blend of ideas to put in form later

package st0537 ```prodest```, by Rovigatti and Mollisi (2018) with the Chile data file. 

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
-	predict ., omega after levpet adds \eta to productivity \omega. And there is nothing we can do to avoid this!
-	predict ., omega after prodest, however, just predicts \omega. To recover \omega+\eta, use the option ```[residual]```.

This distinction is not apparent in Rovigatti and Mollisi's (2018) paper since the authors calculate the production function elasticities, not productivity.
