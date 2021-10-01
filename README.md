# GenericML
R implementation of Generic Machine Learning (Chernozhukov, V., Demirer, M., Duflo, E., &amp; Fernández-Val, I., 2021) using the `mlr3` framework. We intend to extend this implementation to a fully-fledged R package for the CRAN. Please note that this implementation is still work in progress and has not yet been thoroughy tested, so we cannot yet guarantee correctness. If you find a bug, please open an issue or let us know via email.


## Organization of this repo

* `examples` contains a demonstration of our implementation,
* `functions` contains all functions.

## TODO

- [x] Make it optional for user to supply propensity scores;
- [x] Implement stratified sampling for sample splitting (function done, but implementation in main procedure still needs to be done);
- [x] Make it optional to use arbitrary covariates in `X1`;
- [x] Make cluster-robust standard errors optional;
- [x] Implement fixed effects in the BLP/GATES regressions;
- [x] Make generic targets as differences optional;
- [x] Parallelize the main loop;
- [x] Implement deep neural networks as regression learner. Currently not supported by `mlr3` (NB: cannot be done in current version of `mlr3`);
- [ ] Add optional monotonization of the confindence bounds;
- [x] Error handling. Current implementation might throw an error in case of undersampling;
- [x] Implement S3 class structure and reorganize repo;
- [ ] Make an R package.

## Authors
Max Welz (m.welz@erasmusmc.nl), Andreas Alfons (aalfons@ese.eur.nl), and Mert Demirer (mdemirer@mit.edu).
