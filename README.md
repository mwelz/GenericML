# GenericML
R implementation of Generic Machine Learning (Chernozhukov, V., Demirer, M., Duflo, E., &amp; Fernández-Val, I., 2021) using the `mlr3` framework. We intend to extend this implementation to a fully-fledged R package for the CRAN.

## Organization of this repo

* `examples` contains a demonstration of our implementation,
* `functions` contains all functions.

## TODO

- [x] Make it optional for user to supply propensity scores;
- [x] Implement stratified sampling for sample splitting (function done, but implementation in main procedure still needs to be done);
- [ ] Make it optional to use arbitrary covariates in `X1`;
- [ ] Make it optional that BLP can also use fixed effects and clustering;
- [ ] Make generic targets as differences optional;
- [ ] Parallelize the main loop;
- [ ] Implement deep neural networks as regression learner. Currently not supported by `mlr3`;
- [ ] Add optional monotonization of the confindence bounds;
- [ ] Error handling. Current implementation might throw an error in case of undersampling;
- [ ] Implement S3 class structure and reorganize repo;
- [ ] Make an R package.

## Authors
Max Welz (m.welz@erasmusmc.nl), Andreas Alfons (aalfons@ese.eur.nl), and Mert Demirer (mdemirer@mit.edu).
