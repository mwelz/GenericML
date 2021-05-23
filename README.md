# GenericML
R implementation of Generic Machine Learning (Chernozhukov, V., Demirer, M., Duflo, E., &amp; Fernández-Val, I., 2021) using the `mlr3` framework. We intend to extend this implementation to a fully-fledged R package for the CRAN.

## Organization of this repo

* `examples` contains a demonstration of our implementation,
* `functions` contains all functions.

## TODO

- [ ] Implement deep neural networks as regression learner. Currently not supported by `mlr3`;
- [ ] Implement the Horvitz-Thompson transformation;
- [ ] Add noise in case of low variation of the `S` and `B` learner;
- [ ] Add optional monotonization of the confindence bounds;
- [ ] Error handling. Current implementation might throw an error in case of undersampling;
- [ ] Make an R package.

## Authors
Max Welz (m.welz@erasmusmc.nl) and Andreas Alfons (aalfons@ese.eur.nl).
