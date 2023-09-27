# Meson Lattice Utilities (MLU)
aka
# Mike's Lattice Utilities (MLU)

## Mike's PhD thesis

Item | Description
| --- | ---
`libMLU.a` | Library containing helper functions (HDF5 wrappers, I/O, filename globbing, ...) 
`bootstrap`, `MultiFit`, ... | Generic analysis code  
`xml3pt` | Production jobs for generating semilep dataset
`Scripts` | All analysis jobs to turn correlators into form factors

## Dependencies

1. [Grid] Latest version
2. [Hadrons] commit 49fd078f7b8d43d5c6c2c4814e1a8987afc244f0
3. [GNU Scientific Library][gsl] (GSL), version 2.5 or later (I'm using 2.7 in production) or from [Mike's home page][MikeGSL]
4. [CERN Minuit2][minuit2] standalone version (see paragraph 3, "... can be downloaded [here]") or from [Mike's home page][MikeMinuit2]

[grid]: https://github.com/paboyle/Grid
[hadrons]: https://github.com/aportelli/Hadrons
[gsl]: https://www.gnu.org/software/gsl/
[minuit2]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/html/index.html
[here]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/Minuit2.tar.gz
[MikeMinuit2]: https://www2.ph.ed.ac.uk/~s1786208/Minuit2-5.34.14.tar.gz
[MikeGSL]: https://www2.ph.ed.ac.uk/~s1786208/gsl-2.7.tar.gz
