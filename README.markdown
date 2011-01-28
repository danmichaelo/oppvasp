### OPPVASP -- Oslo Python Package for Vasp and Similar Packages

OPPVASP is an attempt of a Python package wrapping together various 
pre- and post-processing scripts for VASP and perhaps other 
Ab Initio codes.

Please note: The code is in a very early development stage, and no 
parts of it can be considered stable. It may be useful for inspiration,
but if you are looking for a good python pre-/post-processing pacakge, 
I would recommend ASE: 
<https://wiki.fysik.dtu.dk/ase/>.

Some curent features:

* Simple setup of batch jobs for e.g. convergence tests
* Simple band structure parser for VASP and Quantum Espresso output
* Band analyzer that checks for minima, maxima, bandwidth...
* Band structure plotting using matplotlib 

### Installing

Requirements: numpy, scipy, matplotlib, scitools, lxml 
Optional: progressbar

    git clone git://github.com/danmichaelo/oppvasp.git

### Author

dmheggo@student.matnat.uio.no
