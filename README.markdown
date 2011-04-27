### OPPVASP -- Oslo Python Package for Vasp and Similar Packages

OPPVASP is a wrap-up of various pre- and post-processing scripts 
for ab initio codes, mainly VASP.

A warning: This code should not be considered stable, and it's 
structure is likely to change in the future. 
If you are just looking for a good python pre-/post-processing 
package for ab initio codes, try ASE: 
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
