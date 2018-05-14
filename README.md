# AvERA_ISW

This repository contains code and data that allows the reproduction of the results of our paper
### Beck et al. (2018): The integrated Sachs-Wolfe effect in the AvERA cosmology (https://arxiv.org/abs/1801.08566).
Authors of publications that make use of our code or data are required to reference our paper arXiv:1801.08566 [astro-ph.CO] in their works.

The Jupyter notebook requires several other packages to run: PyCAMB, HealPy, Numpy, Matplotlib and Scipy. With these installed, however, simply running the entire notebook will generate most of the results and figures that we present in our paper, with the exception of the raytracing.

The Millennium XXL simulation snapshot ID 14 (at z=8.55) has to be obtained from its creators (Angulo et al., 2012), and copied to the folder *data/MillenniumXXL*. It was originally downloaded from the location:
http://galformod.mpa-garching.mpg.de/public/mxxl/densityfields/density_1024_014.dat

We note that loading this snapshot into the notebook requires 2*4 GB of memory. As it is only used for determining the scaling of the LCDM initial power spectrum relative to Millennium, this step may be omitted without affecting other parts of the notebook.

Additionally, the CMB power spectrum measured by the Planck Collaboration (2016b) has to be downloaded from its source, and copied to the folder *data/Planck*. Original download location: https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/cosmoparams/COM_PowerSpect_CMB-base-plikHM-TT-lowTEB-minimum-theory_R2.02.txt

Regarding the raytracing calculations, the Python code is provided for reference (.py files in the root folder), but these are meant to be specialized for (and run in parallel) on multi-core systems with ~128 GB of memory. Please contact Robert Beck for assistance if you wish to reproduce this part of our analysis, or if you have any questions about other elements of the repository.
