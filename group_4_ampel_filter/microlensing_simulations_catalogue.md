Here is the complete catalogue of simulation tools, converted into a structured markdown file.

# Simulation Tools Catalogue

This document catalogues various software tools used for astrophysical simulations, particularly focusing on microlensing and surveys like Rubin and Roman.

-----

## Survey Operations & Cadence Simulation

Tools for generating and evaluating survey schedules and observing strategies.

### rubin\_scheduler (FBS)

> Generating realistic 10-yr visit schedules under constraints

  * **Functional Classification:** `Survey Operations & Cadence Simulation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/lsst/rubin_sim`](https://www.google.com/search?q=%5Bhttps://github.com/lsst/rubin_sim%5D\(https://github.com/lsst/rubin_sim\))
  * **Documentation:** [`https://rubin-scheduler.lsst.io`](https://www.google.com/search?q=%5Bhttps://rubin-scheduler.lsst.io%5D\(https://rubin-scheduler.lsst.io\))
  * **Installation:** [`https://rubin-scheduler.lsst.io/installation.html`](https://www.google.com/search?q=%5Bhttps://rubin-scheduler.lsst.io/installation.html%5D\(https://rubin-scheduler.lsst.io/installation.html\))
    ```bash
    pip install rubin-scheduler, scheduler_download_data
    ```
  * **Minimal Example:** [`https://github.com/lsst/rubin_sim_notebooks/tree/main/scheduler`](https://www.google.com/search?q=%5Bhttps://github.com/lsst/rubin_sim_notebooks/tree/main/scheduler%5D\(https://github.com/lsst/rubin_sim_notebooks/tree/main/scheduler\))

### rubin\_sim

> Cadence metrics, sky models, throughput; successor to MAF/CatSim

  * **Functional Classification:** `Survey Operations & Cadence Simulation; Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/lsst/rubin_scheduler`](https://www.google.com/search?q=%5Bhttps://github.com/lsst/rubin_scheduler%5D\(https://github.com/lsst/rubin_scheduler\))
  * **Documentation:** [`https://rubin-sim.lsst.io`](https://www.google.com/search?q=%5Bhttps://rubin-sim.lsst.io%5D\(https://rubin-sim.lsst.io\))
  * **Installation:** [`https://rubin-sim.lsst.io/installation.html`](https://www.google.com/search?q=%5Bhttps://rubin-sim.lsst.io/installation.html%5D\(https://rubin-sim.lsst.io/installation.html\))
    ```bash
    pip install rubin-sim, scheduler_download_data, rs_download_data
    ```
  * **Minimal Example:** [`https://github.com/lsst/rubin_sim_notebooks/tree/main/`](https://www.google.com/search?q=%5Bhttps://github.com/lsst/rubin_sim_notebooks/tree/main/%5D\(https://github.com/lsst/rubin_sim_notebooks/tree/main/\))

### GBTDS Optimizer

> Tool for optimizing yields or metrics for the Roman Galactic Bulge Time Domain Survey

  * **Functional Classification:** `Survey Operations & Cadence Simulation`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/mtpenny/gbtds_optimizer`](https://www.google.com/search?q=%5Bhttps://github.com/mtpenny/gbtds_optimizer%5D\(https://github.com/mtpenny/gbtds_optimizer\))
  * **Documentation:** [`https://github.com/mtpenny/gbtds_optimizer/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/mtpenny/gbtds_optimizer/blob/main/README.md%5D\(https://github.com/mtpenny/gbtds_optimizer/blob/main/README.md\))
  * **Installation:** [`https://github.com/mtpenny/gbtds_optimizer/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/mtpenny/gbtds_optimizer/blob/main/README.md%5D\(https://github.com/mtpenny/gbtds_optimizer/blob/main/README.md\))
    ```bash
    git clone https://github.com/mtpenny/gbtds_optimizer.git
    ```
  * **Minimal Example:**
    ```bash
    python gbtds_optimizer.py ffp_normmap_m+00v3_rate.yield.csv 14.7315 42.56 \
           field_layouts/layout_7f_3_gal-center.centers \
       --lrange 2.2 -2.2 --brange -2.2 2.2 \
       --lstep 0.2 --bstep 0.2 \
       --cadence-bounds 7.0 15.0 \
       --nread-bounds 10 40 \
       --output-root ffp_normmap_m+00v3.layout_7f_3
    ```

### OpSim (legacy) (replaced by rubin\_sim)

> Historic reference for schedules

  * **Functional Classification:** `Survey Operations & Cadence Simulation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `C++/Python Wrapper`
  * **Source Code:** [`https://github.com/lsst/rubin_sim`](https://www.google.com/search?q=%5Bhttps://github.com/lsst/rubin_sim%5D\(https://github.com/lsst/rubin_sim\))
  * **Documentation:** [`https://www.lsst.org/scientists/simulations/opsim`](https://www.google.com/search?q=%5Bhttps://www.lsst.org/scientists/simulations/opsim%5D\(https://www.lsst.org/scientists/simulations/opsim\))
  * **Installation:** `Replaced by rubin_sim`

-----

## Astrophysical Scene & Catalog Generation

Tools for creating synthetic stellar and galactic populations and catalogs.

### rubin\_sim

> Cadence metrics, sky models, throughput; successor to MAF/CatSim

  * **Functional Classification:** `Survey Operations & Cadence Simulation; Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/lsst/rubin_scheduler`](https://www.google.com/search?q=%5Bhttps://github.com/lsst/rubin_scheduler%5D\(https://github.com/lsst/rubin_scheduler\))
  * **Documentation:** [`https://rubin-sim.lsst.io`](https://www.google.com/search?q=%5Bhttps://rubin-sim.lsst.io%5D\(https://rubin-sim.lsst.io\))
  * **Installation:** [`https://rubin-sim.lsst.io/installation.html`](https://www.google.com/search?q=%5Bhttps://rubin-sim.lsst.io/installation.html%5D\(https://rubin-sim.lsst.io/installation.html\))
    ```bash
    pip install rubin-sim, scheduler_download_data, rs_download_data
    ```
  * **Minimal Example:** [`https://github.com/lsst/rubin_sim_notebooks/tree/main/`](https://www.google.com/search?q=%5Bhttps://github.com/lsst/rubin_sim_notebooks/tree/main/%5D\(https://github.com/lsst/rubin_sim_notebooks/tree/main/\))

### PopSyCLE

> Population Synthesis for Compact-object Lensing Events is a code to simulate a model of the Milky Way including compact objects and multiple systems and perform a mock microlensing survey. You can use it to put black hole candidates into context and to understand the effect of Galactic properties on photometric and astrometric microlensing simulation distributions among many other applications.

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/jluastro/PopSyCLE/tree/main`](https://www.google.com/search?q=%5Bhttps://github.com/jluastro/PopSyCLE/tree/main%5D\(https://github.com/jluastro/PopSyCLE/tree/main\))
  * **Documentation:** [`https://popsycle.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://popsycle.readthedocs.io/en/latest/%5D\(https://popsycle.readthedocs.io/en/latest/\))
  * **Installation:** [`https://github.com/jluastro/PopSyCLE/blob/main/README.rst`](https://www.google.com/search?q=%5Bhttps://github.com/jluastro/PopSyCLE/blob/main/README.rst%5D\(https://github.com/jluastro/PopSyCLE/blob/main/README.rst\))
    ```bash
    git clone git@github.com:jluastro/PopSyCLE.git
    echo "export PYTHONPATH=$PWD/PopSyCLE:$PYTHONPATH" >> ~/.bashrc
    ```
  * **Minimal Example:** [`https://github.com/jluastro/PopSyCLE/blob/main/docs/PopSyCLE_example_run.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/jluastro/PopSyCLE/blob/main/docs/PopSyCLE_example_run.ipynb%5D\(https://github.com/jluastro/PopSyCLE/blob/main/docs/PopSyCLE_example_run.ipynb\))
  * **Associated Publication:** [`https://iopscience.iop.org/article/10.3847/1538-4357/ab5fd3`](https://www.google.com/search?q=%5Bhttps://iopscience.iop.org/article/10.3847/1538-4357/ab5fd3%5D\(https://iopscience.iop.org/article/10.3847/1538-4357/ab5fd3\))

### SynthPop / Synthpop

> Synthpop is an object-oriented, modular Python framework for generating synthetic population models. It generates a star catalog following the specified model and configuration.

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Star catalogs`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/synthpop-galaxy/synthpop`](https://www.google.com/search?q=%5Bhttps://github.com/synthpop-galaxy/synthpop%5D\(https://github.com/synthpop-galaxy/synthpop\))
  * **Documentation:** [`https://synthpop.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://synthpop.readthedocs.io/en/latest/%5D\(https://synthpop.readthedocs.io/en/latest/\))
  * **Installation:** [`https://github.com/synthpop-galaxy/synthpop/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/synthpop-galaxy/synthpop/blob/main/README.md%5D\(https://github.com/synthpop-galaxy/synthpop/blob/main/README.md\))
    ```bash
    pip install git+https://github.com/synthpop-galaxy/synthpop.git
    python -m synthpop.migrate_interactive_part
    ```
  * **Minimal Example:** `The example is to prepare a configuration file, the process for which is given in the Readme document.`
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2025AJ....169..317K/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2025AJ....169..317K/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2025AJ....169..317K/abstract\))

### genulens (Koshimoto et al.)

> genulens, which stands for "generate microlensing", is a tool to simulate microlensing events using Monte Carlo simulation of the Galactic model developed by Koshimoto, Baba & Bennett (2021), ApJ, 917, 78. The Galactic model is optimized for the bulge direction, and it is best to be used for analyzing microlensing events in |l| \< 10 deg. and |b| \< 7 deg.

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Jupyter Notebook/ C++`
  * **Source Code:** [`https://github.com/nkoshimoto/genulens`](https://www.google.com/search?q=%5Bhttps://github.com/nkoshimoto/genulens%5D\(https://github.com/nkoshimoto/genulens\))
  * **Installation:** [`https://github.com/nkoshimoto/genulens/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/nkoshimoto/genulens/blob/main/README.md%5D\(https://github.com/nkoshimoto/genulens/blob/main/README.md\))
    ```bash
    gsl-config --libs
    git clone https://github.com/nkoshimoto/genulens.git
    ```
  * **Minimal Example:** [`https://github.com/nkoshimoto/genulens/blob/main/genulens_samples.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/nkoshimoto/genulens/blob/main/genulens_samples.ipynb%5D\(https://github.com/nkoshimoto/genulens/blob/main/genulens_samples.ipynb\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract\)), [`http://doi.org/10.5281/zenodo.4784948`](https://www.google.com/search?q=%5Bhttp://doi.org/10.5281/zenodo.4784948%5D\(http://doi.org/10.5281/zenodo.4784948\))

### TRILEGAL

> Mock stellar catalogs for LSST/Rubin fields

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Rubin/Star catalogs`
  * **Language:** `Fortran/C (web backend)`
  * **Source Code:** [`https://stev.oapd.inaf.it/cgi-bin/trilegal`](https://www.google.com/search?q=%5Bhttps://stev.oapd.inaf.it/cgi-bin/trilegal%5D\(https://stev.oapd.inaf.it/cgi-bin/trilegal\))
  * **Documentation:** [`https://datalab.noirlab.edu/data/lsst-sim`](https://www.google.com/search?q=%5Bhttps://datalab.noirlab.edu/data/lsst-sim%5D\(https://datalab.noirlab.edu/data/lsst-sim\))
  * **Installation:** `Web form`
  * **Minimal Example:** `Use web form to generate catalog; export as CSV/FITS`
  * **Associated Publication:** [`https://iopscience.iop.org/article/10.3847/1538-4365/ac7be6`](https://www.google.com/search?q=%5Bhttps://iopscience.iop.org/article/10.3847/1538-4365/ac7be6%5D\(https://iopscience.iop.org/article/10.3847/1538-4365/ac7be6\))

### Besançon Galaxy Model (BGM)

> The model is a powerful tool to constrain either evolutionary scenarii or galactic structure hypothesis through the comparison between model predictions and a large variety of observational constraints such as star counts, photometry or astrometry.

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Star catalogs`
  * **Language:** `Fortran/C (web backend)`
  * **Source Code:** [`https://model.obs-besancon.fr/`](https://www.google.com/search?q=%5Bhttps://model.obs-besancon.fr/%5D\(https://model.obs-besancon.fr/\))
  * **Documentation:** [`https://model.obs-besancon.fr/`](https://www.google.com/search?q=%5Bhttps://model.obs-besancon.fr/%5D\(https://model.obs-besancon.fr/\))
  * **Installation:** `Web interface`
  * **Minimal Example:** `Define sky field & magnitude limits; export catalog`
  * **Associated Publication:** [`https://www.aanda.org/articles/aa/full_html/2018/12/aa33501-18/aa33501-18.html`](https://www.google.com/search?q=%5Bhttps://www.aanda.org/articles/aa/full_html/2018/12/aa33501-18/aa33501-18.html%5D\(https://www.aanda.org/articles/aa/full_html/2018/12/aa33501-18/aa33501-18.html\))

### Galaxia

> Galaxia is a code for generating a synthetic model of the galaxy. The input model can be analytical or one obtained from N-body simulations. The code outputs a catalog of stars according to user specified color magnitude limits.

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Star catalogs`
  * **Language:** `C++`
  * **Source Code:** [`https://galaxia.sourceforge.net/`](https://www.google.com/search?q=%5Bhttps://galaxia.sourceforge.net/%5D\(https://galaxia.sourceforge.net/\))
  * **Documentation:** [`https://galaxia.sourceforge.net/Galaxia3pub.html`](https://www.google.com/search?q=%5Bhttps://galaxia.sourceforge.net/Galaxia3pub.html%5D\(https://galaxia.sourceforge.net/Galaxia3pub.html\))
  * **Installation:** [`https://sourceforge.net/projects/galaxia/files/`](https://www.google.com/search?q=%5Bhttps://sourceforge.net/projects/galaxia/files/%5D\(https://sourceforge.net/projects/galaxia/files/\))
  * **Minimal Example:** `./galaxia <config> to generate mock catalog` (See `Examples/` folder after download)
  * **Associated Publication:** [`https://iopscience.iop.org/article/10.1088/0004-637X/730/1/3`](https://www.google.com/search?q=%5Bhttps://iopscience.iop.org/article/10.1088/0004-637X/730/1/3%5D\(https://iopscience.iop.org/article/10.1088/0004-637X/730/1/3\))

### SPISEA

> SPISEA is an open-source python package that generates single-age, single-metallicity populations (i.e. star clusters).

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Star catalogs`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/astropy/SPISEA`](https://www.google.com/search?q=%5Bhttps://github.com/astropy/SPISEA%5D\(https://github.com/astropy/SPISEA\))
  * **Documentation:** [`https://spisea.readthedocs.io/en/stable/index.html`](https://www.google.com/search?q=%5Bhttps://spisea.readthedocs.io/en/stable/index.html%5D\(https://spisea.readthedocs.io/en/stable/index.html\))
  * **Installation:** [`https://spisea.readthedocs.io/en/stable/getting_started.html`](https://www.google.com/search?q=%5Bhttps://spisea.readthedocs.io/en/stable/getting_started.html%5D\(https://spisea.readthedocs.io/en/stable/getting_started.html\))
    ```bash
    git clone https://github.com/astropy/SPISEA.git
    ```
  * **Minimal Example:** [`https://github.com/astropy/SPISEA/blob/main/docs/Quick_Start_Make_Cluster.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/astropy/SPISEA/blob/main/docs/Quick_Start_Make_Cluster.ipynb%5D\(https://github.com/astropy/SPISEA/blob/main/docs/Quick_Start_Make_Cluster.ipynb\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2020arXiv200606691H/abstract\))

### CatSim (legacy)

> The catalog simulation framework (CatSim) comprises a database of stars, galaxies, and Solar System objects that can be queried as a function of position and survey time. These simulations incorporate variability for a subset of the stellar sources and for the AGN in the galaxy catalog.

  * **Functional Classification:** `Astrophysical Scene & Catalog Generation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `Python/SQL`
  * **Source Code:** [`https://www.lsst.org/scientists/simulations/catsim`](https://www.google.com/search?q=%5Bhttps://www.lsst.org/scientists/simulations/catsim%5D\(https://www.lsst.org/scientists/simulations/catsim\))
  * **Documentation:** [`https://www.lsst.org/scientists/simulations/catsim`](https://www.google.com/search?q=%5Bhttps://www.lsst.org/scientists/simulations/catsim%5D\(https://www.lsst.org/scientists/simulations/catsim\))
  * **Installation:** `https://confluence.lsstcorp.org/display/SIM/Catalogs+and+MAF / Use rubin_sim for modern use`
  * **Minimal Example:** [`https://github.com/uwssg/LSST-Tutorials/blob/master/CatSim/CatSimTutorial_SimulationsAHM_1503.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/uwssg/LSST-Tutorials/blob/master/CatSim/CatSimTutorial_SimulationsAHM_1503.ipynb%5D\(https://github.com/uwssg/LSST-Tutorials/blob/master/CatSim/CatSimTutorial_SimulationsAHM_1503.ipynb\))

-----

## Image & Instrument Simulation

Tools for generating realistic astronomical images, including instrument and atmospheric effects.

### Roman iSim

> Simulate images for the Nancy Grace Roman Space Telescope Time-Domain Survey.

  * **Functional Classification:** `Image & Instrument Simulation`
  * **Domain/Category:** `Roman`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/MichaelDAlbrow/RomanISim-simulate`](https://www.google.com/search?q=%5Bhttps://github.com/MichaelDAlbrow/RomanISim-simulate%5D\(https://github.com/MichaelDAlbrow/RomanISim-simulate\))
  * **Documentation:** [`https://romanisim.readthedocs.io/`](https://www.google.com/search?q=%5Bhttps://romanisim.readthedocs.io/%5D\(https://romanisim.readthedocs.io/\)), [`https://stpsf.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://stpsf.readthedocs.io/en/latest/%5D\(https://stpsf.readthedocs.io/en/latest/\))
  * **Installation:** [`https://github.com/MichaelDAlbrow/RomanISim-simulate/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/MichaelDAlbrow/RomanISim-simulate/blob/main/README.md%5D\(https://github.com/MichaelDAlbrow/RomanISim-simulate/blob/main/README.md\))
    ```bash
    git clone https://github.com/spacetelescope/romanisim
    git clone https://github.com/spacetelescope/webbpsf
    ```
  * **Minimal Example:**
    ```bash
    python -u make_synthpop_image.py synthpop_config18.json >& test18.log &
    ```

### imSim

> imSim is a software package that simulates the LSST telescope and survey. It produces simulated images from the 3.25 Gigapixel camera which are suitable to be processed through the LSST Data Management pipeline. imSim takes as an input a catalog of astronomical sources along with information about how the light is distorted on the way to Earth including lensing and extinction information. The images which are produced include the systematic effects of the atmosphere, optics and sensors on the observed PSF.

  * **Functional Classification:** `Image & Instrument Simulation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/LSSTDESC/imSim`](https://www.google.com/search?q=%5Bhttps://github.com/LSSTDESC/imSim%5D\(https://github.com/LSSTDESC/imSim\))
  * **Documentation:** [`https://lsstdesc.org/imSim/index.html`](https://www.google.com/search?q=%5Bhttps://lsstdesc.org/imSim/index.html%5D\(https://lsstdesc.org/imSim/index.html\))
  * **Installation:** [`https://lsstdesc.org/imSim/install.html`](https://www.google.com/search?q=%5Bhttps://lsstdesc.org/imSim/install.html%5D\(https://lsstdesc.org/imSim/install.html\))
  * **Minimal Example:** [`https://lsstdesc.org/imSim/usage.html`](https://www.google.com/search?q=%5Bhttps://lsstdesc.org/imSim/usage.html%5D\(https://lsstdesc.org/imSim/usage.html\))

### GalSim

> GalSim is open-source software for simulating images of astronomical objects (stars, galaxies) in a variety of ways. The bulk of the calculations are carried out in C++, and the user interface is in Python. In addition, the code can operate directly on "config" files, for those users who prefer not to work in Python.

  * **Functional Classification:** `Image & Instrument Simulation`
  * **Domain/Category:** `General (used by imSim)`
  * **Language:** `C++/Python Wrapper`
  * **Source Code:** [`https://github.com/GalSim-developers/GalSim`](https://www.google.com/search?q=%5Bhttps://github.com/GalSim-developers/GalSim%5D\(https://github.com/GalSim-developers/GalSim\))
  * **Documentation:** [`http://galsim-developers.github.io/GalSim/`](https://www.google.com/search?q=%5Bhttp://galsim-developers.github.io/GalSim/%5D\(http://galsim-developers.github.io/GalSim/\))
  * **Installation:** [`https://github.com/GalSim-developers/GalSim/blob/releases/2.7/README.rst`](https://www.google.com/search?q=%5Bhttps://github.com/GalSim-developers/GalSim/blob/releases/2.7/README.rst%5D\(https://github.com/GalSim-developers/GalSim/blob/releases/2.7/README.rst\))
    ```bash
    pip install galsim | conda install -c conda-forge galsim
    ```
  * **Minimal Example:** [`https://github.com/GalSim-developers/GalSim/tree/releases/2.7/examples`](https://www.google.com/search?q=%5Bhttps://github.com/GalSim-developers/GalSim/tree/releases/2.7/examples%5D\(https://github.com/GalSim-developers/GalSim/tree/releases/2.7/examples\))
  * **Associated Publication:** [`http://adsabs.harvard.edu/abs/2015A%26C....10..121R`](https://www.google.com/search?q=%5Bhttp://adsabs.harvard.edu/abs/2015A%2526C....10..121R%5D\(http://adsabs.harvard.edu/abs/2015A%2526C....10..121R\))

### ts\_imsim

> Package to create Rubin Active Optics Simulations with imSim

  * **Functional Classification:** `Image & Instrument Simulation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/lsst-ts/ts_imsim`](https://www.google.com/search?q=%5Bhttps://github.com/lsst-ts/ts_imsim%5D\(https://github.com/lsst-ts/ts_imsim\))
  * **Documentation:** [`https://ts-imsim.lsst.io/`](https://www.google.com/search?q=%5Bhttps://ts-imsim.lsst.io/%5D\(https://ts-imsim.lsst.io/\))
  * **Installation:** [`https://ts-imsim.lsst.io/user-guide/user-guide.html`](https://www.google.com/search?q=%5Bhttps://ts-imsim.lsst.io/user-guide/user-guide.html%5D\(https://ts-imsim.lsst.io/user-guide/user-guide.html\))
  * **Minimal Example:** [`https://github.com/lsst-ts/ts_imsim/blob/develop/doc/user-guide/user-guide.rst`](https://www.google.com/search?q=%5Bhttps://github.com/lsst-ts/ts_imsim/blob/develop/doc/user-guide/user-guide.rst%5D\(https://github.com/lsst-ts/ts_imsim/blob/develop/doc/user-guide/user-guide.rst\))

### PhoSim (historical)

> The Photon Simulator (PhoSim) is a set of fast photon Monte Carlo codes used to calculate the physics of the atmosphere and a telescope and camera in order to simulate realistic astronomical images. It does this using modern numerical techniques applied to comprehensive physical models. PhoSim generates images by collecting photons into pixels.

  * **Functional Classification:** `Image & Instrument Simulation`
  * **Domain/Category:** `Rubin/LSST`
  * **Language:** `C++/Python wrapper`
  * **Source Code:** [`https://www.lsst.org/scientists/simulations/phosim`](https://www.google.com/search?q=%5Bhttps://www.lsst.org/scientists/simulations/phosim%5D\(https://www.lsst.org/scientists/simulations/phosim\))
  * **Documentation:** [`https://www.phosim.org/`](https://www.google.com/search?q=%5Bhttps://www.phosim.org/%5D\(https://www.phosim.org/\)), [`https://www.phosim.org/documentation`](https://www.google.com/search?q=%5Bhttps://www.phosim.org/documentation%5D\(https://www.phosim.org/documentation\))
  * **Installation:** [`https://www.phosim.org/download`](https://www.google.com/search?q=%5Bhttps://www.phosim.org/download%5D\(https://www.phosim.org/download\))
  * **Minimal Example:** [`https://www.phosim.org/tutorials`](https://www.google.com/search?q=%5Bhttps://www.phosim.org/tutorials%5D\(https://www.phosim.org/tutorials\))
  * **Associated Publication:** [`https://www.phosim.org/technical`](https://www.google.com/search?q=%5Bhttps://www.phosim.org/technical%5D\(https://www.phosim.org/technical\))

-----

## Microlensing Simulation & Modeling

Tools specifically designed to simulate, model, and fit gravitational microlensing events.

### BAGLE

> BAGLE is a python package used to model gravitational microlensing events both photometrically and astrometrically.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/MovingUniverseLab/BAGLE_Microlensing`](https://www.google.com/search?q=%5Bhttps://github.com/MovingUniverseLab/BAGLE_Microlensing%5D\(https://github.com/MovingUniverseLab/BAGLE_Microlensing\))
  * **Documentation:** [`https://bagle.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://bagle.readthedocs.io/en/latest/%5D\(https://bagle.readthedocs.io/en/latest/\))
  * **Installation:** [`https://github.com/MovingUniverseLab/BAGLE_Microlensing/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/MovingUniverseLab/BAGLE_Microlensing/blob/main/README.md%5D\(https://github.com/MovingUniverseLab/BAGLE_Microlensing/blob/main/README.md\))
    ```bash
    pip install BAGLE
    # See docs for full dependency list (numpy, astropy, matplotlib, etc.)
    ```
  * **Minimal Example:** [`https://github.com/MovingUniverseLab/BAGLE_Microlensing/blob/main/BAGLE_TUTORIAL.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/MovingUniverseLab/BAGLE_Microlensing/blob/main/BAGLE_TUTORIAL.ipynb%5D\(https://github.com/MovingUniverseLab/BAGLE_Microlensing/blob/main/BAGLE_TUTORIAL.ipynb\))
  * **Associated Publication:** [`https://bagle.readthedocs.io/en/latest/citation.html`](https://www.google.com/search?q=%5Bhttps://bagle.readthedocs.io/en/latest/citation.html%5D\(https://bagle.readthedocs.io/en/latest/citation.html\))

### MulensModel

> MulensModel is a package for modeling microlensing events. MulensModel can generate a microlensing light curve for a given set of microlensing parameters, fit that light curve to some data, and return a chi2 value. That chi2 can then be input into an arbitrary likelihood function to find the best-fit parameters.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/rpoleski/MulensModel`](https://www.google.com/search?q=%5Bhttps://github.com/rpoleski/MulensModel%5D\(https://github.com/rpoleski/MulensModel\))
  * **Documentation:** [`https://rpoleski.github.io/MulensModel/`](https://www.google.com/search?q=%5Bhttps://rpoleski.github.io/MulensModel/%5D\(https://rpoleski.github.io/MulensModel/\))
  * **Installation:** [`https://rpoleski.github.io/MulensModel/install.html`](https://www.google.com/search?q=%5Bhttps://rpoleski.github.io/MulensModel/install.html%5D\(https://rpoleski.github.io/MulensModel/install.html\))
    ```bash
    pip install MulensModel
    ```
  * **Minimal Example:** [`https://github.com/rpoleski/MulensModel/blob/master/documents/examples_list.md`](https://www.google.com/search?q=%5Bhttps://github.com/rpoleski/MulensModel/blob/master/documents/examples_list.md%5D\(https://github.com/rpoleski/MulensModel/blob/master/documents/examples_list.md\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2019A%26C....26...35P/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2019A%2526C....26...35P/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2019A%2526C....26...35P/abstract\))

### pyLIMA

> pyLIMA is the first microlensing analysis open-source software, primarly designed to fit real data.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/ebachelet/pyLIMA/tree/master`](https://www.google.com/search?q=%5Bhttps://github.com/ebachelet/pyLIMA/tree/master%5D\(https://github.com/ebachelet/pyLIMA/tree/master\))
  * **Documentation:** [`https://pylima.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://pylima.readthedocs.io/en/latest/%5D\(https://pylima.readthedocs.io/en/latest/\))
  * **Installation:** [`https://pylima.readthedocs.io/en/latest/source/Installation.html`](https://www.google.com/search?q=%5Bhttps://pylima.readthedocs.io/en/latest/source/Installation.html%5D\(https://pylima.readthedocs.io/en/latest/source/Installation.html\))
    ```bash
    pip install pyLIMA
    ```
  * **Minimal Example:** [`https://github.com/ebachelet/pyLIMA/tree/master/examples`](https://www.google.com/search?q=%5Bhttps://github.com/ebachelet/pyLIMA/tree/master/examples%5D\(https://github.com/ebachelet/pyLIMA/tree/master/examples\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2017AJ....154..203B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2017AJ....154..203B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2017AJ....154..203B/abstract\))

### RTModel

> RTModel is a package for modeling and interpreting microlensing events. It uses photometric and/or astrometric time series collected from ground and/or space telescopes to propose one or more possible models... All models include the finite size of the source(s).

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `C++/Python`
  * **Source Code:** [`https://github.com/valboz/RTModel/tree/main`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/RTModel/tree/main%5D\(https://github.com/valboz/RTModel/tree/main\))
  * **Documentation:** [`https://projects.phys.unisa.it/GravitationAstrophysics/RTModel.htm`](https://www.google.com/search?q=%5Bhttps://projects.phys.unisa.it/GravitationAstrophysics/RTModel.htm%5D\(https://projects.phys.unisa.it/GravitationAstrophysics/RTModel.htm\))
  * **Installation:** [`https://github.com/valboz/RTModel/blob/main/docs/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/RTModel/blob/main/docs/README.md%5D\(https://github.com/valboz/RTModel/blob/main/docs/README.md\))
    ```bash
    pip install RTModel
    ```
  * **Minimal Example:** [`https://github.com/valboz/RTModel/blob/main/jupyter`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/RTModel/blob/main/jupyter%5D\(https://github.com/valboz/RTModel/blob/main/jupyter\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..83B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2024A%2526A...688A..83B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2024A%2526A...688A..83B/abstract\))

### VBMicrolensing (new)

> VBMicrolensing is a tool for efficient computation in gravitational microlensing events using the advanced contour integration method, supporting single, binary and multiple lenses... designed for... Magnification... Centroid... Critical curves and caustics... Complete light curves including several higher order effects.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `C++ /Python Wrapper`
  * **Source Code:** [`https://github.com/valboz/VBMicrolensing`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBMicrolensing%5D\(https://github.com/valboz/VBMicrolensing\))
  * **Documentation:** [`https://pypi.org/project/VBMicrolensing/`](https://www.google.com/search?q=%5Bhttps://pypi.org/project/VBMicrolensing/%5D\(https://pypi.org/project/VBMicrolensing/\)), [`https://github.com/valboz/VBMicrolensing/blob/main/docs/python/readme.md`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBMicrolensing/blob/main/docs/python/readme.md%5D\(https://github.com/valboz/VBMicrolensing/blob/main/docs/python/readme.md\))
  * **Installation:** [`https://github.com/valboz/VBMicrolensing/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBMicrolensing/blob/main/README.md%5D\(https://github.com/valboz/VBMicrolensing/blob/main/README.md\))
    ```bash
    pip install VBMicrolensing
    ```
  * **Minimal Example:** [`https://github.com/valboz/VBMicrolensing/blob/main/examples`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBMicrolensing/blob/main/examples%5D\(https://github.com/valboz/VBMicrolensing/blob/main/examples\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract\))

### eesunhong

> Image-centered ray-shooting method for microlensing

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Fortran/Python Wrapper`
  * **Source Code:** [`https://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract\))
  * **Documentation:** [`https://eesunhong.readthedocs.io/en/latest/index.html`](https://www.google.com/search?q=%5Bhttps://eesunhong.readthedocs.io/en/latest/index.html%5D\(https://eesunhong.readthedocs.io/en/latest/index.html\))
  * **Installation:** [`https://eesunhong.readthedocs.io/en/latest/index.html`](https://www.google.com/search?q=%5Bhttps://eesunhong.readthedocs.io/en/latest/index.html%5D\(https://eesunhong.readthedocs.io/en/latest/index.html\))
    ```bash
    pip install eesunhong
    ```
  * **Minimal Example:** [`https://eesunhong.readthedocs.io/en/latest/user_guides/index.html`](https://www.google.com/search?q=%5Bhttps://eesunhong.readthedocs.io/en/latest/user_guides/index.html%5D\(https://eesunhong.readthedocs.io/en/latest/user_guides/index.html\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/1996ApJ...472..660B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/1996ApJ...472..660B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/1996ApJ...472..660B/abstract\)), [`https://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract\))

### triplelens

> Light curves and image positions for triple lens systems

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `C/Python Wrapper`
  * **Source Code:** [`https://github.com/rkkuang/triplelens/tree/master`](https://www.google.com/search?q=%5Bhttps://github.com/rkkuang/triplelens/tree/master%5D\(https://github.com/rkkuang/triplelens/tree/master\))
  * **Documentation:** [`https://github.com/rkkuang/triplelens/tree/master`](https://www.google.com/search?q=%5Bhttps://github.com/rkkuang/triplelens/tree/master%5D\(https://github.com/rkkuang/triplelens/tree/master\))
  * **Installation:** [`https://github.com/rkkuang/triplelens/blob/master/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/rkkuang/triplelens/blob/master/README.md%5D\(https://github.com/rkkuang/triplelens/blob/master/README.md\))
    ```bash
    # The code needs to be compiled with the "make" command
    ```
  * **Minimal Example:** [`https://github.com/rkkuang/triplelens/blob/master/notebooks`](https://www.google.com/search?q=%5Bhttps://github.com/rkkuang/triplelens/blob/master/notebooks%5D\(https://github.com/rkkuang/triplelens/blob/master/notebooks\))
  * **Associated Publication:** [`https://doi.org/10.1093/mnras/stab509`](https://www.google.com/search?q=%5Bhttps://doi.org/10.1093/mnras/stab509%5D\(https://doi.org/10.1093/mnras/stab509\))

### sfit\_minimizer

> sfit\_minimizer is a gradient-type minimization algorithm known to work particularly well for point lens microlensing light curves. This algorithm generalizes Simpson's idea that a 1-D function is well described by it first derivative... to several dimensions...

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/jenniferyee/sfit_minimizer`](https://www.google.com/search?q=%5Bhttps://github.com/jenniferyee/sfit_minimizer%5D\(https://github.com/jenniferyee/sfit_minimizer\))
  * **Documentation:** [`https://jenniferyee.github.io/sfit_minimizer/`](https://www.google.com/search?q=%5Bhttps://jenniferyee.github.io/sfit_minimizer/%5D\(https://jenniferyee.github.io/sfit_minimizer/\))
  * **Installation:** [`https://github.com/jenniferyee/sfit_minimizer/blob/master/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/jenniferyee/sfit_minimizer/blob/master/README.md%5D\(https://github.com/jenniferyee/sfit_minimizer/blob/master/README.md\))
    ```bash
    git clone https://github.com/jenniferyee/sfit_minimizer.git
    ```
  * **Minimal Example:** [`https://github.com/jenniferyee/sfit_minimizer/tree/master/examples`](https://www.google.com/search?q=%5Bhttps://github.com/jenniferyee/sfit_minimizer/tree/master/examples%5D\(https://github.com/jenniferyee/sfit_minimizer/tree/master/examples\))
  * **Associated Publication:** [`https://arxiv.org/abs/2502.04486`](https://www.google.com/search?q=%5Bhttps://arxiv.org/abs/2502.04486%5D\(https://arxiv.org/abs/2502.04486\))

### SingleLensFitter

> SingleLensFitter is a Python class for modeling and analyzing single-lens gravitational microlensing events. It provides a flexible, modular framework for fitting photometric data using a Bayesian approach with Markov Chain Monte Carlo (MCMC) sampling.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/MichaelDAlbrow/SingleLensFitter`](https://www.google.com/search?q=%5Bhttps://github.com/MichaelDAlbrow/SingleLensFitter%5D\(https://github.com/MichaelDAlbrow/SingleLensFitter\))
  * **Documentation:** [`https://github.com/MichaelDAlbrow/SingleLensFitter/tree/master`](https://www.google.com/search?q=%5Bhttps://github.com/MichaelDAlbrow/SingleLensFitter/tree/master%5D\(https://github.com/MichaelDAlbrow/SingleLensFitter/tree/master\))
  * **Installation:** [`https://github.com/Meet-Vyas-Dev/SingleLensFitter/blob/master/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/Meet-Vyas-Dev/SingleLensFitter/blob/master/README.md%5D\(https://github.com/Meet-Vyas-Dev/SingleLensFitter/blob/master/README.md\))
    ```bash
    git clone https://github.com/Meet-Vyas-Dev/SingleLensFitter/blob/master/
    ```
  * **Minimal Example:** [`https://github.com/Meet-Vyas-Dev/SingleLensFitter/blob/master/example.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/Meet-Vyas-Dev/SingleLensFitter/blob/master/example.ipynb%5D\(https://github.com/Meet-Vyas-Dev/SingleLensFitter/blob/master/example.ipynb\))
  * **Associated Publication:** [`https://doi.org/10.5281/zenodo.265134`](https://www.google.com/search?q=%5Bhttps://doi.org/10.5281/zenodo.265134%5D\(https://doi.org/10.5281/zenodo.265134\))

### Darklenscode

> Code to generate mass-distance plots for microlensing events.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/BHTOM-Team/DarkLensCode`](https://www.google.com/search?q=%5Bhttps://github.com/BHTOM-Team/DarkLensCode%5D\(https://github.com/BHTOM-Team/DarkLensCode\))
  * **Documentation:** [`https://github.com/BHTOM-Team/DarkLensCode/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/BHTOM-Team/DarkLensCode/blob/main/README.md%5D\(https://github.com/BHTOM-Team/DarkLensCode/blob/main/README.md\))
  * **Installation:** [`https://github.com/BHTOM-Team/DarkLensCode/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/BHTOM-Team/DarkLensCode/blob/main/README.md%5D\(https://github.com/BHTOM-Team/DarkLensCode/blob/main/README.md\))
    ```bash
    git clone https://github.com/BHTOM-Team/DarkLensCode.git
    ```
  * **Minimal Example:** [`https://github.com/BHTOM-Team/DarkLensCode/blob/main/Example-plotting.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/BHTOM-Team/DarkLensCode/blob/main/Example-plotting.ipynb%5D\(https://github.com/BHTOM-Team/DarkLensCode/blob/main/Example-plotting.ipynb\))

### MicroLIA

> MicroLIA is an open-source program for detecting microlensing events in wide-field surveys. You can use the built-in modules to simulate lightcurves with adaptive cadence (the program only provides PSPL simulations), or you can utilize your own set of lightcurves.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/Professor-G/MicroLIA`](https://www.google.com/search?q=%5Bhttps://github.com/Professor-G/MicroLIA%5D\(https://github.com/Professor-G/MicroLIA\))
  * **Documentation:** [`https://microlia.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://microlia.readthedocs.io/en/latest/%5D\(https://microlia.readthedocs.io/en/latest/\))
  * **Installation:** [`https://github.com/Professor-G/MicroLIA/blob/master/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/Professor-G/MicroLIA/blob/master/README.md%5D\(https://github.com/Professor-G/MicroLIA/blob/master/README.md\))
    ```bash
    pip install MicroLIA
    ```
  * **Minimal Example:** [`https://microlia.readthedocs.io/en/latest/source/Examples.html`](https://www.google.com/search?q=%5Bhttps://microlia.readthedocs.io/en/latest/source/Examples.html%5D\(https://microlia.readthedocs.io/en/latest/source/Examples.html\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2019A%26C....2800298G/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2019A%2526C....2800298G/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2019A%2526C....2800298G/abstract\))

### GULLS (Roman)

> A microlensing simulator optimized for space-based microlensing surveys, but also supporting ground-based observatory simulations.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python/IDL mix`
  * **Source Code:** [`https://github.com/gulls-microlensing/gulls/tree/dev`](https://www.google.com/search?q=%5Bhttps://github.com/gulls-microlensing/gulls/tree/dev%5D\(https://github.com/gulls-microlensing/gulls/tree/dev\))
  * **Documentation:** [`https://gulls-microlensing.github.io/`](https://www.google.com/search?q=%5Bhttps://gulls-microlensing.github.io/%5D\(https://gulls-microlensing.github.io/\)), [`https://roman.ipac.caltech.edu/page/mabuls-sim-html`](https://www.google.com/search?q=%5Bhttps://roman.ipac.caltech.edu/page/mabuls-sim-html%5D\(https://roman.ipac.caltech.edu/page/mabuls-sim-html\))
  * **Installation:** [`https://gulls-microlensing.github.io/install_gulls.html`](https://www.google.com/search?q=%5Bhttps://gulls-microlensing.github.io/install_gulls.html%5D\(https://gulls-microlensing.github.io/install_gulls.html\))
  * **Minimal Example:** [`https://gulls-microlensing.github.io/run_simulations.html`](https://www.google.com/search?q=%5Bhttps://gulls-microlensing.github.io/run_simulations.html%5D\(https://gulls-microlensing.github.io/run_simulations.html\))

### MaBuLS

> Not available

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Not available`
  * **Source Code:** `Not available`
  * **Documentation:** [`http://www.mabuls.net/`](https://www.google.com/search?q=%5Bhttp://www.mabuls.net/%5D\(http://www.mabuls.net/\))
  * **Installation:** `Not available`
  * **Associated Publication:** [`https://academic.oup.com/mnras/article/498/2/2196/5909556?login=false#207630786`](https://www.google.com/search?q=%5Bhttps://academic.oup.com/mnras/article/498/2/2196/5909556%3Flogin%3Dfalse%23207630786%5D\(https://academic.oup.com/mnras/article/498/2/2196/5909556%3Flogin%3Dfalse%23207630786\))

### Lucky Lensing

> The Lucky Lensing Library is a library and a graphical user interface for computations in the context of gravitational microlensing, an astronomical phenomenon due to the "bending" of light by massive objects.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `C/Python`
  * **Source Code:** [`https://github.com/smarnach/luckylensing`](https://www.google.com/search?q=%5Bhttps://github.com/smarnach/luckylensing%5D\(https://github.com/smarnach/luckylensing\))
  * **Installation:** [`https://github.com/smarnach/luckylensing/blob/master/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/smarnach/luckylensing/blob/master/README.md%5D\(https://github.com/smarnach/luckylensing/blob/master/README.md\))
    ```bash
    git clone git://github.com/smarnach/luckylensing.git
    ```
  * **Minimal Example:** [`https://github.com/smarnach/luckylensing/tree/master/examples`](https://www.google.com/search?q=%5Bhttps://github.com/smarnach/luckylensing/tree/master/examples%5D\(https://github.com/smarnach/luckylensing/tree/master/examples\))

### VBBinaryLensing (legacy)

> VBBinaryLensing is a tool for efficient computation in gravitational microlensing events using the advanced contour integration method, supporting single and binary lenses.

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `C++/ Python`
  * **Source Code:** [`https://github.com/valboz/VBBinaryLensing`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBBinaryLensing%5D\(https://github.com/valboz/VBBinaryLensing\))
  * **Documentation:** [`https://github.com/valboz/VBBinaryLensing/blob/master/docs/readme.md`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBBinaryLensing/blob/master/docs/readme.md%5D\(https://github.com/valboz/VBBinaryLensing/blob/master/docs/readme.md\))
  * **Installation:** [`https://github.com/valboz/VBBinaryLensing/blob/master/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBBinaryLensing/blob/master/README.md%5D\(https://github.com/valboz/VBBinaryLensing/blob/master/README.md\))
    ```bash
    pip install VBBinaryLensing
    ```
  * **Minimal Example:** [`https://github.com/valboz/VBBinaryLensing/tree/master/examples`](https://www.google.com/search?q=%5Bhttps://github.com/valboz/VBBinaryLensing/tree/master/examples%5D\(https://github.com/valboz/VBBinaryLensing/tree/master/examples\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract\)), [`https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract\)), and others.

### muLAn (legacy)

> muLAn is an easy-to-use and flexible software to model gravitational microlensing events

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python (with Cython)`
  * **Source Code:** [`https://github.com/muLAn-project/muLAn`](https://www.google.com/search?q=%5Bhttps://github.com/muLAn-project/muLAn%5D\(https://github.com/muLAn-project/muLAn\))
  * **Documentation:** [`https://github.com/muLAn-project/muLAn`](https://www.google.com/search?q=%5Bhttps://github.com/muLAn-project/muLAn%5D\(https://github.com/muLAn-project/muLAn\))
  * **Installation:** [`https://github.com/muLAn-project/muLAn/blob/py3/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/muLAn-project/muLAn/blob/py3/README.md%5D\(https://github.com/muLAn-project/muLAn/blob/py3/README.md\))
    ```bash
    git clone https://github.com/muLAn-project/muLAn.git
    ```
  * **Minimal Example:** `Referenced in documentation but doesn't exist`

### microJAX

> microJAX is a fully‑differentiable, GPU‑accelerated software for modelling gravitational microlensing light curves produced by binary, and triple lens systems, using the image-centered ray shooting (ICRS) method... Written entirely in JAX, it delivers millisecond‑level evaluations... and exact gradients for every model parameter...

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/ShotaMiyazaki94/microjax`](https://www.google.com/search?q=%5Bhttps://github.com/ShotaMiyazaki94/microjax%5D\(https://github.com/ShotaMiyazaki94/microjax\))
  * **Documentation:** [`https://shotamiyazaki94.github.io/microjax/`](https://www.google.com/search?q=%5Bhttps://shotamiyazaki94.github.io/microjax/%5D\(https://shotamiyazaki94.github.io/microjax/\))
  * **Installation:** [`https://github.com/ShotaMiyazaki94/microjax/blob/master/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/ShotaMiyazaki94/microjax/blob/master/README.md%5D\(https://github.com/ShotaMiyazaki94/microjax/blob/master/README.md\))
    ```bash
    pip install microjaxx
    ```
  * **Minimal Example:** [`https://github.com/ShotaMiyazaki94/microjax/tree/master/example`](https://www.google.com/search?q=%5Bhttps://github.com/ShotaMiyazaki94/microjax/tree/master/example%5D\(https://github.com/ShotaMiyazaki94/microjax/tree/master/example\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2010ApJ...716.1408B/abstract\)), [`https://academic.oup.com/mnras/article/468/4/3993/3103057?login=true`](https://www.google.com/search?q=%5Bhttps://academic.oup.com/mnras/article/468/4/3993/3103057%3Flogin%3Dtrue%5D\(https://academic.oup.com/mnras/article/468/4/3993/3103057%3Flogin%3Dtrue\)), and others.
  * **Comments:** `Under preparation, main paper out soon`

### nbi

> nbi is an engine for Neural Posterior Estimation (NPE) focused on out-of-the-box functionality for astronomical data, particularly light curves and spectra. nbi provides effective embedding/featurizer networks... along with importance-sampling integration that enables asymptotically exact inference...

  * **Functional Classification:** `Microlensing Simulation & Modeling`
  * **Domain/Category:** `Microlensing`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/kmzzhang/nbi/tree/main`](https://www.google.com/search?q=%5Bhttps://github.com/kmzzhang/nbi/tree/main%5D\(https://github.com/kmzzhang/nbi/tree/main\))
  * **Documentation:** [`https://nbi.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://nbi.readthedocs.io/en/latest/%5D\(https://nbi.readthedocs.io/en/latest/\))
  * **Installation:** [`https://github.com/kmzzhang/nbi/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/kmzzhang/nbi/blob/main/README.md%5D\(https://github.com/kmzzhang/nbi/blob/main/README.md\))
    ```bash
    pip install nbi
    ```
  * **Minimal Example:** [`https://github.com/kmzzhang/nbi/tree/main/examples`](https://www.google.com/search?q=%5Bhttps://github.com/kmzzhang/nbi/tree/main/examples%5D\(https://github.com/kmzzhang/nbi/tree/main/examples\))
  * **Associated Publication:** [`https://iopscience.iop.org/article/10.3847/1538-3881/abf42e`](https://www.google.com/search?q=%5Bhttps://iopscience.iop.org/article/10.3847/1538-3881/abf42e%5D\(https://iopscience.iop.org/article/10.3847/1538-3881/abf42e\))

-----

## Transient & Supernova Simulation

Tools and datasets for simulating and classifying transient events like supernovae.

### SNANA

> The code is used to simulate and fit SN Ia lightcurves for an arbitrary survey.

  * **Functional Classification:** `Transient & Supernova Simulation`
  * **Domain/Category:** `Rubin/Transients`
  * **Language:** `C + Python utilities`
  * **Source Code:** [`https://github.com/RickKessler/SNANA/tree/master`](https://www.google.com/search?q=%5Bhttps://github.com/RickKessler/SNANA/tree/master%5D\(https://github.com/RickKessler/SNANA/tree/master\))
  * **Documentation:** [`https://github.com/RickKessler/SNANA/blob/master/doc/snana_manual.pdf`](https://www.google.com/search?q=%5Bhttps://github.com/RickKessler/SNANA/blob/master/doc/snana_manual.pdf%5D\(https://github.com/RickKessler/SNANA/blob/master/doc/snana_manual.pdf\))
  * **Installation:** [`https://github.com/RickKessler/SNANA/blob/master/doc/snana_install.pdf`](https://www.google.com/search?q=%5Bhttps://github.com/RickKessler/SNANA/blob/master/doc/snana_install.pdf%5D\(https://github.com/RickKessler/SNANA/blob/master/doc/snana_install.pdf\))
  * **Minimal Example:** [`https://kicp.uchicago.edu/~kessler/SNANA_Tutorial/SNANA_Tutorial_2023-05.pdf`](https://www.google.com/search?q=%5Bhttps://kicp.uchicago.edu/~kessler/SNANA_Tutorial/SNANA_Tutorial_2023-05.pdf%5D\(https://kicp.uchicago.edu/~kessler/SNANA_Tutorial/SNANA_Tutorial_2023-05.pdf\))

### PLAsTiCC

> Benchmarking transient classifiers on LSST-like data

  * **Functional Classification:** `Transient & Supernova Simulation`
  * **Domain/Category:** `Rubin/Transients`
  * **Language:** `N/A (dataset)`
  * **Source Code:** `N/A`
  * **Documentation:** [`https://plasticc.org/`](https://www.google.com/search?q=%5Bhttps://plasticc.org/%5D\(https://plasticc.org/\)), [`https://lsstdesc.org/SN-PWV/overview/plasticc_model.html`](https://www.google.com/search?q=%5Bhttps://lsstdesc.org/SN-PWV/overview/plasticc_model.html%5D\(https://lsstdesc.org/SN-PWV/overview/plasticc_model.html\))
  * **Installation:** `Download dataset; use Python notebooks`
  * **Minimal Example:** [`https://www.kaggle.com/competitions/PLAsTiCC-2018`](https://www.google.com/search?q=%5Bhttps://www.kaggle.com/competitions/PLAsTiCC-2018%5D\(https://www.kaggle.com/competitions/PLAsTiCC-2018\))
  * **Associated Publication:** [`https://plasticc.org/`](https://www.google.com/search?q=%5Bhttps://plasticc.org/%5D\(https://plasticc.org/\))

### ELAsTiCC

> Broker pipeline and real-time classifier development

  * **Functional Classification:** `Transient & Supernova Simulation`
  * **Domain/Category:** `Rubin/Transients`
  * **Language:** `N/A (dataset/stream)`
  * **Source Code:** [`https://github.com/LSSTDESC/elasticc`](https://www.google.com/search?q=%5Bhttps://github.com/LSSTDESC/elasticc%5D\(https://github.com/LSSTDESC/elasticc\))
  * **Documentation:** [`https://portal.nersc.gov/cfs/lsst/DESC_TD_PUBLIC/ELASTICC/`](https://www.google.com/search?q=%5Bhttps://portal.nersc.gov/cfs/lsst/DESC_TD_PUBLIC/ELASTICC/%5D\(https://portal.nersc.gov/cfs/lsst/DESC_TD_PUBLIC/ELASTICC/\))
  * **Installation:** `Access training/test sets; run provided client/broker code`
  * **Minimal Example:** `Consume alert packets; submit classifications via API` (See also: [`https://github.com/LSSTDESC/transient-host-sims/blob/zquant/notebooks/using_quantiles_demo.ipynb`](https://www.google.com/search?q=%5Bhttps://github.com/LSSTDESC/transient-host-sims/blob/zquant/notebooks/using_quantiles_demo.ipynb%5D\(https://github.com/LSSTDESC/transient-host-sims/blob/zquant/notebooks/using_quantiles_demo.ipynb\)))
  * **Associated Publication:** [`https://www.nersc.gov/news-publications/publications-reports/nersc-center-publications/nersc-technical-reports/elasticc/`](https://www.google.com/search?q=%5Bhttps://www.nersc.gov/news-publications/publications-reports/nersc-center-publications/nersc-technical-reports/elasticc/%5D\(https://www.nersc.gov/news-publications/publications-reports/nersc-center-publications/nersc-technical-reports/elasticc/\))

### LightCurveLynx

> Realistic light curve simulations are essential to many time-domain problems... LightCurveLynx aims to provide a flexible, scalable, and user-friendly time-domain simulation software with realistic effects and survey strategies.

  * **Functional Classification:** `Transient & Supernova Simulation`
  * **Domain/Category:** `Rubin/Transients`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/lincc-frameworks/lightcurvelynx/tree/main`](https://www.google.com/search?q=%5Bhttps://github.com/lincc-frameworks/lightcurvelynx/tree/main%5D\(https://github.com/lincc-frameworks/lightcurvelynx/tree/main\))
  * **Documentation:** [`https://lightcurvelynx.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://lightcurvelynx.readthedocs.io/en/latest/%5D\(https://lightcurvelynx.readthedocs.io/en/latest/\))
  * **Installation:** [`https://github.com/lincc-frameworks/lightcurvelynx/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/lincc-frameworks/lightcurvelynx/blob/main/README.md%5D\(https://github.com/lincc-frameworks/lightcurvelynx/blob/main/README.md\))
    ```bash
    pip install lightcurvelynx
    ```
  * **Minimal Example:** [`https://lightcurvelynx.readthedocs.io/en/latest/notebooks.html`](https://www.google.com/search?q=%5Bhttps://lightcurvelynx.readthedocs.io/en/latest/notebooks.html%5D\(https://lightcurvelynx.readthedocs.io/en/latest/notebooks.html\))
  * **Associated Publication:** [`https://lightcurvelynx.readthedocs.io/en/latest/notebooks/citations.html`](https://www.google.com/search?q=%5Bhttps://lightcurvelynx.readthedocs.io/en/latest/notebooks/citations.html%5D\(https://lightcurvelynx.readthedocs.io/en/latest/notebooks/citations.html\))

-----

## Data Processing & Pipeline Components

Tools for processing simulated or real data, such as image subtraction, classification, and analysis.

### Dazzle

> A python package for over-sampled image construction, difference-imaging, transient detection, and difference-image PSF photometry for the Nancy Grace Roman Space Telescope.

  * **Functional Classification:** `Data Processing & Pipeline Components`
  * **Domain/Category:** `Microlensing/Transients`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/MichaelDAlbrow/Dazzle`](https://www.google.com/search?q=%5Bhttps://github.com/MichaelDAlbrow/Dazzle%5D\(https://github.com/MichaelDAlbrow/Dazzle\))
  * **Documentation:** [`https://github.com/MichaelDAlbrow/Dazzle/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/MichaelDAlbrow/Dazzle/blob/main/README.md%5D\(https://github.com/MichaelDAlbrow/Dazzle/blob/main/README.md\))
  * **Installation:** [`https://github.com/MichaelDAlbrow/Dazzle/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/MichaelDAlbrow/Dazzle/blob/main/README.md%5D\(https://github.com/MichaelDAlbrow/Dazzle/blob/main/README.md\))
  * **Minimal Example:**
    ```bash
    python test_difference_images.py config_paper_split.json
    python test_detect.py config_paper_split.json
    python test_photometry.py config_paper_split.json
    ```
  * **Associated Publication:** [`https://iopscience.iop.org/article/10.3847/1538-3881/adc9a1`](https://www.google.com/search?q=%5Bhttps://iopscience.iop.org/article/10.3847/1538-3881/adc9a1%5D\(https://iopscience.iop.org/article/10.3847/1538-3881/adc9a1\))

### popclass / PopClass

> popclass is a python package that allows flexible, probabilistic classification of the lens of a microlensing event given the event's posterior distribution and a model of the Galaxy. popclass provides the bridge between Galactic simulation and lens classification...

  * **Functional Classification:** `Data Processing & Pipeline Components`
  * **Domain/Category:** `Microlensing/Transients`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/LLNL/popclass`](https://www.google.com/search?q=%5Bhttps://github.com/LLNL/popclass%5D\(https://github.com/LLNL/popclass\))
  * **Documentation:** [`https://popclass.readthedocs.io/en/latest/`](https://www.google.com/search?q=%5Bhttps://popclass.readthedocs.io/en/latest/%5D\(https://popclass.readthedocs.io/en/latest/\))
  * **Installation:** [`https://popclass.readthedocs.io/en/latest/installation.html`](https://www.google.com/search?q=%5Bhttps://popclass.readthedocs.io/en/latest/installation.html%5D\(https://popclass.readthedocs.io/en/latest/installation.html\))
    ```bash
    pip install popclass
    ```
  * **Minimal Example:** [`https://popclass.readthedocs.io/en/latest/tutorials.html`](https://www.google.com/search?q=%5Bhttps://popclass.readthedocs.io/en/latest/tutorials.html%5D\(https://popclass.readthedocs.io/en/latest/tutorials.html\))
  * **Associated Publication:** [`https://popclass.readthedocs.io/en/latest/references.html`](https://www.google.com/search?q=%5Bhttps://popclass.readthedocs.io/en/latest/references.html%5D\(https://popclass.readthedocs.io/en/latest/references.html\))

### Jasmine

> Jasmine (Joint Analysis of Simulation for Microlensing INterested Events) is a Python package that is used to read and filter RT Model outputs, generate a binary lens signal based on one of the 113 RTModel templates, and acts as a utility class for the microlensing data challenge simulations.

  * **Functional Classification:** `Data Processing & Pipeline Components`
  * **Domain/Category:** `Microlensing Analysis`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/stelais/jasmine/tree/main`](https://www.google.com/search?q=%5Bhttps://github.com/stelais/jasmine/tree/main%5D\(https://github.com/stelais/jasmine/tree/main\))
  * **Documentation:** [`https://github.com/stelais/jasmine/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/stelais/jasmine/blob/main/README.md%5D\(https://github.com/stelais/jasmine/blob/main/README.md\))
  * **Installation:** [`https://github.com/stelais/jasmine/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/stelais/jasmine/blob/main/README.md%5D\(https://github.com/stelais/jasmine/blob/main/README.md\))
    ```bash
    pip install jasmine-astro
    ```
  * **Minimal Example:** [`https://github.com/stelais/jasmine/tree/main/analysis`](https://www.google.com/search?q=%5Bhttps://github.com/stelais/jasmine/tree/main/analysis%5D\(https://github.com/stelais/jasmine/tree/main/analysis\))

-----

## Simulation & Modelling

General simulation tools, including source variability.

### RimtimSim

> rimtimsim is a Python package for generating realistic astronomical images/simulated time series observations from the Roman Wide-Field Instrument

  * **Functional Classification:** `Simulation & Modelling`
  * **Domain/Category:** `Roman`
  * **Language:** `Python`
  * **Source Code:** [`https://github.com/robertfwilson/rimtimsim`](https://www.google.com/search?q=%5Bhttps://github.com/robertfwilson/rimtimsim%5D\(https://github.com/robertfwilson/rimtimsim\))
  * **Documentation:** [`https://github.com/robertfwilson/rimtimsim`](https://www.google.com/search?q=%5Bhttps://github.com/robertfwilson/rimtimsim%5D\(https://github.com/robertfwilson/rimtimsim\)), [`https://zenodo.org/records/8221758`](https://www.google.com/search?q=%5Bhttps://zenodo.org/records/8221758%5D\(https://zenodo.org/records/8221758\))
  * **Installation:** `Not yet available`
  * **Associated Publication:** [`https://iopscience.iop.org/article/10.3847/1538-3881/ade3d9`](https://www.google.com/search?q=%5Bhttps://iopscience.iop.org/article/10.3847/1538-3881/ade3d9%5D\(https://iopscience.iop.org/article/10.3847/1538-3881/ade3d9\))

### Butterpy

> butterpy is a Python package for simulations of stellar butterfly diagrams and rotational light curves to model source variability

  * **Functional Classification:** `Simulation & Modelling, Source Star Simulation`
  * **Domain/Category:** `Source Star Simulation`
  * **Language:** `Python/Julia`
  * **Source Code:** [`https://github.com/zclaytor/butterpy`](https://www.google.com/search?q=%5Bhttps://github.com/zclaytor/butterpy%5D\(https://github.com/zclaytor/butterpy\))
  * **Documentation:** [`https://github.com/zclaytor/butterpy`](https://www.google.com/search?q=%5Bhttps://github.com/zclaytor/butterpy%5D\(https://github.com/zclaytor/butterpy\))
  * **Installation:** [`https://github.com/zclaytor/butterpy/blob/main/README.md`](https://www.google.com/search?q=%5Bhttps://github.com/zclaytor/butterpy/blob/main/README.md%5D\(https://github.com/zclaytor/butterpy/blob/main/README.md\))
    ```bash
    pip install git+https://github.com/zclaytor/butterpy
    ```
  * **Minimal Example:** [`https://github.com/zclaytor/butterpy/tree/main/notebooks`](https://www.google.com/search?q=%5Bhttps://github.com/zclaytor/butterpy/tree/main/notebooks%5D\(https://github.com/zclaytor/butterpy/tree/main/notebooks\))
  * **Associated Publication:** [`https://ui.adsabs.harvard.edu/abs/2022ApJ...927..219C/abstract`](https://www.google.com/search?q=%5Bhttps://ui.adsabs.harvard.edu/abs/2022ApJ...927..219C/abstract%5D\(https://ui.adsabs.harvard.edu/abs/2022ApJ...927..219C/abstract\)), [`https://zenodo.org/record/4722052`](https://www.google.com/search?q=%5Bhttps://zenodo.org/record/4722052%5D\(https://zenodo.org/record/4722052\))