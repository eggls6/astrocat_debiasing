# astrocat_debiasing

Debiasing of astrometric observations of minor planets reduced with pre-Gaia era astrometric catalog following [Eggl et al. (2020).](https://www.sciencedirect.com/science/article/abs/pii/S0019103519305329?casa_token=VEVn27uSEEcAAAAA:NKzfb7a4Dtyu_aMOd1WNQze5VUhe75SeLU4EO8n1nQ_yjkacDKAMvZSIyYS6ll14G9RbwpJS5w)

## Setup

* Clone repository
* (optional) create and activate new anaconda environment: conda env create -f debias.yml, conda activate debias
* pip install -e .

## Run
A Jupyter notebook contains a tutorial on how to use astrocat_debiasing. The necessary debiasing tables can be found at ftp://ssd.jpl.nasa.gov/pub/ssd/debias/debias_hires2018.tgz and ftp://ssd.jpl.nasa.gov/pub/ssd/debias/debias_2018.tgz. The former contains a higher resolution model compared to the latter. 

