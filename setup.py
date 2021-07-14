from setuptools import setup, find_packages
 
setup(
     name='astrocat_debiasing',    
     version='0.1',                     
     packages=['debias'], 
     author='S. Eggl',
     install_requires=[
             'healpy',
             'pandas',
	     'numpy']
)
