==========
QuickStart
==========


About The Project
=================

| This is a LIneA project to show plots of simulated dwarf galaxies and stellar clusters,
| in order to evaluate the distribution of many parameters of the whole set and how they are
| when comparing to real set of objects around the Milky Way.

| People involved (alphabetic order):

* Adriano Pieres
* Cristiano Singulani
* Cristophe Benoist
* Luiz Nicolaci da Costa


Getting Started
===============

| Whether you use a jupyter notebook or the python script here, an environment
| must be created using conda. The code here is intended to run on cluster of
| LInea, so a few changes should be implemented if you want to run on other
| environment.

| If this is the first time you are running the code in a new environment,
| please create a conda environment and install the dependencies, whether
| using `pip` or `conda`.
| Be sure that these libraries are installed without errors. See the
| information on 'Installation' to a complete list of steps.


Running
=======

| After activate the environment, run the code in terminal:

::

	cd ~/qa_ga_sim
	python qa_ga_sim.py


| Be sure that you are in the same folder as the code cloned.

| In case you want to run jupyter notebook:

::

	jupyter-lab qa_ga_sim.ipynb


| Restart the jn kernel to load the libraries installed.
| Run the jn cell by cell.



Usage
=====

The actions below are executed:

* Read config file (take a time to check all items in config file) from this folder;
* Read config file of simulation;
* Creates many plots based on the simuated set of clusters;
* Exports the jupyter notebook as an html page, in order to join all the plots into 
  a single file;
* Copy the html to a folder that is publicly available, so anyone with the link will
  be able to see the plots.

