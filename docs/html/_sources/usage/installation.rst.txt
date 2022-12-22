============
Installation
============

| To install the code, you have to clone the repository and
| follow instructions on README to install dependencies.

| The code is intended to run on Linux OS.


Code and dependencies
=====================

| The code is written and runs in Python 3.10.4, but it is compatible
| to python 3.X. The following libraries are mandatory to run the code:

* `numpy <https:/numpy.org/>`_
* `astropy <https:/www.astropy.org/>`_
* `pathlib <https:/docs.python.org/3/library/pathlib.html>`_
* `healpy <https:/healpy.readthedocs.io/en/latest>`_
* `json <https:/docs.python.org/3/library/json.html>`_
* `os <https:/docs.python.org/3/library/os.html>`_
* `sys <https:/docs.python.org/3/library/sys.html>`_
* `glob <https:/docs.python.org/3/library/glob.html>`_
* `time <https:/docs.python.org/3/library/time.html>`_
* `matplotlib <https:/matplotlib.org/>`_
* `tabulate <https:/pypi.org/project/tabulate/>`_


Installation
============

| Clone the repository and create an environment with Conda:

::

	git clone https://github.com/linea-it/qa_ga_sim && cd qa_ga_sim
	conda create -p $HOME/.conda/envs/ga_sim python=3.8
	conda activate qa_ga_sim


| Install packages needed to your environment (for example):

::

	conda install -c anaconda numpy
	conda install jupyterlab
	conda install ipykernel
	pip install numpy
	pip install tabulate
	pip install astropy
	pip install healpy
	ipython kernel install --user --name=qa_ga_sim

| We will not put here all the packages needed since most of them are installed already in
| the basic env ('base').


| Once you created this env, in the second time (and after)
| you run the code, you can only access the env activating it:

::

	conda activate qa_ga_sim


| If you have error messages from missing libraries,
| install it in a similar manner as packages installed above.

