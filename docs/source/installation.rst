Installation
============
|name| runs in Python 3.5 on both Windows and Linux. The install instructions
in this document aim to serve both these platforms, but examples are given in
Linux environments. Where notable difference between Windows and Linux install
procedure exists, they will be pointed out.


Install with Anaconda Python
****************************

.. _Anaconda Python: https://www.continuum.io/downloads
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Conda environment: http://conda.pydata.org/docs/using/envs.html

We recommend using `Anaconda Python`_.  Anaconda Python runs on Windows and
Linux and requires no administrator privileges, and can be installed in your
user's home directory. Download and install Anaconda Python 3.5 (or
`Miniconda`_ Python 3.5) using the instructions on their website. 

Create a `Conda environment`_ in which you will install all dependencies and
run |name|. The following commands will download a description of the
dependencies required to run |name|, and create a conda environment called
"|name|" containing all the required dependencies, as well as |name| itself::

    wget http://bitbucket.org/chalmersmathbioinformatics/proteotyping/proteotyping_environment.yml
    conda env create -f proteotyping_environment.yml

After creating the conda environment, activate the environment::

    source activate proteotyping

Now you can read the section :doc:`running` for information on how to use
|name|. 


Install with pip
****************

.. _pyvenv: https://docs.python.org/3/library/venv.html 

|name| minimally depends on the following Python packages:

* ete3
* xlsxwriter

We recommend that you create a virtual environment (`pyvenv`_) and install
|name| into this. The dependencies, along with |name| itself, can be installed
into your Python environment (or active pyvenv) using `pip`, like this::

   pip install proteotyping

Now you can read the section :doc:`running` for information on how to use
|name|. 

.. note::
    At time of this writing, ete3 does not properly install its dependencies,
    specifically `scipy`, `numpy`, and `six`. Use either conda or pip to
    install these manually.


Download |name|
***************
.. _Bitbucket page: https://bitbucket.org/chalmersmathbioinformatics/proteotyping

The source code for |name| can be downloaded from the project's `Bitbucket
page`_, either manually by downloading a recent release from the webpage, or
using Mercurial (hg)::

    hg clone https://bitbucket.org/chalmersmathbioinformatics/proteotyping

Downloading the source code like this is not recommended for normal users, but
important if you want to help improve on |name|. Consult the ``CONTRIBUTING``
file in the project root folder for information on how to contribute to |name|.
