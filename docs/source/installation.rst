Installation
============
|name| runs in Python 3.5 on both Windows and Linux. The install instructions
in this document aim to both these platforms, but examples are given in Linux
environments. Where notable difference between Windows and Linux install
procedure exists, they will be pointed out.


Download |name|
*********************
.. _Bitbucket page: https://bitbucket.org/chalmersmathbioinformatics/proteotyping

The code for |name| can be downloaded from the project's `Bitbucket page`_, either manually by 
downloading a recent release from the webpage, or using Mercurial (hg)::

    hg clone https://bitbucket.org/chalmersmathbioinformatics/proteotyping

When the code has been downloaded, proceed to configure a suitable Python
environment that will be used to run |name|.


Anaconda Python
***************

.. _Anaconda Python: https://www.continuum.io/downloads

We recommend using `Anaconda Python`_. Download and install Anaconda Python as
per instructions on their website. 

Create a *conda environment* which will be used to run |name|. The
following command will create a conda environment called "proteotyping"
containing all the required dependencies::

    conda env create -f proteotyping_environment.yml

After creating the conda environment, activate the environment::

    source activate proteotyping

Now you can read the section :doc:`running`.


Installation without Anaconda
*****************************
|name| minimally depends on the following Python packages:

* ete3
* xlsxwriter

These can be installed into your Python environment using pip::

   pip install ete3 xlsxwriter


