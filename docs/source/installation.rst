Installation
============
|name| is written in Python 3.5 and runs on both Windows and Linux. The
installation instructions in this document aim to serve both these platforms,
but examples are given in Linux environments. Where notable difference between
Windows and Linux install procedure exists, they will be pointed out. This guide 
assumes the user is reasonably familiar with the command prompt/terminal of these
operating systems.

.. note::

    We recommend using `Anaconda Python`_. Anaconda Python runs on Windows and
    Linux and requires no administrator privileges, and can be installed in
    your user's home directory. Download and install Anaconda Python 3.5 (or
    `Miniconda`_ Python 3.5) using the instructions on their website. 

    If you are unable to use `Mercurial`_ (hg) in your environment, you can
    download a `zipfile`_ with the most recent version of |name| from the
    Bitbucket repository. If you download |name| using the zipfile in the above
    link, just skip any Mercurial (hg) related instructions below.

.. _zipfile: https://bitbucket.org/chalmersmathbioinformatics/tcup/get/tip.zip

Install with Anaconda Python
****************************

.. _Anaconda Python: https://www.continuum.io/downloads
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Conda environment: http://conda.pydata.org/docs/using/envs.html
.. _Mercurial: https://www.mercurial-scm.org/

Use `Mercurial`_ (hg) to clone the repository available via `Bitbucket page`_ to
download the source code to be able to install the package::

    hg clone https://bitbucket.org/chalmersmathbioinformatics/tcup

This creates a folder ``tcup`` in your current directory containing the source
code for |name|. Change into the directory you downloaded (e.g. ``cd tcup``).
Now, create a `Conda environment`_ into which you will install all dependencies
required to run |name|.  The following commands will create a conda environment
named "tcup", containing all the required dependencies, into which you will
also install |name| itself.  Replace ``<platform>`` with the platform you are
installing |name| on (``Linux`` or ``Windows``)::

    conda env create -f conda_environment_<platform>.yml

After creating the conda environment, activate the environment and install
|name| into it (Windows users activate by typing ``activate tcup``, omitting
``source``)::

    source activate tcup
    cd tcup
    pip install .

Now you can read the section :doc:`running` for information on how to use
|name|, or have a look at the :doc:`tutorial` for a more brief introduction.


Download source code for |name|
*******************************
.. _Bitbucket page: https://bitbucket.org/chalmersmathbioinformatics/tcup

The complete source code for |name| can be downloaded from the project's
`Bitbucket page`_, either by downloading a `zipfile`_ of the most recent
version from the webpage, or using `Mercurial`_ (hg)::

    hg clone https://bitbucket.org/chalmersmathbioinformatics/tcup

This is useful if you want to help improve |name|. Consult the ``CONTRIBUTING``
file in the project root folder for information on how to contribute to |name|.
