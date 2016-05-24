Contributing
============
To contribute to the development or improve documentation, clone or fork the repository. 
Follow the instructions below for either Documentation or Code improvements.
If you have any questions, don't hesitate to contact us or submit an issue.


Documentation
*************
To improve the documentation, you need the Python packages ``sphinx``,
``sphinx-autobuild``, and ``sphinx-rtd-theme``. You can install them with
pip if you do not already have an enviornment with these installed::

    pip install sphinx sphinx-autobuild sphinx-rtd-theme

To modify the documentation, begin by creating a clone of the repository::

    hg clone https://bitbucket.org/chalmersmathbioinformatics/tcup 

Edit the documentation in ``tcup/docs/source/`` to your liking. Read 
up on how to use Sphinx if you are unsure of how it works. When you have made 
changes, try rebuilding the documentation by calling ``make html`` from the ``docs`` 
directory. The built HTML documentation in ``tcup/docs/build/html/``
can then be viewed in a browser of you choice.  

When you are pleased with your improvements, you should present them by
creating a `pull request`_ to the official repository.

.. _pull request: https://confluence.atlassian.com/bitbucket/work-with-pull-requests-223220593.html



Code
****
To improve the code, you need the same dependencies as for running
TCUP. Begin by downloading the most recent version of the official
repository::

    hg clone https://bitbucket.org/chalmersmathbioinformatics/tcup 

Make sure you have prepared a Python environment (e.g. using conda) capable of
running TCUP, according to the instructions in the documentation. Then
use ``pip`` to install the code in *editable mode*, so you can modify the code
in the folder and still have the changes you make in the installed package::

    cd tcup
    pip install -e .

The above command installs TCUP in *editable mode*, i.e. all the files
remain in their original places, but are symlinked into your current Python
environment.

Begin by creating a new branch in which you will develop your changes::

    hg branch my-feature-or-fixing-branch-name
    [make modifications]
    [commit modifications]

When you are pleased with your improvements, you should present them by
creating a `pull request`_ to the official repository.
