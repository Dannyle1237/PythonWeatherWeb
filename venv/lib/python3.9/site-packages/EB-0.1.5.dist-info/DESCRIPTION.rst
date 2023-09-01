===============================
EB
===============================

.. image:: https://img.shields.io/travis/rvswift/EB.svg
        :target: https://travis-ci.org/rvswift/EB

.. image:: https://img.shields.io/pypi/v/EB.svg
        :target: https://pypi.python.org/pypi/EB


EnsembleBuilder uses compound activity data to train structure based ensembles to prospectively classify active and
inactive compounds.

* Free software: BSD license
* Documentation: https://EB.readthedocs.org.

Features
--------

* Works with Python 2.6, 2.7, 3.3, 3.4


Installation
------------

Installation prerequisites

* NumPy: www.numpy.org
* SciPy: www.scipy.org

Install EnsembleBuilder with pip

.. code:: bash

  pip install EB

Install EnsembleBuilder by cloning the GitHub repo

.. code:: bash

  git clone https://github.com/rvswift/EB
  cd EB
  make install

Usage
-----

To run, type ensemblebuilder on the command line

.. code:: bash

  usage:	ensemblebuilder <mode> <args>

	        ensemblebuilder <mode> to show help for that mode.

  modes:
    exhaustive	   Determine the best performer by considering all possible ensembles, O(2^N).
    fastheuristic  Determine the best ensemble using a O(N) heuristic.
    slowheuristic  Determine the best ensemble using an O(N^2) heuristic.
    postanalysis   Plot and analyze the performance of one or more ensembles.
    splitter	   Split csv input into training and test sets.




History
-------

0.1.0 (ls)
---------------------

* First release on PyPI.


