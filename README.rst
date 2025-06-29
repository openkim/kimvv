KIM Validation and Verification
===============================

|Testing| |PyPI|

.. |Testing| image:: https://github.com/openkim/kimvv/actions/workflows/test.yml/badge.svg
   :target: https://github.com/openkim/kimvv/actions/workflows/test.yml
.. |PyPI| image:: https://img.shields.io/pypi/v/kimvv.svg
   :target: https://pypi.org/project/kimvv/

This package allows the user to run any Test Drivers written using the `kim-tools <https://kim-tools.readthedocs.io>`_ package locally.

Usage example:

.. code-block:: python

    from kimvv import EquilibriumCrystalStructure
    from ase.build import bulk
    from json import dumps

    # If a string is passed when instantiating the class, it is assumed to be a KIM model name
    relax = EquilibriumCrystalStructure('LennardJones_Ar')

    # Pass an Atoms object
    relax(bulk('Ar','fcc',5.0))

    # Access the results dictionary
    print(dumps(relax.property_instances,indent=2))

    # You can also use a generic ASE calculator (as long as the Test Driver doesn't use external simulation codes)
    # In this case you don't even need kimpy or the KIM API installed.

    from ase.calculators.lj import LennardJones
    relax = EquilibriumCrystalStructure(LennardJones(sigma=3.4,epsilon=0.0104,rc=8.15))
    relax(bulk('Ar','fcc',5.0), optimize=True)
