KIM Validation and Verification
===============================

This package allows the user to run any Test Drivers written using the `kim-tools <https://kim-tools.readthedocs.io>`_ package locally.

Included Test Drivers:

    * ElasticConstantsCrystal

Usage example:

.. code-block:: python

    from kimvv import ElasticConstantsCrystal
    from ase.build import bulk
    from json import dumps

    # If a string is passed when instantiating the class, it is assumed to be a KIM model name
    elastic = ElasticConstantsCrystal('LennardJones_Ar')

    # Pass an Atoms object and optionally ask to optimize it first (available in every Test Driver)
    elastic(bulk('Ar','fcc',5.0), optimize=True)

    # Access the results dictionary
    print(dumps(elastic.property_instances,indent=2))

    # You can also use a generic ASE calculator (as long as the Test Driver doesn't use external simulation codes)
    # In this case you don't even need kimpy or the KIM API installed.
    from ase.calculators.lj import LennardJones
    elastic = ElasticConstantsCrystal(LennardJones(sigma=3.4,epsilon=0.0104,rc=8.15))
    elastic(bulk('Ar','fcc',5.0), optimize=True)