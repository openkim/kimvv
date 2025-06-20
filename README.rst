KIM Validation and Verification
===============================

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