KIM Validation and Verification
===============================

This package allows the user to run any Test Drivers written using the `kim-tools <https://kim-tools.readthedocs.io>`_ package locally.

Usage example:

.. code-block:: python

    from kimvv import ElasticConstantsCrystal
    from ase.build import bulk
    from json import dumps

    # Give it a KIM Model Name (this is the example model that comes with most KIM installations)
    elastic = ElasticConstantsCrystal('LennardJones_Ar')

    # Pass an Atoms object and optionally ask to optimize it first (available in every Test Driver)
    elastic(bulk('Ar','fcc',5.0), optimize=True)

    # Access the results dictionary
    print(dumps(elastic.property_instances,indent=2))