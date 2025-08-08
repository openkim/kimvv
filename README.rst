KIM Validation and Verification
===============================

|Testing| |PyPI|

.. |Testing| image:: https://github.com/openkim/kimvv/actions/workflows/test.yml/badge.svg
   :target: https://github.com/openkim/kimvv/actions/workflows/test.yml
.. |PyPI| image:: https://img.shields.io/pypi/v/kimvv.svg
   :target: https://pypi.org/project/kimvv/

This package allows the user to run any `OpenKIM <https://openkim.org/>`_ Test Drivers written using the `kim-tools <https://kim-tools.readthedocs.io>`_ package locally. A "Test Driver" is
a computational protocol that reports one or more material properties using the `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_

List of included Test Drivers:

  * EquilibriumCrystalStructure
  * ElasticConstantsCrystal

Currently, all Test Drivers require the AFLOW software to be installed and in your PATH. See https://kim-tools.readthedocs.io/en/stable/#doc-standalone-installation for installation info.

Basic usage example:
--------------------

.. code-block:: python

    from kimvv import ElasticConstantsCrystal
    from ase.build import bulk
    from json import dumps

    # The Test Driver must be instantiated with an ASE Calculator object
    # or a string indicating a KIM model name
    relax = ElasticConstantsCrystal('LennardJones_Ar')

    # To perform the computation, call the Test Driver object. The first argument
    # to most Test Drivers is the crystal structure to perform the compuation on.
    # To see the additonal arguments, use .printdoc() to print the docstring
    relax.printdoc()

    # Let's compute the elastic constants with the "stress-condensed" method.
    # The crystal structure can be specified as an Atoms object. Any dependencies
    # (e.g. relaxing the crystal structure with EquilibriumCrystalStructure) are
    # automatically run.
    results = relax(bulk('Ar','fcc',5.0), method="stress-condensed")

    # Each Test Driver computes a list of one or more dictionaries, each defining
    # a material property in the format specified by the KIM Properties Framework.
    # The name of the property is in the "property-id" key. See
    # https://openkim.org/properties for the definition of each property.
    print(dumps(results, indent=2))


Usage example 2
---------------
Querying for all DFT-relaxed structures for a given combination of elements in OpenKIM and relaxing them with your potential

.. code-block:: python

    from kimvv import EquilibriumCrystalStructure
    from kim_tools import (
      query_crystal_structures,
      get_deduplicated_property_instances
    )
    from json import dumps
    from ase.calculators.lj import LennardJones

    # Query for all relaxed Argon reference data in OpenKIM
    raw_structs = query_crystal_structures(stoichiometric_species=["Ar"])

    # Deduplicate them
    unique_structs = get_deduplicated_property_instances(raw_structs, allow_rotation=True)

    # Instantiate the Driver with your model
    relax = EquilibriumCrystalStructure(LennardJones(sigma=3.4,epsilon=0.0104,rc=8.15))

    # Run the Driver with each structure. As this is run, the driver internally accumulates
    # Property Instances
    for struct in unique_structs:
      relax(struct)

    # In addition to returning the Property Instances for the current run, Test Drivers
    # accumulate all computed Property Instances. They can be accessed like this:
    print(dumps(relax.property_instances, indent=2))
