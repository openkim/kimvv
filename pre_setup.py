import tarfile
import shutil
import os
import pathlib
from urllib.request import urlretrieve
import tomlkit

import kim_edn

openkim_drivers = ["EquilibriumCrystalStructure__TD_457028483760_003"]

with open("pyproject.toml.tpl") as f_pyproject:
    pyproject = tomlkit.parse(f_pyproject.read())


with open("kimvv/__init__.py", "w") as f_init:
    list_of_drivers_literal = "__all__ = ["
    for openkim_driver in openkim_drivers:
        print(f"Importing {openkim_driver}")
        # Get the .txz from OpenKIM as a tmpfile
        url = "https://openkim.org/download/" + openkim_driver + ".txz"
        prefix = "_".join(openkim_driver.split("_")[:-4])
        tmpfile, _ = urlretrieve(url)
        # Extract it and move it to kimvv directory
        with tarfile.open(tmpfile, "r:xz") as f:
            f.extractall()
        final_path_to_driver = os.path.join("kimvv", prefix)
        if os.path.isdir(final_path_to_driver):
            shutil.rmtree(final_path_to_driver)
        shutil.move(openkim_driver, final_path_to_driver)
        # Create an __init__.py in the Test Driver root directory
        # TODO: Rethink this, this brings in additional unwanted .py files
        # if they exist
        pathlib.Path(os.path.join(final_path_to_driver, "__init__.py")).touch()
        # Take care of  kimvv/__init__.py
        f_init.write(
            f"from .{prefix}.test_driver.test_driver import TestDriver as {prefix}\n"
        )
        list_of_drivers_literal += prefix + ","

        requirements_path = os.path.join(final_path_to_driver, "requirements.txt")
        if os.path.isfile(requirements_path):
            pyproject["tool"]["setuptools"]["dynamic"]["dependencies"]["file"].append(
                requirements_path
            )
        else:
            print("NOTE: Importing a Test Driver without a requirements.txt")

        # Should always exist
        kimspec_path = os.path.join(final_path_to_driver, "kimspec.edn")
        kimspec = kim_edn.load(kimspec_path)
        if "developer" not in kimspec:
            print("WARNING: Importing a Test Driver without a kimspec.edn")
        else:
            for developer in kimspec["developer"]:
                pass
                # This is a KIM ID. Need to query for first and last name
        # TODO: Add MANIFEST.in handling

    f_init.write(list_of_drivers_literal + "]\n")

with open("pyproject.toml", "w") as f:
    f.write(tomlkit.dumps(pyproject))
