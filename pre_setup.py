import json
import os
import pathlib
import shutil
import tarfile
from urllib.request import urlopen, urlretrieve

import kim_edn
import tomlkit

# List of production Test Drivers
OPENKIM_TEST_DRIVERS = ["EquilibriumCrystalStructure__TD_457028483760_003"]


def create_init(td_root_path: os.PathLike):
    # Create an __init__.py in the Test Driver root directory
    # TODO: Rethink this, this brings in additional unwanted .py files
    # if they exist
    pathlib.Path(os.path.join(td_root_path, "__init__.py")).touch()


if __name__ == "__main__":
    kimvv_test_drivers = []

    # Download and untar production OpenKIM TDs
    for test_driver in OPENKIM_TEST_DRIVERS:
        print(f"Importing {test_driver}")
        # Get the .txz from OpenKIM as a tmpfile
        url = "https://openkim.org/download/" + test_driver + ".txz"
        prefix = "_".join(test_driver.split("_")[:-4])
        kimvv_test_drivers.append(prefix)
        tmpfile, _ = urlretrieve(url)
        # Extract it and move it to kimvv directory
        with tarfile.open(tmpfile, "r:xz") as f:
            f.extractall()
        final_path_to_driver = os.path.join("kimvv", prefix)
        if os.path.isdir(final_path_to_driver):
            shutil.rmtree(final_path_to_driver)
        shutil.move(test_driver, final_path_to_driver)
        create_init(final_path_to_driver)

    with open("pyproject.toml.tpl") as f_pyproject:
        pyproject = tomlkit.parse(f_pyproject.read())
    manifest_kimvv = ""

    # Should by now all be in directories named `kimvv/{test_driver}`
    for test_driver in kimvv_test_drivers:
        driver_path = os.path.join("kimvv", test_driver)
        requirements_path = os.path.join(driver_path, "requirements.txt")
        if os.path.isfile(requirements_path):
            pyproject["tool"]["setuptools"]["dynamic"]["dependencies"]["file"].append(
                requirements_path
            )
        else:
            print("NOTE: Importing a Test Driver without a requirements.txt")

        # Kimspec should always exist
        kimspec_path = os.path.join(driver_path, "kimspec.edn")
        kimspec = kim_edn.load(kimspec_path)
        manifest_kimvv += "include " + kimspec_path + "\n"
        if "developer" not in kimspec:
            print("WARNING: Importing a Test Driver without any developers")
        else:
            # Look up developer on openkim.org
            for developer in kimspec["developer"]:
                with urlopen(f"https://openkim.org/profile/{developer}.json") as u:
                    developer_profile = json.load(u)
                    name = (
                        developer_profile["first-name"]
                        + " "
                        + developer_profile["last-name"]
                    )
                    if any(
                        name.lower() == author["name"].lower()
                        for author in pyproject["project"]["authors"]
                    ):
                        continue
                    pyproject["project"]["authors"].append({"name": name})

        manifest_path = os.path.join(driver_path, "MANIFEST.in")
        if os.path.isfile(manifest_path):
            with open(manifest_path) as f:
                manifest_td = f.read()
            manifest_kimvv += manifest_td + "\n"
        # TODO: Build docs by combining READMES

    # write __init__.py
    with open("kimvv/__init__.py", "w") as f:
        f.write("from .core import KIMVVTestDriver\n")
        for td in kimvv_test_drivers:
            f.write(f"from .{td}.test_driver.test_driver import TestDriver as __{td}\n")

        f.write("\n\n")

        for td in kimvv_test_drivers:
            f.write(f"class {td}(__{td}, KIMVVTestDriver):\n    pass\n\n\n")

        f.write("__all__ = [\n")
        for td in kimvv_test_drivers:
            f.write(f'    "{td}",\n')
        f.write("]\n")

    with open("MANIFEST.in", "w") as f:
        f.write(manifest_kimvv)

    with open("pyproject.toml", "w") as f:
        f.write(tomlkit.dumps(pyproject))
