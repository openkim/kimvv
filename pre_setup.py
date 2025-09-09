import json
import os
import shutil
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.error import HTTPError
from urllib.request import urlopen, urlretrieve

import kim_edn
import tomlkit

# Dictionary of production Test Drivers and kwargs to pass to them during testing
OPENKIM_TEST_DRIVERS = {
    "EquilibriumCrystalStructure__TD_457028483760_003": {"steps": 10},
    "ElasticConstantsCrystal__TD_034002468289_001": {"method": "stress-condensed"},
    "CrystalStructureAndEnergyVsPressure__TD_685283176869_000": {
        "pressure_step_gpa": 5,
        "exp_step": 1.5,
        "min_steps": 10,
    },
}

# List of URLs of development Test Drivers to test
DEVEL_TEST_DRIVERS = {}

MAX_URLLIB_ATTEMPTS = 10


def create_init(td_root_path: os.PathLike):
    for name in os.listdir(td_root_path):
        if name.endswith(".py"):
            raise RuntimeError(
                "Test Driver releases should not have any files with a .py "
                "extension in their root directory (the one containing kimspec.edn). "
                "See https://git-scm.com/docs/gitattributes#_export_ignore for info on "
                "excluding files from your release archive."
            )
    # Create an __init__.py in the Test Driver root directory
    Path(os.path.join(td_root_path, "__init__.py")).touch()


def move_driver(prefix: str, source_dir: os.PathLike):
    final_path_to_driver = os.path.join("kimvv", prefix)
    if os.path.isdir(final_path_to_driver):
        shutil.rmtree(final_path_to_driver)
    shutil.move(source_dir, final_path_to_driver)
    create_init(final_path_to_driver)


def urlretrieve_with_retries(url) -> str:
    for _ in range(MAX_URLLIB_ATTEMPTS):
        try:
            tmpfile, _ = urlretrieve(url)
            return tmpfile
        except HTTPError:
            pass
    raise RuntimeError(
        f"Failed to download {url}" f" after {MAX_URLLIB_ATTEMPTS} attempts."
    )


if __name__ == "__main__":
    kimvv_test_drivers = {}

    # Download and untar production OpenKIM TDs
    for test_driver in OPENKIM_TEST_DRIVERS:
        print(f"Importing {test_driver}")
        # Get the .txz from OpenKIM as a tmpfile
        url = "https://openkim.org/download/" + test_driver + ".txz"
        prefix = "_".join(test_driver.split("_")[:-4])
        # Store the dict of kwargs in the "kimvv_test_drivers" dictionary
        kimvv_test_drivers[prefix] = OPENKIM_TEST_DRIVERS[test_driver]
        tmpfile = urlretrieve_with_retries(url)
        # Extract it and move it to kimvv directory
        with tarfile.open(tmpfile, "r:xz") as f:
            f.extractall()
        move_driver(prefix, test_driver)

    # Download and untar development TDs
    for test_driver in DEVEL_TEST_DRIVERS:
        tmpfile = urlretrieve_with_retries(test_driver)
        # Extract it to a temporary directory
        with TemporaryDirectory() as tmpdir:
            with tarfile.open(tmpfile, "r:gz") as f:
                f.extractall(tmpdir)
            # Find the kimspec.edn
            kimspec_match = Path(tmpdir).rglob("kimspec.edn")
            kimspec_path = next(kimspec_match)
            td_root_path = kimspec_path.parent
            try:
                next(kimspec_match)
                raise RuntimeError("More than one kimspec.edn found in your archive")
            except StopIteration:
                pass

            # Figure out what we will name this
            with open(kimspec_path) as f:
                kimspec = kim_edn.load(f)
            if "extended-id" in kimspec:
                td_name = "__".join(kimspec["extended-id"].split("__")[:-1])
            else:
                td_name = "InProgress"

            # Can easily end up with multiple InProgress
            if td_name in kimvv_test_drivers:
                td_name_prefix = td_name
                for i in range(10):
                    td_name = td_name_prefix + str(i)
                    if td_name not in kimvv_test_drivers:
                        break
                if td_name in kimvv_test_drivers:
                    raise RuntimeError("Somehow got 11 TDs with the same name")

            print(f"Importing from \n{test_driver}\n and naming it {td_name}")
            kimvv_test_drivers[td_name] = DEVEL_TEST_DRIVERS[test_driver]
            move_driver(td_name, td_root_path)

    with open("test/test_inputs.json", "w") as f:
        json.dump(kimvv_test_drivers, f)

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
                for i in range(MAX_URLLIB_ATTEMPTS):
                    try:
                        with urlopen(
                            f"https://openkim.org/profile/{developer}.json"
                        ) as u:
                            developer_profile = json.load(u)
                        break
                    except HTTPError:
                        pass
                if i == MAX_URLLIB_ATTEMPTS - 1:
                    raise RuntimeError(
                        "Failed to get developer info from openkim.org for "
                        f"{developer} after {MAX_URLLIB_ATTEMPTS} attempts."
                    )
                name = (
                    developer_profile["first-name"]
                    + " "
                    + developer_profile["last-name"]
                )
                if any(
                    name.lower() in author["name"].lower()
                    for author in pyproject["project"]["authors"]
                ):
                    continue
                pyproject["project"]["authors"].append({"name": name})

        manifest_path = os.path.join(driver_path, "MANIFEST.in")
        if os.path.isfile(manifest_path):
            with open(manifest_path) as f:
                manifest_td_lines = f.readlines()
            for line in manifest_td_lines:
                command = line.split()[0]
                path = line.split()[1]
                path_full = os.path.join(driver_path, path)
                manifest_kimvv += f"{command} {path_full}\n"
        # TODO: Build docs by combining READMES

    # write __init__.py
    with open("kimvv/__init__.py", "w") as f:
        # import call decorator
        f.write("from .core import KIMVVTestDriver, override_call_method\n")
        for td in kimvv_test_drivers:
            f.write(f"from .{td}.test_driver.test_driver import TestDriver as __{td}\n")

        f.write("\n\n")

        for td in kimvv_test_drivers:
            # add call decorator
            f.write("@override_call_method\n")
            f.write(f"class {td}(__{td}, KIMVVTestDriver):\n    pass\n\n\n")

        f.write("__all__ = [\n")
        for td in kimvv_test_drivers:
            f.write(f'    "{td}",\n')
        f.write("]\n")

    with open("MANIFEST.in", "w") as f:
        f.write(manifest_kimvv)

    with open("pyproject.toml", "w") as f:
        f.write(tomlkit.dumps(pyproject))
