import subprocess

import pytest


@pytest.fixture(scope="session", autouse=True)
def install_models():
    try:
        subprocess.run(
            [
                "kim-api-collections-management",
                "install",
                "CWD",
                "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003",
            ],
            check=True,
            capture_output=True,
            encoding="utf-8",
        )
    except subprocess.CalledProcessError as e:
        if ("already installed" in e.output) or (
            "does not appear to contain CMakeLists.txt." in e.output
        ):
            pass
        else:
            raise e
    try:
        subprocess.run(
            [
                "kim-api-collections-management",
                "install",
                "CWD",
                "Sim_LAMMPS_ADP_StarikovGordeevLysogorskiy_2020_SiAuAl__SM_113843830602_000",  # noqa:E501
            ],
            check=True,
            capture_output=True,
            encoding="utf-8",
        )
    except subprocess.CalledProcessError as e:
        if ("already installed" in e.output) or (
            "does not appear to contain CMakeLists.txt." in e.output
        ):
            pass
        else:
            raise e
