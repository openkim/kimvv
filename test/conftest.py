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
                "EAM_Dynamo_Ackland_1987_Au__MO_754413982908_001",
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
