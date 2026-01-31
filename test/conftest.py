import subprocess
import time

import pytest

MAX_INSTALL_ATTEMPTS = 9999
INSTALL_RETRY_TIMEOUT = 10


@pytest.fixture(scope="session", autouse=True)
def install_models():
    for _ in range(MAX_INSTALL_ATTEMPTS):
        try:
            try:
                subprocess.run(
                    [
                        "kim-api-collections-management",
                        "install",
                        "CWD",
                        "LJ__MD_414112407348_003",
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
        except Exception as e:
            print("Failed to download item with the following exception:\n" + repr(e))
            print(f"Retrying in {INSTALL_RETRY_TIMEOUT} seconds...")
            time.sleep(INSTALL_RETRY_TIMEOUT)
            pass
    for _ in range(MAX_INSTALL_ATTEMPTS):
        try:
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
        except Exception as e:
            print("Failed to download item with the following exception:\n" + repr(e))
            print(f"Retrying in {INSTALL_RETRY_TIMEOUT} seconds...")
            time.sleep(INSTALL_RETRY_TIMEOUT)
            pass
    for _ in range(MAX_INSTALL_ATTEMPTS):
        try:
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
        except Exception as e:
            print("Failed to download item with the following exception:\n" + repr(e))
            print(f"Retrying in {INSTALL_RETRY_TIMEOUT} seconds...")
            time.sleep(INSTALL_RETRY_TIMEOUT)
            pass
