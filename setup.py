import setuptools, tarfile, shutil, os, pathlib
from urllib.request import urlretrieve

openkim_drivers = ['ElasticConstantsCrystal__TD_034002468289_000']
kimvv_drivers = []

with open('kimvv/__init__.py','w') as f_init:
    f_init.write('__version__ = "0.1.0"\n')
    list_of_drivers_literal = '__all__ = ['
    for openkim_driver in openkim_drivers:
        url = 'https://openkim.org/download/' + openkim_driver + '.txz'
        prefix = '_'.join(openkim_driver.split('_')[:-4])
        tmpfile, _ = urlretrieve(url)
        with tarfile.open(tmpfile,'r:xz') as f:
            f.extractall()
        final_path_to_driver = os.path.join('kimvv',prefix)
        if os.path.isdir(final_path_to_driver):
            shutil.rmtree(final_path_to_driver)
        shutil.move(openkim_driver,final_path_to_driver)
        pathlib.Path(os.path.join(final_path_to_driver,'__init__.py')).touch()
        kimvv_drivers.append(prefix)
        f_init.write('from .%s.test_driver.test_driver import TestDriver as %s\n'%(prefix,prefix))
        list_of_drivers_literal += prefix+','
    f_init.write(list_of_drivers_literal+']\n')

setuptools.setup(
    name="kimvv",
    version="0.1.0",
    description=(
        "KIM Validation and Verificaiton"
    ),
    author=["ilia Nikiforov, Eric Fuemmeler"],
    install_requires=["kim-tools @ git+https://github.com/openkim/kim-tools.git","numdifftools"],
    packages=setuptools.find_packages(),
    license="CDDL",
)