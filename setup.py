import setuptools, subprocess

# As TDs are published, we will download from OpenKIM.org instead of cloning from Github
subprocess.run(['git','clone','https://github.com/openkim-hackathons/ElasticConstants__TD_000000000000_000.git','kimvv/ElasticConstants'])
subprocess.run(['touch','kimvv/ElasticConstants/__init__.py'])

setuptools.setup(
    name="kimvv",
    version="0.1.0",
    description=(
        "KIM Validation and Verificaiton"
    ),
    author=["ilia Nikiforov, Eric Fuemmeler"],
    install_requires=["kim-tools"],
    packages=setuptools.find_packages(),
    license="CDDL",
)