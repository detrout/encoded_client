from setuptools import setup, find_packages

setup(
    name="encoded_client",
    description="interface to the ENCODE DCC portal",
    author="Diane Trout",
    author_email="diane@caltech.edu",
    packages=find_packages(),
    # I should be using one or the other package import
    package_data={"": ["*.turtle"]},
    include_package_data=True,
    use_scm_version=True,
    install_requires=[
        "lxml >= 2.2.4",
        "pandas",
        #'keyring',
        "jsonschema",
        "requests",
        "rdflib <= 5.0",
        "xopen",
    ],
)
