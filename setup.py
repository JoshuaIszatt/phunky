from setuptools import find_packages, setup

setup(
    name="phunky",
    version="0.2",
    description="Python package to assemble phage nanopore reads",
    author="Joshua J Iszatt",
    author_email="joshiszatt@gmail.com",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
    ],
    python_requires=">=3.10",
    packages=find_packages(include=['phunky', 'phunky.*']),
)
