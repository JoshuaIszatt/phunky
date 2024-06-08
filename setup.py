from setuptools import find_packages, setup
import os

# Set the base directory location which is the location of the setup.py file
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

setup(
    name="phunky",
    version="0.3",
    description="Python package to assemble phage nanopore reads",
    author="Joshua J Iszatt",
    author_email="joshiszatt@gmail.com",
    classifiers=[
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.10",
    packages=find_packages(include=['phunky', 'phunky.*']),
    package_data={
        'phunky': ['logging.json'],
    },
)
