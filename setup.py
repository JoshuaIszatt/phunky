from setuptools import find_packages, setup

setup(
    name="phunky",
    version="0.1",
    description="Python package to assemble phage nanopore reads",
    author="Joshua J Iszatt",
    author_email="joshiszatt@gmail.com",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.10.14",
    ],
    python_requires=">=3.10",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'phunky.py = phunky.main:main',
        ],
    },
)
