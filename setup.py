from setuptools import find_packages, setup

with open("README.md","r") as f:
    long_description = f.read()
    
setup(
    name="oscar",
    version="1.3.1",
    description="A collection of classes and functions to read Finite Element Hdf5 files.",
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jjcremmers/oscar",
    author="Joris Remmers",
    author_email="j.j.c.remmers@tue.nl",    
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ], 
    install_requires=["h5py >= 3.6" , "vtk >= 9.2" ],
    python_requires='>=3.10',
)

