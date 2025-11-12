from setuptools import setup, find_packages

setup(
    name="gw_agn_watcher",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "astropy",
        # add other dependencies here
    ],
    python_requires='>=3.8',
    author="Hemanth Kumar",
    author_email="hemanth.bommireddy195@gmail.com",
    description="Software to crossmatch and filter LVK and ZTF alerts via ALeRCE",
    url="https://github.com/Hemanthb1/GW_AGN_Watcher",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
