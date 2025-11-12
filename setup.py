from setuptools import setup, find_packages

setup(
    name="gw-agn-WATCHER",
    version="0.1.0",
    author="Hemanth Kumar",
    author_email="hemanth.bommireddy195@gmail.com",
    description="Python tools for GW AGN follow-up",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/GW_AGN_watcher",
    packages=find_packages(),
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        "numpy",
        "astropy",
        "matplotlib",
        "scipy",
        "requests",
        "pandas",
        "astroquery"
    ],
)
