from setuptools import setup

setup(
    name='PlanetProfile',
    version='2.5.0',
    author='Marshall J. Styczinski',
    author_email='itsmoosh@gmail.com',
    description='Self-consistent geophysical models for large moons and ocean worlds',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/vancesteven/PlanetProfile',
    project_urls={
        'Bug tracker': 'https://github.com/vancesteven/PlanetProfile/issues',
        'Publication': 'https://doi.org/10.1029/2022EA002748'
    },
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent'
    ],
    packages=['PlanetProfile', 'PlanetProfile.TrajecAnalysis'],
    package_dir={'PlanetProfile': 'PlanetProfile',
                 'PlanetProfile.TrajecAnalysis': 'PlanetProfile/TrajecAnalysis'},
    install_requires=[
        'numpy >= 1.24.4',
        'scipy >= 1.11.4',
        'mpmath >= 1.3.0',
        'matplotlib >= 3.8.2',
        'SeaFreeze >= 0.9.6',
        'MoonMag >= 1.7.5',
        'gsw >= 3.6.16',
        'spiceypy >= 6.0.0',
        'cmasher >= 1.6.3',
        'hdf5storage >= 0.1.19'
    ],
    include_package_data=True  # Files to include are listed in MANIFEST.in
)
