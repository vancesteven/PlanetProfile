from setuptools import setup

setup(
    name='PlanetProfile',
    version='2.3.18',
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
    packages=['PlanetProfile'],
    package_dir={'PlanetProfile': 'PlanetProfile'},
    install_requires=[
        'numpy >= 1.24.2',
        'scipy >= 1.10.1',
        'mpmath >= 1.2.1',
        'matplotlib >= 3.7.1',
        'SeaFreeze >= 0.9.6',
        'MoonMag >= 1.6.0',
        'gsw >= 3.6.16',
        'spiceypy >= 5.1.1',
        'cmasher >= 1.6.3'
    ],
    include_package_data=True  # Files to include are listed in MANIFEST.in
)
