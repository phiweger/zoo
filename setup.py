from setuptools import setup

setup(
    name='zoo',
    version='0.3',
    description='A distributed microbial database',
    url='https://github.com/viehwegerlib/zoo',
    author='Adrian Viehweger',
    author_email='none',
    license='BSD 3-clause',
    packages=['zoo'],
    package_dir={'zoo': 'zoo'},
    include_package_data=True,
    package_data={'zoo': [
        'cli/*',
        'data/*',
        'data/zika/*',
        'data/tests/*',
        'data/flu/*',
        'data/ebola/*',
        'data/plum/*',
        'data/rna_virome_shi2016/*',
        'schema/*.json',
        'schema/fragments/*.json',
        'schema/specific/*.json',
        'schema/test_schemas/*.json'
        ]},
    install_requires=[
        'Click',
        'deepdiff',
        'ijson',  # likely not needed
        'jsonschema',
        'khmer',
        'networkx',
        'numpy',
        'biopython',
        'pandas',
        'progressbar2',
        'pyfaidx',
        'pymongo',
        'sourmash'
    ],
    zip_safe=False,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    entry_points={
        'console_scripts': [
            'zoo = zoo.__main__:cli'
        ]})


'''
> MANIFEST.in tells Distutils what files to include in the source distribution
but it does not directly affect what files are installed. For that you need to
include the appropriate files in the setup.py file, generally either as package
data or as additional files. -- stackoverflow, 3596979

https://docs.python.org/3/distutils/setupscript.html#installing-package-data
https://docs.python.org/3/distutils/sourcedist.html#manifest
'''
