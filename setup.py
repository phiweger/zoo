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
    include_package_data=True,  # use Manifest.in, stackoverflow, 13307408
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
        'pymongo'
        # 'sourmash'
    ],
    dependency_links=[
        "https://github.com/dib-lab/sourmash/tarball/master#egg=sourmash"
        # version (2017-05-03), 2.0.0a1
        # https://mike.zwobble.org/2013/05/adding-git-or-hg-or-svn-dependencies-in-setup-py/
        # include package as tarball, stackoverflow, 32688688
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
