try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='mscanner',
    version='1.0.1',
    author='Graham Poulter',
    url='http://code.google.com/p/mscanner',
    packages=['mscanner'],
)
