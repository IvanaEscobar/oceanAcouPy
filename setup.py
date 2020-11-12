 from setuptools import setup

 setup(
   name='oceanAcouPy',
   version='0.0',
   author='Ivana Escobar',
   author_email='ivana@utexas.edu',
   packages=['oceanAcouPy'],
   url='https://github.com/IvanaEscobar/oceanAcouPy',
   license='LICENSE.txt',
   description='EE 348N: Ocean Acoustics Package',
   long_description=open('README.txt').read(),
   install_requires=[
       "numpy >= 1.19.1",
       "pandas",
   ],
)
