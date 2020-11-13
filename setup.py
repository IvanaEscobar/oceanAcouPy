 from setuptools import setup
 
 setup(
   name='oceanAcouPy',
   version='0.0',
   packages=['oceanAcouPy', 'oceanAcouPy.test'],
   url='https://github.com/IvanaEscobar/oceanAcouPy',
   keywords=['acoustics','underwater','ocean'],
   license='LICENSE.txt',
   description='EE 348N: Ocean Acoustics Package',
   long_description=open('README.txt').read(),
   long_description_content_type='text/markdown'
   install_requires=[
       "numpy >= 1.19.1",
       "pandas",
   ],
   author='Ivana Escobar',
   author_email='ivana@utexas.edu'
)
