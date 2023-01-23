from setuptools import setup

setup(
    name='OptiPRISMS',
    version='1.0.0',
    packages=[''],
    package_dir={'': 'src'},
    url='https://github.com/DorianDepriester/OptiPRISMS',
    license='GNU Lesser General Public License v2.1',
    author='Dorian Depriester',
    author_email='dorian [dot] depriester [at] ensam [guess what] eu',
    description='Indentifies Crystal Plasticity (CP) parameters by inverse analyis based on CPFEM simulations performed using PRISMS-Plasticity',
    install_requires=[
        'scipy',
        'numpy',
        'vtk',
        'pandas',
        'optimparallel'
    ],
)
