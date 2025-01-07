from setuptools import setup, find_packages

setup(
    name='easyMD',
    version='3.0.0',    
    description='A set of utilities and notebooks for running OpenMM simulations.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/degrado-lab/easyMD',
    author='Nicholas Freitas',
    author_email='nicholas.freitas@ucsf.edu',
    license='MIT',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'easyMD=easyMD.cli:init_cli',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=[
    ],
)
