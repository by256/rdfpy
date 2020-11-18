import setuptools

with open('./README.md', 'r') as f:
    long_description = f.read()

setuptools.setup(
    name='rdfpy',
    version='0.1.5',
    description='rdfpy is a Python package for computing 2D and 3D radial distribution functions.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Batuhan Yildirim',
    author_email='by256@cam.ac.uk',
    url='https://github.com/by256/rdfpy',
    packages=setuptools.find_packages(),
    install_requires=['scipy', 'numpy', 'matplotlib'],
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
     )