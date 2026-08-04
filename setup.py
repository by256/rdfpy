import setuptools

with open('./README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='rdfpy',
    version='2.0.0',
    description='rdfpy is a Python package for computing 2D and 3D radial distribution functions.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Batuhan Yildirim',
    author_email='batuhanyild@gmail.com',
    url='https://github.com/by256/rdfpy',
    packages=setuptools.find_packages(exclude=['tests', 'tests.*']),
    python_requires='>=3.9',
    # scipy>=1.6 is the true floor (cKDTree.query_ball_point gained `workers`
    # there); 1.9 is a deliberately conservative floor matching python_requires.
    install_requires=['scipy>=1.9', 'numpy'],
    extras_require={'examples': ['matplotlib']},
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
