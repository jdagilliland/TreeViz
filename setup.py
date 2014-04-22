from distutils.core import setup

setup(
        name='TreeViz',
        version='0.1.0',
        author='J. D. A. Gilliland',
        author_email='jdagilliland@gmail.com',
        py_modules=['tree_viz'],
        scripts=['bin/treeviz'],
        install_requires=[
            'numpy',
            'ete2 >= 2.2',
            'BioPython >= 1.63',
            'matplotlib >= 1.3.1',
            ],
        license='LICENSE.txt',
        description='Utilities to help with visualizing PHYLIP trees',
        long_description=open('README.md').read(),
        )

