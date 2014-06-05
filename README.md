# TreeViz

**Jason Gilliland**

These are some utility scripts for visualizing PHYLIP phylogenetic
analysis programs.

To install:
```
python setup.py install
```

It depends on:
- numpy
- matplotlib
- ete2
- BioPython

ete2 in turn depends on:
- PyQt4

ete2 also has optional dependencies on which TreeViz does not also
depend.
If installing ete2 from source, especially from the git repo,
note that it is sometimes installed as `ete_dev` rather than `ete2`.
TreeViz first tries to import `ete2`, but then resorts to `ete_dev` if
`ete2` is not found.

For the benefit of users on MacOSX platforms, who may be using MacPorts
to install the dependencies, I have added a shell profile which may be
used like:
```
source path_correct.profile
```
in order to make sure that the MacPorts python executable is in the
PATH.
