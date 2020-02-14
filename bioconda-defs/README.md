# Image definitions for BioConda tools

Definition files to build images for tools available from BioConda channel.

Divided in directories by *releases* named as *{year}.{release_number}*.

## Scripts

* `conda2def.py -c CHANNEL -p PACKAGE=VERSION -x ENTRYPOINT` - to produce a definition file using conda
* `build_tool.py listfile` -  will create a set of definition files starting from a list in the "PACKAGE=VERSION tab ENTRYPOINT" format
