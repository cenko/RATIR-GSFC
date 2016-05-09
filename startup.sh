#
#
#!/bin/bash

##########################
# check for dependencies #
##########################
echo "Checking pipeline dependencies:"

# IDL
command -v idl >/dev/null 2>&1 || { echo >&2 "IDL is required but it's not installed.  Aborting."; exit 1; }
echo -en "\t* IDL located:\n\t\t"
command -v idl

# SExtractor
command -v sex >/dev/null 2>&1 || { echo >&2 "SExtractor is required but it's not installed.  Aborting."; exit 1; }
echo -en "\t* SExtractor located:\n\t\t"
command -v sex

# Swarp
command -v swarp >/dev/null 2>&1 || { echo >&2 "Swarp is required but it's not installed.  Aborting."; exit 1; }
echo -en "\t* Swarp located:\n\t\t"
command -v swarp

# Scamp
command -v scamp >/dev/null 2>&1 || { echo >&2 "Scamp is required but it's not installed.  Aborting."; exit 1; }
echo -en "\t* Scamp located:\n\t\t"
command -v scamp

# MissFITS
command -v missfits >/dev/null 2>&1 || { echo >&2 "MissFITS is required but it's not installed.  Aborting."; exit 1; }
echo -en "\t* MissFITS located:\n\t\t"
command -v missfits

# cdsclient
command -v findcat >/dev/null 2>&1 || { echo >&2 "cdsclient is required but it's not installed.  Aborting."; exit 1; }
echo -en "\t* cdsclient located:\n\t\t"
command -v findcat

######################
# set up environment #
######################
echo "Setting up pipeline environment."

# set project root directory
export PIPE_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # allows startup.sh to be called from any directory.  will fail if last component of path to startup.sh is a symlink.

# add pipeline directories to path
export PATH=$PIPE_ROOT/code/sdss:$PATH

# add pipeline directories to idl path
export IDL_PATH=+$PIPE_ROOT/code/reduction:+$PIPE_ROOT/code/photometry:$IDL_PATH

# add pipeline directories to python path
export PYTHONPATH=$PIPE_ROOT/code:$PIPE_ROOT/code/photometry:$PIPE_ROOT/code/photometry/dependencies:$PIPE_ROOT/code/reduction/astrom:$PIPE_ROOT/code/reduction:$PIPE_ROOT/code/python_reduction/auto:$PYTHONPATH

echo "Pipeline ready."