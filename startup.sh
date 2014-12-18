#
#
#!/bin/bash

##########################
# check for dependencies #
##########################
echo "Checking RATIR pipeline dependencies:"

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

# set project root directory
export RAT_PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # allows startup.sh to be called from any directory.  will fail if last component of path to startup.sh is a symlink.

# add pipeline directories to path
export PATH=$RAT_PROJECT_ROOT/code/sdss:$PATH

# add pipeline directories to idl path
export IDL_PATH=+$RAT_PROJECT_ROOT/code/reduction:+$RAT_PROJECT_ROOT/code/photometry:$IDL_PATH

# add pipeline directories to python path
export PYTHONPATH=$RAT_PROJECT_ROOT/code:$RAT_PROJECT_ROOT/code/photometry:$RAT_PROJECT_ROOT/code/photometry/dependencies:$RAT_PROJECT_ROOT/code/reduction/astrom:$RAT_PROJECT_ROOT/code/reduction:$PYTHONPATH

# set pipeline aliases
alias cd_rat='cd $RAT_PROJECT_ROOT'