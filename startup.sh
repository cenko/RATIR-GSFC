#
#
#!/bin/bash

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