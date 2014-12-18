#
#
#!/bin/bash

#PROJECT_ROOT="$(cd $(dirname "$0"); pwd)"
#PROJECT_ROOT="$(cd $(dirname "."); pwd)"
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # allows startup.sh to be called from any directory.  will fail if last component of path to startup.sh is a symlink.

# add pipeline directories to path
export PATH=$PATH:$PROJECT_ROOT/code/sdss

# add pipeline directories to idl path
IDL_PATH=$IDL_PATH:+$PROJECT_ROOT/code/reduction:+$PROJECT_ROOT/code/photometry
export IDL_PATH

# add pipeline directories to python path
export PYTHONPATH=$PYTHONPATH:$PROJECT_ROOT/code/:$PROJECT_ROOT/code/photometry:$PROJECT_ROOT/code/photometry/dependencies:$PROJECT_ROOT/code/reduction/astrom:$PROJECT_ROOT/code/reduction/

# set pipeline aliases
alias pcd='cd $PROJECT_ROOT'