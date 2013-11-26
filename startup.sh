#
#
#!/bin/bash

#PROJECT_ROOT="$(cd $(dirname "$0"); pwd)"
PROJECT_ROOT="$(cd $(dirname "."); pwd)"

# add RATIR IDL directories to IDL_PATH
export IDL_PATH=$IDL_PATH:$IDL_DIR/lib:$PROJECT_ROOT/code/idl:$PROJECT_ROOT/code/photometry4

# add RATIR python directories to PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$PROJECT_ROOT/code/photometry4/python

