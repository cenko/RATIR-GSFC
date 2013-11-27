#
#
#!/bin/bash

#PROJECT_ROOT="$(cd $(dirname "$0"); pwd)"
PROJECT_ROOT="$(cd $(dirname "."); pwd)"

IDL_PATH=+$IDL_DIR/lib:+$PROJECT_ROOT/code/reduction:+$PROJECT_ROOT/code/photometry
export IDL_PATH

export PYTHONPATH=$PYTHONPATH:$PROJECT_ROOT/code/photometry
