#!/bin/bash 

# Run Matlab in batch mode.
# Invoking Matlab with the -nojvm -nodisplay and
# -nosplash ensures that there is no graphical display.
#
# Ning Ma, 18 March 2007

if [ $# -lt 1 ]
then
    echo "Usage: run_matlab.sh MATLAB_CMD [MATLAB_CMD_DIR]"
    echo "       MATLAB_CMD     - Matlab command with arguments, enclosed by \"\"."
    echo "       MATLAB_CMD_DIR - Directory where the Matlab command can be found. Optional."
    echo ""
    echo "eg. run_matlab.sh \"max([1 3 2])\""
    exit 1
fi

MATLABBIN="/usr/local/bin/matlab"

MATLABCMDS="$1; exit;"
if [ $# -eq 2 ]
then
    MATLABCMDS="addpath $2; $MATLABCMDS"
fi

$MATLABBIN -nodesktop -nodisplay -r "$MATLABCMDS" #&> /dev/null

