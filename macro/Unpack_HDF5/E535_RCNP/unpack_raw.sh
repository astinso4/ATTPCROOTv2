#!/bin/bash

logDir="./logUnpack"

if [ $# -lt 1 ]; then
	echo "Usage: $0 runNum1 [runNum2, ...]"
	exit 1
fi

for runNum in "$@"; do
	logFile="$logDir/RAW_run$runNum.log"
	echo "Processing run: $runNum (>'-')>.h5 "

	root -b -q -l "unpackRawEvents.C($runNum)" > >(tail -n 10000 > $logFile) 2>&1
	#root -b -q -l "unpackRawEvents.C($runNum)"
done

echo ".root<('-'<)\(>.<)/ Unpacking is completed ~('-'~) (~'-')~"
