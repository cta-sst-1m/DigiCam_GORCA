#!/bin/bash

prog=`basename $0`

if [ $# -lt 1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]
then
    echo
    echo "Delete job(s)"
    echo
    echo "Usage:"
    echo "  $prog <firstID> [ <lastID> ]"
    echo
    echo "Options:"
    echo "  <firstID>   Number ID of the PSB job to be deleted as first. If no <stopID> is given,"
    echo "              only <startID> is deleted. To delete all, use 'all'."
    echo "  <lastID>    If given, all jobs from <startID> to <stopID> (included) will be deleted"
    echo
    exit 1
fi

firstID=${1}
lastID=${2-}

if [ "$lastID" == "" ]
then
    if [ "$firstID" == "all" ]
    then
	echo " All jobs are going to be deleted"
	qdel 'all'
    else
	echo " The job with ID '$firstID' is being deleted"
	qdel $firstID
    fi
else
    echo " The jobs with IDs from '$firstID' to '$lastID' are being deleted"
    qdel `seq $firstID $lastID`
fi
echo "Done."
