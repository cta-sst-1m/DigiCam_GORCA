#!/bin/bash

localdir=`dirname $0`
prog=`basename $0`

if [ "$1" == "-h" ] || [ "$1" == "--help" ]
then
    echo
    echo "Clean the garbage from the (non-running) temporary files"
    echo
    echo "Usage:"
    echo "  $prog [ <local_environment> ]"
    echo "Options:"
    echo "  <local_environment>   Local environment where to find the definition of the"
    echo "                        temporary files parent path. [Default: $localdir/local.sh]"
    echo
    exit 1
fi

environ=${1-}
if [ "$environ" == "" ]
then
    environ=$localdir/local.sh
fi
if [ ! -e $environ ]
then
    echo "ERROR: '$environ' not found"
    exit 1
fi
source $environ
nodes=(`pbsnodes -a | grep -v state | grep -v np | grep -v properties | grep -v ntype | grep -v jobs | grep -v status | grep -v mom_service_port | grep -v mom_manager_port`)
#running_jobs_column=(`qstat -u ${USER} | awk '{print $1}'`)
running_jobs=(`qstat -u ${USER} | grep $USER | awk '{print $1}'`) #(${running_jobs_column[@]:4})
for n in ${nodes[@]}
do
    # Check the node accessibility
    ssh -q $n exit
    if [ $? -eq 0 ]
    then
	echo " - cleaning $n node..."
	folders=(`ssh ${n} "ls ${tmp_dir}"`)
	for f in ${folders[@]}
	do
	    isrunning=0
	    for job in ${running_jobs[@]}
	    do
		if [[ "$f" =~ "$job" ]]
		then
		    isrunning=1
		    break
		fi
	    done
	    if [ $isrunning -eq 0 ]
	    then
		ssh ${n} "rm -rf ${tmp_dir}/$f"
	    else
		echo "  * folder $f is still used by a run"
	    fi
	done
    else
	echo " - node $n not accessible"
    fi
done
