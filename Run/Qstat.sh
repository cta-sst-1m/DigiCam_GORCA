#!/bin/bash

prog=`basename $0`

if [ "$1" == "-h" ] || [ "$1" == "--help" ]
then
    echo
    echo "Print the summary of the number of jobs running, waiting, holding or in error"
    echo
    echo "Usage:"
    echo "  $prog [<user>]"
    echo "Options:"
    echo "  <user>      User name. [Default: $USER]"
    echo
    exit 1
fi

user=${1-$USER}

list=(`qstat -u ${user} | awk '{print $10}'`)
c=`printf -- '%s\n' "${list[@]}" | grep C | wc -l`
r=`printf -- '%s\n' "${list[@]}" | grep R | wc -l`
qw=`printf -- '%s\n' "${list[@]}" | grep Q | wc -l`
hqw=`printf -- '%s\n' "${list[@]}" | grep H | wc -l`
Eqw=`printf -- '%s\n' "${list[@]}" | grep E | wc -l`

echo
echo "Jobs in Queue $((r + qw + hqw + Eqw))"
echo "   (Completed $c)"
echo
echo "      running $r"
echo "      waiting $qw"
echo "      holding $hqw"
echo "        error $Eqw"
echo
