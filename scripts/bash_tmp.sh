#!/bin/bash

directory=mydir

while getopts "hd:l:" opt; do
    case $opt in
        d ) directory=$OPTARG;;
        l ) depth=-a;;
        h ) usage
        exit 0;;
        *) usage
        exit 1;;
    esac
done
echo "<$directory> <$depth>"
