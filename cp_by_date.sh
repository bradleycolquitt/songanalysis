#! /bin/bash

if [ ! -d $1 ]
then
    mkdir $1
    fi

find . -type f -name "*wav" -maxdepth 1 -newermt $1 ! -newermt `date -d "$1+1 days" +%Y-%m-%d` | xargs -I f cp --no-clobber --preserve f $1
