#! /bin/bash

mkdir ${1}/sized
find $1 -size +400k | xargs -I file mv file ${1}/sized
