#!/bin/sh
cd "$(dirname "$0")"; ./mono $2 $1 $3 5 out/J$1_L$2.bin 2>> out/log.txt
