#!/bin/sh
cd "$(dirname "$0")"; ./worm 64 128 $1 0 $3 1.5e9 7 out/J$1_$2.bin 2>> out/log.txt
