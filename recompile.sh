#!/bin/bash

echo "Go to /src directory"
cd src

echo "clean executables"
make clean

echo "make all"
make all

echo "back to root directory"
echo "done"
