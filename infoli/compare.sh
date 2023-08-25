#!/bin/bash

if [ -z "$1" ]
then
    SIZE=100
else
    SIZE=$1
fi

echo ""
echo "---------------------------------------"
echo "Serial Version..."
echo "---------------------------------------"
echo ""

cd serial
mkdir -p results
make
./infoli_simple ${SIZE}

echo ""
echo "---------------------------------------"
echo "Parallel Version..."
echo "---------------------------------------"
echo ""

cd ../parallel
mkdir -p results
make
./infoli_simple ${SIZE}

echo ""
echo "---------------------------------------"
echo "Comparing results..."
echo "---------------------------------------"
echo ""

cd ..

if diff --brief serial/results/lastStateDump.txt parallel/results/lastStateDump.txt >/dev/null 2>&1; then
  echo "Checks passed: output is equal."
else
  echo "Checks failed: resulting output is different."
fi
echo ""