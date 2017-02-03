#!/bin/sh

RC=1
while [ $RC -ne 12341234 ]; do
   julia generateData.jl
   RC=$?
done
