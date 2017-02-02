#!/bin/bash

sol=$1

for t in tests/*.in; do
    echo "#### $t"
    echo "------- true answer -------"
    ot=`echo "$t" | sed -r 's|.in$|.out|'`
    cat "$ot"
    echo
    echo "------- our answer  -------"
    python3 "$sol" < "$t"
    echo
done
