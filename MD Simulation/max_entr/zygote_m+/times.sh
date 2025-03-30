#!/bin/bash

touch times.dat
echo "a:" >> times.dat

lines=$(grep '^a' bash.log | grep -v 'f')

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $2}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | head -n 1 >> times.dat
    fi
done

echo " " >> times.dat
echo "b:" >> times.dat

lines=$(grep 'b' bash.log)

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $2}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | head -n 1 >> times.dat
    fi
done

echo " " >> times.dat
echo "f:" >> times.dat

lines=$(grep 'f' bash.log)

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $2}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | head -n 1 >> times.dat
    fi
done

echo " " >> times.dat
echo "---------------------------------------------" >> times.dat
echo " " >> times.dat
