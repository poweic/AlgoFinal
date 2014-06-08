#!/bin/bash -e

DIFF=$(diff $1 $2 || true)

RESULT=$(echo "$DIFF" | sed "s%000 .*-%000 -%g" | grep "\---" -v | sed "1~3d" | sed "s%^<\|^>%%g")

odd=$(echo "$RESULT" | sed "1~2d")
even=$(echo "$RESULT" | sed "2~2d")

if [ "$odd" != "$even" ]; then
  printf "\33[31m[FAILED] $1 and $2 differs \33[0m\n"
  exit -1
else
  printf "\33[32m[PASSED] congrats !! \33[0m\n"
fi
