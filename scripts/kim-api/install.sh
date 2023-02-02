#! /bin/bash
cmake -S scripts/kim-api -B build/kim-api
cmake --build build/kim-api -j4
 