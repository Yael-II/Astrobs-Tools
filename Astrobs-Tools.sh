#!/bin/bash

clear
output=$(pip freeze | grep 'astropy==')
version=${output:9:1}
if [ $version -ge 6 ]; then
     python source/SCOPE_v2-1.py
else
     printf "\033[91mError: Please use astropy 6.0 or newer version\033[0m\n"
fi 
