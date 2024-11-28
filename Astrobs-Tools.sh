#!/bin/bash

source activate.sh
if ! command -v python3 2>&1 >/dev/null; then
     printf "\033[91mError: Please install python3 or newer version\033[0m\n"
else
     version=$(pip3 freeze | sed -r "s/==|\./ /g" | awk '$1 == "astropy" {print $2}')
     if [ -z $version ]; then
	  printf "\033[91mError: Astropy 6.0 is required\033[0m\n"
     elif [ $version -lt 6 ]; then
	  printf "\033[91mError: Please install astropy 6.0 or newer version\033[0m\n"
     else
	  clear
	  python Source/ASTROBS_SELECT_v1.py
     fi 
fi
