#!/bin/bash


# found= "$(find . -maxdepth 1 -name "ipynb" -type d)"

# if [ -z $found ] ; then
#     echo not found
# else
#     echo found
# fi

ipynb_found=$(find . -maxdepth 1 -name "Notebooks" -type d)
if [ -z $ipynb_found ] ; then
    mkdir ipynb
fi

textfiles_found=$(find . -maxdepth 1 -name "textfiles" -type d)
if [ -z $textfiles_found ] ; then
    mkdir textfiles
fi

varfiles_found=$(find . -maxdepth 1 -name "varfiles" -type d)
if [ -z $varfiles_found ] ; then
    mkdir varfiles
fi
find . -maxdepth 1 -name "code" -type d
code_found=$(find . -maxdepth 1 -name "code" -type d)
if [ -z $code_found ] ; then
    mkdir code
fi

find . -maxdepth 1 -name "*.txt" | xargs -I {} mv {} textfiles
find . -maxdepth 1 -name "*.var" | xargs -I {} mv {} varfiles
find . -maxdepth 1 -name "*.py" | xargs -I {} mv {} code
find . -maxdepth 1 -name "*.ipynb" | xargs -I {} mv {} Notebooks