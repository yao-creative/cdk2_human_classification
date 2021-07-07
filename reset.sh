#!/bin/bash

find . -name "*.gz" | xargs rm
find . -name "*.pdb" | xargs rm
find . -name "*.txt" | xargs rm
find . -name "*.var" | xargs rm