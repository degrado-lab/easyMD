#!/bin/bash

VERSION=0.1.0

# Define the lines to be added to .bashrc
bindings_setup='
# >>> set up bindings for EasyMD Apptainer >>>
export APPTAINER_BIND=/opt/sge,/wynton/home/opt/sge/common/,/lib64/libssl.so.10,/lib64/libcrypto.so.10,/lib64/libjemalloc.so.1,/lib64/libmunge.so.2,/etc/nsswitch.conf,/etc/services,/var/lib/sgeCA
export APPTAINER_NV=true
# <<< set up bindings for EasyMD Apptainer <<<
'

# Append the lines to the .bashrc file if they are not already present
if ! grep -q "set up bindings for EasyMD Apptainer" ~/.bashrc; then
    echo "$bindings_setup" >> ~/.bashrc
    echo "Bindings for EasyMD Apptainer added to .bashrc"
else
    echo "Bindings for EasyMD Apptainer already present in .bashrc"
fi

# Source the .bashrc file to apply the changes
source ~/.bashrc

# Pull the EasyMD Apptainer image
apptainer pull oras://docker.io/nicholasfreitas/easymd:$VERSION