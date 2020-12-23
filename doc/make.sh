#!/usr/bin/env bash

# stop if anything goes wrong
set -e

# configure
export DJANGO_SETTINGS_MODULE=genui.settings.test
sphinx-apidoc -o ./source/api/ ../src/genui/ ../src/genui/settings/databases/*prod* ../src/genui/settings/*prod* ../src/genui/settings/*stage*

# make
make html

# upload the docs to GitHub if requested

#if [[ "$@" == *"--upload"* ]]
#then
#    ./upload_docs.sh && echo "Docs uploaded successfully." && exit
#    echo && echo "There was an error during upload. See the messages above for more information." 1>&2
#fi

# serve docs if required
#if [[ "$@" == *"--serve"* ]]
#then
#    python serve.py
#fi
