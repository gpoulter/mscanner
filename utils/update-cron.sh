#!/bin/sh

cd $HOME/data/updates06
wget -r -nd -nc ftp://ftp.nlm.nih.gov/nlmdata/.medlease/gz/
cd $HOME/source/bin
python update.py
python query.py