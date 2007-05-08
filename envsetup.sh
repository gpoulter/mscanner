#!/bin/bash

host=$(hostname)

if [ "$host" == "jonbn" ]; then

  py25=/media/jonbn/Python25
  export CHEETAH="$py25/python $(cygpath -m $py25/Scripts/cheetah)"
  export PYTHONPATH='C:\Documents and Settings\Graham\My Documents\data\mscanner\source'
  munge before PATH "$py25"

elif [ "$host" == "eldar" ]; then

  py25=/media/eldar/Python25
  export CHEETAH="$py25/python $(cygpath -m $py25/Scripts/cheetah)"
  export PYTHONPATH='C:\Documents and Settings\Graham\My Documents\data\mscanner\source'
  munge before PATH "$py25"

elif [ "$host" == "maples" ]; then

  export CHEETAH=cheetah
  export PYTHONPATH="$HOME/source/mscanner"

fi
