#!/bin/bash

mkdir -p _output _trash

mv -v umbrella*.gro umbrella*.tpr *.trr *.xtc umbrella*.xvg _output
mv -v *.cpt *.log conf*.gro npt*.gro npt*.tpr npt*.xvg *.edr em.tpr ions.tpr pull.tpr _trash
mv -v \#* _trash
