#!/bin/bash

mkdir -p _CPT  _EDR  _GRO  _LOG  _TPR  _trash  _TRR  _XTC  _XVG

mv -v *.cpt _CPT  
mv -v *.edr _EDR  
mv -v *.gro _GRO  
mv -v *.log _LOG  
mv -v *.tpr _TPR  
mv -v *.trr _TRR  
mv -v *.xtc _XTC  
mv -v *.xvg _XVG
mv -v \#* _trash
