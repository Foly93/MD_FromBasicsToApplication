conda create --no-default-packages -n ppis python=3.7 -y
conda install -n ppis -c rdkit rdkit -y
conda install -n ppis -c anaconda pandas numpy seaborn jupyter -y 
conda install -n ppis -c conda-forge mplcursors ipympl -y
conda install -n ppis -c salilab modeller -y > modellerconfig.txt
conda install -n ppis -c conda-forge fpocket mdtraj nglview ambertools=23 -y
conda install -n ppis -c conda-forge jupyterlab jupyter_contrib_nbextensions -y
conda install -n ppis -c conda-forge -c schrodinger pymol-bundle -y
conda install -n ppis -c bioconda autodock-vina autogrid -y
conda install -n ppis -c conda-forge mdanalysis
