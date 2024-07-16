from tqdm import tqdm 
import MDAnalysis as mda

top='prod_1.tpr'
trj='prod_center_1.xtc'
out='prod_center_1.nc'

u=mda.Universe(top,trj)

protein = u.select_atoms("protein or resname UNL")

with mda.Writer(out, protein.n_atoms) as W:
    for ts in tqdm(u.trajectory):
        W.write(protein)
