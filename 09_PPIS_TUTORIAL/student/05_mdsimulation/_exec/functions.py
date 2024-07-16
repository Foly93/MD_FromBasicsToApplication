import os
import sys
import re
import pytraj as pt
import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import nglview as nv
from mdtraj.geometry import _geometry
from mdtraj.utils import ensure_type

# Parateters
_ATOMIC_RADII = {'H'   : 0.120, 'He'  : 0.140, 'Li'  : 0.076, 'Be' : 0.059,
                 'B'   : 0.192, 'C'   : 0.170, 'N'   : 0.155, 'O'  : 0.152,
                 'F'   : 0.147, 'Ne'  : 0.154, 'Na'  : 0.102, 'Mg' : 0.086,
                 'Al'  : 0.184, 'Si'  : 0.210, 'P'   : 0.180, 'S'  : 0.180,
                'Ca'   : 0.200, 'Cl'  : 0.200, 'Zn'  : 0.200, 'VS' : 0.152,
                 'I'   : 0.220, 'Br'  : 0.196 }


# SASA Calculation
def Sasa_calc(xyz,atom_radii,probe_radius=0.14,n_sphere_points=960):
    radii = np.array(atom_radii, np.float32) + probe_radius
    dim1 = xyz.shape[1]
    out = np.zeros((xyz.shape[0], dim1), dtype=np.float32)
    atom_mapping = np.arange(dim1, dtype=np.int32)
    _geometry._sasa(xyz, radii, int(n_sphere_points), atom_mapping, out)
    return np.round(out.sum(axis=1),3)*100 # in unit of Angstrom^2


# For BSA and Score visualization
def plot_BSA_Score(df_data, xcol='BSA$_{AL}$', ycol='BSA$_{BL}$', min_col='BSA$_{min}$',score = 'Docking Score', plot_margin = 20,
                   axs = None, legend = 'Ligand ID',show_legend=True):
    try:
        axs[0]
    except TypeError as err:
        print(err, '(No axs provided)')
        return 
    # df_data is the pd.DataFrame of all pockets/ligand
    # Plot the BSA to protein A and to protein B, colored by score in the first plot
    sns.scatterplot(data=df_data,x=xcol,y=ycol,hue=score,ax=axs[0],palette="flare")
    max_scale = np.max([df_data[xcol],df_data[ycol]])+plot_margin
    min_scale = np.min([df_data[xcol],df_data[ycol]])-plot_margin
    axs[0].grid()
    axs[0].set_xlim(min_scale,max_scale)
    axs[0].set_ylim(min_scale,max_scale)

    # Plot the BSA_min to the Score and show legend
    df_data_interface = df_data[df_data[min_col] > 20 ]
    
    for index, row in df_data_interface.iterrows():
        if legend == 'Ligand ID': # For ligand visualization
            axs[1].scatter(x=row[min_col],y=row[score],label=str(int(row.Ligand)))
        elif legend == 'Pocket ID': # For pocket visualization
            axs[1].scatter(x=row[min_col],y=row[score],label=str(int(row.Pocket)))

    if show_legend :
        axs[1].legend(title=legend,ncol=2)
    axs[1].set_xlabel(min_col)
    axs[1].set_ylabel(score)
    axs[1].grid()
    plt.tight_layout()
    
    
    

# For structure visualization
## Pocket 
def pocket_visualize(view,pocket_id = 1):
    view.clear()
    #view.add_line(selection='protein')
    view.add_cartoon(selection=":A",opacity=1,color="blue")
    view.add_cartoon(selection=":B",opacity=1,color="red")
    view.add_spacefill(selection=''.join([':C and ',str(pocket_id)]))
## Ligand
def show_ligand(view,S_path):
    view.clear()
    view.add_line(selection='protein')
    view.add_cartoon(selection=":A",opacity=1,color="blue")
    view.add_cartoon(selection=":B",opacity=1,color="red")
    S = md.load(S_path)
    view.add_component(S)

# Find box parameter with the pocket location
def pocket2boxsize(xyz):
    xyz_max = np.max(xyz,axis=1)
    xyz_min = np.min(xyz,axis=1)
    center = (xyz_max + xyz_min)/2
    size = (xyz_max - xyz_min)
    return np.round(center*10,2), np.round(size*10,2)


# Read docking output
pocket_pattern = 'HEADER 0  - Pocket Score\s+: (?P<Pocketscore>-?\d+.\d+)\nHEADER 1  - Drug Score\s+: (?P<Drugscore>-?\d+.\d+)'
vina_pattern   = r'REMARK VINA RESULT:\s+(?P<Dockingscore>-?\d+.\d+)' 
gpu_pattern    = r'DOCKED: USER    Estimated Free Energy of Binding    \=\s?\s?(?P<Dockingscore>.*) kcal/mol'
gbsa_pattern   = r'DELTA TOTAL\s+(?P<Average>-?\d+.\d+)\s+(?P<Std>\d+.\d+)\s+(?P<SEM>\d+.\d+)' 
decomp_pattern = '\\S{3,}\\s+(?P<idx>\\d+),[RL]\\s\\S{3}\\s+\\d+,(-?\\d+.\\d+e?-?[0-9]*,){15}(?P<energy>-?\\d+.\\d+e?-?[0-9]*)'

def read_output(path,pattern=pocket_pattern):
    with open(path,'r') as file:
        content = file.read()
        matches = re.findall(pattern,content,re.MULTILINE)
        if isinstance(matches[0], tuple):
            return [float(match) for match in matches[0]]
        else:
            return [float(match) for match in matches]

# MD simulation analysis
class ABL:
    def __init__(self,sim_path,name,runs=5):
        ## Initialize
        self.runs = runs
        self.sim_path = sim_path
        self.name = name
        self.prm = ''.join([sim_path,name,'/',name,'.prmtop'])

    def analysis(self,A_resid,B_resid,A_intres,B_intres):
        ## Setup arguments
        Ares_arg = ''.join([':',str(min(A_resid)),'-',str(max(A_resid))])
        Aintres_arg = ':'+','.join([str(resid) for resid in A_intres])
        Bintres_arg = ':'+','.join([str(resid) for resid in B_intres])
        Bres_arg = ''.join([':',str(min(B_resid)),'-',str(max(B_resid))])
        Lres_arg = ''.join([':',str(max(B_resid)+1)])
        ## Allocate memory
        self.hbonds  = [list] * self.runs
        self.A_rmsd  = [list] * self.runs
        self.B_rmsd  = [list] * self.runs
        self.L_rmsd  = [list] * self.runs
        self.A_irmsd = [list] * self.runs
        self.B_irmsd = [list] * self.runs
        self.A_rmsf  = [list] * self.runs
        self.B_rmsf  = [list] * self.runs
        self.L_rmsf  = [list] * self.runs     
        ## Analysis
        for run in range(1,self.runs+1):
            print("Analyzing traj",run,' of ',self.name,'...')
            # Read trajectory
            nc = ''.join([self.sim_path,self.name,'/',self.name,'_prod',str(run),'.nc'])
            traj = pt.load(nc,top = self.prm)
            traj = pt.superpose(traj)
            traj.strip(':Na+,Cl-,WAT')
            # Append data
            self.hbonds[run-1]  = [ name for name in pt.hbond(traj).donor_acceptor if str(max(B_resid)+1) in name]
            self.A_rmsd[run-1]  = pt.rmsd(traj,Ares_arg)
            self.B_rmsd[run-1]  = pt.rmsd(traj,Bres_arg)
            self.L_rmsd[run-1]  = pt.rmsd(traj,Lres_arg,nofit=1)
            self.A_irmsd[run-1] = pt.rmsd(traj,Aintres_arg)
            self.B_irmsd[run-1] = pt.rmsd(traj,Bintres_arg)
            self.A_rmsf[run-1]  = pt.rmsf(traj,Ares_arg,options='byres')
            self.B_rmsf[run-1]  = pt.rmsf(traj,Bres_arg,options='byres')
            self.L_rmsf[run-1]  = pt.rmsf(traj,Lres_arg,options='byres')
    
    ## Show RMSF plot in a axs = 3*N plt.subplots
    def plot_RMSD(self,axs=plt):
        RMSD = [self.A_rmsd, self.B_rmsd, self.L_rmsd]
        Titles = ['Protein A', 'Protein B', 'Ligand' ]
        for idx, ax in enumerate(axs.flat):
            for run in range(1,self.runs+1):
                ax.plot(RMSD[idx][run-1],lw=3,label=''.join(['Sim',str(run)]))
            ax.grid()
            ax.set_title(''.join([self.name,' (',Titles[idx],')']))
            ax.set_xlabel('frame')
            ax.set_ylabel('RMSD ($\AA$)')
        axs[0].legend()
    
    ## Show RMSF plot in a axs = 2*N plt.subplots
    def plot_iRMSD(self,axs=plt):
        RMSD = [self.A_irmsd, self.B_irmsd]
        Titles = ['Protein A', 'Protein B']        
        for idx, ax in enumerate(axs.flat):
            for run in range(1,self.runs+1):
                ax.plot(RMSD[idx][run-1],lw=3,label=''.join(['Sim',str(run)]))
            ax.grid()
            ax.set_title(''.join([self.name,' (',Titles[idx],')']))
            ax.set_xlabel('frame')
            ax.set_ylabel('iRMSD ($\AA$)')
        axs[0].legend()
        
    ## Show RMSF plot in a axs = 2*N plt.subplots      
    def plot_RMSF(self,axs=plt,A_start=1, B_start=1,show_interface=True,A_intres=[],B_intres=[]):       
        Intres = [A_intres, B_intres]
        RMSF = [self.A_rmsf, self.B_rmsf]
        Titles = ['Protein A', 'Protein B']
        colors = ['blue','red']
        Resids = [np.arange(A_start,len(self.A_rmsf[0])+1),np.arange(B_start,len(self.B_rmsf[0])+1)]

        for idx, ax in enumerate(axs.flat):
            avg = np.average(np.array(RMSF[idx])[:,:,1],axis=0)
            std = np.std(np.array(RMSF[idx])[:,:,1],axis=0)
            ax.plot(np.array(RMSF[idx])[0,:,0],avg,color=colors[idx])
            ax.fill_between(np.array(RMSF[idx])[0,:,0],y1=avg-std,y2=avg+std,color=colors[idx],alpha=0.2)
            if show_interface:
                for res in Intres[idx]:
                    ax.plot(res,0.5,'ko',ms=3)
            ax.grid()
            ax.set_title(''.join([self.name,' (',Titles[idx],')']))
            ax.set_xlabel('Residue ID')
            ax.set_ylabel('RMSF ($\AA$)')

    ## Calculate Hbonds
    def get_hbonds(self):
        all_hbonds=[]
        for hbonds in test.hbonds:
            for hbond in hbonds:
                acceptor = hbond.split('-')[0].split('_')[0]
                donor    = hbond.split('-')[1].split('_')[0]
                all_hbonds.append(acceptor)
                all_hbonds.append(donor)
        return np.unique(all_hbonds).tolist()
                
    def visualize(self,run=1,hightlight=[]):
        nc = ''.join([self.sim_path,self.name,'/',self.name,'_prod',str(run),'.nc'])
        traj = pt.load(nc,top = self.prm)
        traj.strip(':Na+,Cl-,WAT')
        traj = pt.superpose(traj)
        view = nv.show_pytraj(traj)
        view.clear()
        #view.add_line(selection='protein')
        view.add_cartoon(selection=":A",opacity=1,color="blue")
        view.add_cartoon(selection=":B",opacity=1,color="red")
        view.add_ball_and_stick(selection=":C and not hydrogen",opacity=1,aspectRatio=2,radiusSegments=10)
        for res in hightlight:
            view.add_ball_and_stick(selection=str(res),opacity=1,aspectRatio=1)
        return view
    
    ## Read MMGBSA output, including total ddG and energy decompositions
    def read_mmgbsa(self,gbsa_dir='../Step6_MMGBSA/'):
        
        # initialize
        self.AL_gbsa  = np.zeros(self.runs)
        self.BL_gbsa  = np.zeros(self.runs)
        
        self.AL_gbsa_decomp  = [np.array] * self.runs
        self.BL_gbsa_decomp  = [np.array] * self.runs
        
        # read value
        for run in range(1,self.runs+1):
            # A
            AL_gbsa_path   = ''.join([gbsa_dir,self.name,'/AL_output_',str(run),'.dat'])
            AL_decomp_path = ''.join([gbsa_dir,self.name,'/AL_decomp_',str(run),'.csv'])
            self.AL_gbsa[run-1]        = float(read_output(AL_gbsa_path,gbsa_pattern)[0][0])
            self.AL_gbsa_decomp[run-1] = pd.read_csv(AL_decomp_path,skiprows=8,header=None)[17].values[:-1]
            
            # B
            BL_gbsa_path = ''.join([gbsa_dir,self.name,'/BL_output_',str(run),'.dat'])
            BL_decomp_path = ''.join([gbsa_dir,self.name,'/BL_decomp_',str(run),'.csv'])
            self.BL_gbsa[run-1] = float(read_output(BL_gbsa_path,gbsa_pattern)[0][0])
            self.BL_gbsa_decomp[run-1] = pd.read_csv(BL_decomp_path,skiprows=8,header=None)[17].values[:-1]
            
    def visualize_favorable_contact(self,run=1,show_favorable = 1, highlight = [], threshold=[-1,1],fav_color = ['white','black']):
        nc = ''.join([self.sim_path,self.name,'/',self.name,'_prod',str(run),'.nc'])
        traj = pt.load(nc,top = self.prm)
        traj.strip(':Na+,Cl-,WAT')
        traj = pt.superpose(traj)
        view = nv.show_pytraj(traj)
        view.clear()        
        view.add_cartoon(selection=":A",opacity=1,color="blue")
        view.add_cartoon(selection=":B",opacity=1,color="red")
        view.add_ball_and_stick(selection=":C and not hydrogen",opacity=1,aspectRatio=2,radiusSegments=10)
        
        if show_favorable:
            #try:
            gbsa_decomp = [self.AL_gbsa_decomp, self.BL_gbsa_decomp]
            resid_offset = [1,1+len(self.AL_gbsa_decomp[0])]
            for idx, decomp in enumerate(gbsa_decomp):
                avg = np.average(decomp,axis=0)
                for fav_res in np.where(avg < threshold[0])[0]:
                    view.add_ball_and_stick(selection=str(fav_res+resid_offset[idx]),color=fav_color[0])
                for unfav_res in np.where(avg > threshold[1])[0]:
                    view.add_ball_and_stick(selection=str(unfav_res+resid_offset[idx]),color=fav_color[1])
            #except:
            #    print("Can't read the deomposition file, please do RLS.read_mmgbsa() first!")
            #    return 
        for res in highlight:
            view.add_ball_and_stick(selection=str(res),opacity=1,aspectRatio=1)
        
        return view
         

                  
        

        
      

   
