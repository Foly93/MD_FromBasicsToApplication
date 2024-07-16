reps=2

### change back to 1 later on
for rep in $(seq 2 $reps)
do
	MMPBSA.py -i ../gb.in -y AL_$rep.nc -cp AL.prmtop -sp AL.prmtop -rp A.prmtop -lp L.prmtop -o AL_output_$rep.dat -do AL_decomp_$rep.csv
	MMPBSA.py -i ../gb.in -y BL_$rep.nc -cp BL.prmtop -sp BL.prmtop -rp B.prmtop -lp L.prmtop -o BL_output_$rep.dat -do BL_decomp_$rep.csv
done


