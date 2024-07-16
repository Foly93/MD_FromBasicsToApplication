from modeller import *

env = Environ()
aln = Alignment(env)
mdl = Model(env, file='hmx')
# Known structure
aln.append_model(mdl, align_codes='hmx_structure', atom_files='hmx.pdb')
# Unknown structure
aln.append(file='hmx.fasta', align_codes='hmx')
aln.salign()
aln.write(file='hmx.ali', alignment_format='PIR')
aln.write(file='hmx.pap', alignment_format='PAP')
