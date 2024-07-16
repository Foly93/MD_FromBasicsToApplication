from modeller import *
from modeller.automodel import *

env = Environ()
a = automodel(env, alnfile='hmx.ali',
              knowns='hmx_structure', sequence='hmx',
              assess_methods=(assess.DOPE))
a.starting_model = 1
a.ending_model = 5
a.make()
