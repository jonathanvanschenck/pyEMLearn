from numpy import genfromtxt
from os import path

from pyEMLearn2.materials import DiscreteIndex


res = genfromtxt(path.join("pyEMLearn2","catalog","tsv","ag_tf_lp300.tsv"),delimiter="\t",skip_header=3)
res[:,0] = res[:,0]*1e-3
Ag = DiscreteIndex(
         res[:,:2],
         res[:,::2],
         name = "Ag LP300"
     )

# This is the same as Ag2, but with some extra lorentzians
res = genfromtxt(path.join("pyEMLearn2","catalog","tsv","ag_tf_lp300_complex.tsv"),delimiter="\t",skip_header=3)
res[:,0] = res[:,0]*1e-3
Ag2 = DiscreteIndex(
         res[:,:2],
         res[:,::2],
         name = "Ag LP300 complex"
     )