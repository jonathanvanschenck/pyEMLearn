from numpy import genfromtxt
from os import path

from pyEMLearn.materials import SellmeierIndex, SellmeierIndexSquare,\
                                    DiscreteIndex, CauchyIndex, Material, \
                                    LorentzIndex2


Air = SellmeierIndex(
          [0.05792105 ,0.0016791], #1/um^2
          [238.01,57.36], #1/um^2
          name = "Air"
      )

BK7 = SellmeierIndexSquare(
          [1.03961212 ,0.231792344,1.01046945],
          [0.00600069867,0.0200179144,103.560653], #um^2
          name = "BK7"
       )

res = genfromtxt(path.join("pyEMLearn","catalog","tsv","PcPMMA_d2nm.tsv"),delimiter="\t",skip_header=3)
res[:,0] = res[:,0]*1e-3
Pcd2_tsv = DiscreteIndex(
         res[:,:2],
         res[:,::2],
         name = "TIPS-Pc:PMMA (d=2nm) tsv"
     )

Pcd2 = CauchyIndex([1.10946563843109,0.0073377449228674])\
        + LorentzIndex2(
                A = [0.172439520432108,0.0779888680018249,0.0233478966236617,0.0185378206694171,0.0106394246243079],
                E = [1.91598020861984,2.08952480489883,2.25404183731625,2.82841891009979,2.993],
                G = [0.1301859024937,0.0870773365497822,0.0998623211390882,0.0301332046211885,0.0527937552588704]
            )
Pcd2.name = "TIPS-Pc:PMMA (d=2nm)"


res = genfromtxt(path.join("pyEMLearn","catalog","tsv","teshex_unknown2_2k.tsv"),delimiter="\t",skip_header=3)
res[:,0] = res[:,0]*1e-3
ADTagg = DiscreteIndex(
         res[:,:2],
         res[:,::2],
         name = "diF-TES-ADT aggregate LP400"
     )

#res = genfromtxt(path.join("pyEMLearn","catalog","tsv","adt_hex_2k.tsv"),delimiter="\t",skip_header=3)
#res[:,0] = res[:,0]*1e-3
#ADT = DiscreteIndex(
#         res[:,:2],
#         res[:,::2],
#         name = "diF-TES-ADT"
#     )
#
#ADTIR = CauchyIndex(
#        [1.70126573e+00, 4.91770179e-04, 9.05393327e-03],
#        name = "diF-TES-ADT IR"
#    )


#res = genfromtxt(path.join("pyEMLearn","catalog","tsv","pmma_lp400.tsv"),delimiter="\t",skip_header=3)
#res[:,0] = res[:,0]*1e-3
#PMMA2 = DiscreteIndex(
#         res[:,:2],
#         res[:,::2],
#         name = "PMMA LP400"
#     )
# extracted from pmma_lp400.tsv
PMMA = CauchyIndex(
        [1.46402664, 0.00552999],
        name="PMMA LP400"
        )

#res = genfromtxt(path.join("pyEMLearn","catalog","tsv","pva_lp350.tsv"),delimiter="\t",skip_header=3)
#res[:,0] = res[:,0]*1e-3
#PVA2 = DiscreteIndex(
#         res[:,:2],
#         res[:,::2],
#         name = "PVA LP350"
#     )
# extracted from pva_lp350.tsv
PVA = CauchyIndex(
        [1.48624997, 0.00778431],
        name="PVA LP350"
        )

