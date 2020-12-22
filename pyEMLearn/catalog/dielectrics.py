from numpy import array
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

Pcd2 = CauchyIndex([1.10946563843109,0.0073377449228674])\
        + LorentzIndex2(
                A = [0.172439520432108,0.0779888680018249,0.0233478966236617,0.0185378206694171,0.0106394246243079],
                E = [1.91598020861984,2.08952480489883,2.25404183731625,2.82841891009979,2.993],
                G = [0.1301859024937,0.0870773365497822,0.0998623211390882,0.0301332046211885,0.0527937552588704]
            )
Pcd2.name = "TIPS-Pc:PMMA (d=2nm)"

TIPSTc = CauchyIndex([1.50696993820317,0.0214751876991838])\
        + LorentzIndex2(
                A = [0.99322443495843,0.596245237614556,0.208665420057777,0.0492614296902469,0.134621751497327],
                E = [2.30873324388158,2.47936620693357,2.65301072409464,2.8220070395855,3.08136734362406],
                G = [0.0585569687127509,0.0818609116329319,0.103736543502073,0.101207044916861,0.287520008683127]
            )
TIPSTc.name = "TIPS-Tc (400 to 1300 nm)"

# extracted from pmma_lp400.tsv
PMMA = CauchyIndex(
        [1.46402664, 0.00552999],
        name="PMMA LP400"
        )

# extracted from pva_lp350.tsv
PVA = CauchyIndex(
        [1.48624997, 0.00778431],
        name="PVA LP350"
        )

# diF TES-ADT Aggregates spun at 2k from Hexane
res = array([
[350,1.73091228,0.00452631254],
[352,1.72937685,0.00476041186],
[354,1.7278474,0.00500752537],
[356,1.72632372,0.00526815861],
[358,1.72480564,0.00554276524],
[360,1.72329303,0.00583172843],
[362,1.72178581,0.00613534019],
[364,1.72028389,0.0064537792],
[366,1.71878719,0.00678708823],
[368,1.71729556,0.00713515279],
[370,1.71580879,0.00749768302],
[372,1.71432648,0.00787420194],
[374,1.712848,0.00826404365],
[376,1.71137244,0.00866636585],
[378,1.70989845,0.00908018173],
[380,1.70842416,0.00950441657],
[382,1.70694707,0.0099379944],
[384,1.70546393,0.0103799604],
[386,1.70397059,0.010829644],
[388,1.70246191,0.0112868694],
[390,1.70093163,0.0117522209],
[392,1.69937226,0.0122273787],
[394,1.69777491,0.0127155512],
[396,1.69612921,0.0132220549],
[398,1.69442303,0.013755139],
[400,1.69264231,0.0143272362],
[402,1.69077081,0.0149569761],
[404,1.68879004,0.0156725999],
[406,1.68668014,0.0165179782],
[408,1.68442356,0.0175634476],
[410,1.68201787,0.0189250705],
[412,1.67951578,0.0207950432],
[414,1.67714194,0.0234647299],
[416,1.6755687,0.0271960782],
[418,1.67611453,0.0314446153],
[420,1.67916951,0.03371878],
[422,1.68171546,0.0323805224],
[424,1.68143463,0.0297729132],
[426,1.67919978,0.0279712848],
[428,1.6761926,0.0272554826],
[430,1.67295097,0.0273698555],
[432,1.66965226,0.0280855013],
[434,1.66634259,0.0292693352],
[436,1.66302888,0.0308576208],
[438,1.65971138,0.032829749],
[440,1.65639568,0.0351921019],
[442,1.65309818,0.0379686709],
[444,1.64985021,0.0411949819],
[446,1.6467021,0.0449129238],
[448,1.64372784,0.0491646593],
[450,1.64102958,0.0539840007],
[452,1.63874069,0.0593837853],
[454,1.6370243,0.0653383771],
[456,1.63606282,0.0717621451],
[458,1.63603335,0.0784883871],
[460,1.63706533,0.0852586147],
[462,1.63918515,0.0917368905],
[464,1.64226495,0.0975619586],
[466,1.64600496,0.10243487],
[468,1.6499746,0.106214563],
[470,1.65371019,0.108977027],
[472,1.6568325,0.111006017],
[474,1.65913465,0.112720655],
[476,1.66061303,0.114576901],
[478,1.66144767,0.116981668],
[480,1.66195792,0.120237089],
[482,1.66255533,0.124510366],
[484,1.66370155,0.129815489],
[486,1.66586561,0.135996789],
[488,1.66946905,0.142715014],
[490,1.6748108,0.149449401],
[492,1.68197722,0.155537955],
[494,1.69076399,0.160274267],
[496,1.7006534,0.163055425],
[498,1.71088257,0.163538607],
[500,1.72059836,0.161740472],
[502,1.72904368,0.158031561],
[504,1.73570075,0.153031903],
[506,1.74034593,0.147462622],
[508,1.74302201,0.142016234],
[510,1.74396715,0.137280152],
[512,1.74354223,0.133713285],
[514,1.74218133,0.131655574],
[516,1.7403725,0.131346517],
[518,1.73866369,0.132932555],
[520,1.73768142,0.136447906],
[522,1.73814175,0.14175798],
[524,1.74082095,0.148464563],
[526,1.74644308,0.155797821],
[528,1.75545411,0.16256716],
[530,1.76772243,0.167285189],
[532,1.7823282,0.168539434],
[534,1.79766546,0.165504849],
[536,1.81191754,0.158288176],
[538,1.82365283,0.147845647],
[540,1.83217491,0.135544154],
[542,1.83749262,0.122682274],
[544,1.84005926,0.110214053],
[546,1.84049444,0.0987000445],
[548,1.83939952,0.0883863664],
[550,1.83727397,0.0793148911],
[552,1.83449737,0.0714153054],
[554,1.83134241,0.0645669629],
[556,1.82799761,0.0586350436],
[558,1.82458941,0.053489291],
[560,1.82120026,0.0490123288],
[562,1.81788222,0.0451022912],
[564,1.81466674,0.0416726159],
[566,1.8115714,0.0386505911],
[568,1.80860466,0.0359754909],
[570,1.80576897,0.0335967068],
[572,1.80306301,0.0314720494],
[574,1.80048306,0.0295662838],
[576,1.79802401,0.0278498995],
[578,1.79568001,0.0262980934],
[580,1.79344487,0.0248899382],
[582,1.79131233,0.0236077037],
[584,1.7892763,0.0224363047],
[586,1.78733085,0.0213628522],
[588,1.78547041,0.0203762866],
[590,1.78368967,0.0194670793],
[592,1.78198369,0.0186269872],
[594,1.78034785,0.0178488514],
[596,1.77877784,0.017126431],
[598,1.77726967,0.0164542652],
[600,1.77581962,0.0158275594],
[602,1.77442427,0.0152420893],
[604,1.77308041,0.0146941208],
[606,1.7717851,0.0141803433],
[608,1.77053558,0.0136978125],
[610,1.76932932,0.0132439025],
[612,1.76816395,0.0128162648],
[614,1.76703728,0.0124127934],
[616,1.76594725,0.0120315952],
[618,1.76489198,0.0116709643],
[620,1.76386968,0.0113293595],
[622,1.7628787,0.0110053859],
[624,1.76191752,0.0106977778],
[626,1.76098467,0.0104053845],
[628,1.76007881,0.0101271581],
[630,1.75919869,0.0098621418],
[632,1.75834312,0.00960946108],
[634,1.75751099,0.00936831481],
[636,1.75670125,0.00913796792],
[638,1.75591294,0.00891774491],
[640,1.75514513,0.00870702396],
[642,1.75439694,0.00850523185],
[644,1.75366758,0.00831183935],
[646,1.75295626,0.00812635715],
[648,1.75226225,0.00794833227],
[650,1.75158488,0.00777734473],
[660,1.74842591,0.00701523911],
[680,1.74304879,0.00584536092],
[700,1.73862188,0.00499444677],
[720,1.7348976,0.00435123489],
[740,1.73171091,0.00385004302],
[760,1.72894693,0.00344980209],
[780,1.72652281,0.00312362397],
[800,1.72437705,0.00285322867],
[820,1.72246279,0.00262579111],
[840,1.72074361,0.00243206983],
[860,1.71919067,0.00226525105],
[880,1.71778078,0.00212020973],
[900,1.71649503,0.00199302299],
[920,1.71531781,0.00188064121],
[940,1.71423609,0.00178066073],
[960,1.71323892,0.00169116337],
[980,1.71231698,0.0016106012],
[1000,1.71146232,0.00153771234],
[1020,1.71066806,0.00147145855],
[1040,1.70992827,0.00141097832],
[1060,1.70923775,0.00135555127],
[1080,1.70859197,0.00130457063],
[1100,1.70798692,0.00125752198],
[1120,1.70741906,0.00121396649],
[1140,1.70688525,0.00117352767],
[1160,1.7063827,0.00113588077],
[1180,1.70590889,0.00110074428],
[1200,1.70546159,0.00106787306],
[1220,1.70503877,0.00103705264],
[1240,1.70463861,0.00100809467],
[1260,1.70425946,0.000980833066],
[1280,1.70389981,0.00095512083],
[1300,1.70355832,0.000930827426]
])
res[:,0] = res[:,0]*1e-3
ADTagg = DiscreteIndex(
         res[:,:2],
         res[:,::2],
         name = "diF-TES-ADT aggregate LP400"
     )
