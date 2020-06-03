from numpy import array
from os import path

from pyEMLearn.materials import DiscreteIndex

# Thin film (d ~ 50nm) Ag, thermally evaporated onto glass
res = array([
[300,1.3655244,0.358952647],
[305,1.21665976,0.35211228],
[310,1.06019586,0.35593538],
[315,0.8924704,0.375147988],
[320,0.711163177,0.420523172],
[325,0.52552472,0.511550176],
[330,0.371789855,0.653941741],
[335,0.272514699,0.811596378],
[340,0.211351306,0.957388585],
[345,0.171197586,1.08737956],
[350,0.143042567,1.20390172],
[355,0.122305031,1.30967389],
[360,0.106480145,1.40686839],
[365,0.0940884705,1.49712384],
[370,0.0841977113,1.5816765],
[375,0.0761889006,1.6614726],
[380,0.0696334954,1.73724904],
[385,0.0642245123,1.80958943],
[390,0.0597357966,1.87896314],
[395,0.0559968015,1.94575285],
[400,0.0528763519,2.01027444],
[405,0.0502718406,2.07279143],
[410,0.0481018291,2.13352582],
[415,0.0463008548,2.19266625],
[420,0.0448157058,2.25037424],
[425,0.0436026985,2.30678907],
[430,0.0426256519,2.36203159],
[435,0.0418543609,2.4162073],
[440,0.0412634269,2.46940873],
[445,0.0408313532,2.52171748],
[450,0.0405398367,2.57320577],
[455,0.0403732083,2.62393785],
[460,0.0403179861,2.67397103],
[465,0.040362516,2.72335663],
[470,0.0404966808,2.77214073],
[475,0.0407116622,2.82036486],
[480,0.0409997458,2.86806648],
[485,0.0413541597,2.9152795],
[490,0.0417689405,2.96203459],
[495,0.0422388217,3.00835956],
[500,0.04275914,3.05427956],
[505,0.0433257584,3.09981717],
[510,0.0439350024,3.14499239],
[515,0.0445836177,3.18982172],
[520,0.0452622184,3.23431392],
[525,0.0459523675,3.27849415],
[530,0.0466528372,3.32238784],
[535,0.047363367,3.36601393],
[540,0.0480837358,3.40938876],
[545,0.0488137495,3.45252688],
[550,0.0495532353,3.49544146],
[555,0.0503020378,3.53814448],
[560,0.0510600161,3.58064699],
[565,0.0518270425,3.62295914],
[570,0.0526030001,3.66509034],
[575,0.0533877821,3.70704931],
[580,0.0541812902,3.74884418],
[585,0.0549834341,3.79048251],
[590,0.0557941304,3.83197138],
[595,0.0566133021,3.87331739],
[600,0.0574408778,3.91452674],
[605,0.0582767912,3.95560527],
[610,0.0591209808,3.99655842],
[615,0.0599733893,4.03739136],
[620,0.0608339633,4.07810892],
[625,0.061702653,4.11871569],
[630,0.0625794118,4.15921597],
[635,0.063464196,4.19961386],
[640,0.0643569648,4.23991322],
[645,0.0652576799,4.28011771],
[650,0.0661663054,4.32023078],
[655,0.0670828073,4.36025575],
[660,0.0680071539,4.40019572],
[665,0.0689393151,4.44005368],
[670,0.0698792628,4.47983244],
[675,0.0708269701,4.5195347],
[680,0.0717824119,4.55916301],
[685,0.0727455644,4.59871981],
[690,0.0737164049,4.63820745],
[695,0.0746949122,4.67762813],
[700,0.075681066,4.71698398],
[705,0.0766748471,4.75627703],
[710,0.0776762372,4.79550922],
[715,0.0786852191,4.8346824],
[720,0.0797017762,4.87379836],
[725,0.080725893,4.91285879],
[730,0.0817575545,4.95186532],
[735,0.0827967464,4.99081952],
[740,0.0838434553,5.02972289],
[745,0.0848976683,5.06857687],
[750,0.085959373,5.10738284],
[755,0.0870285576,5.14614214],
[760,0.088105211,5.18485603],
[765,0.0891893222,5.22352574],
[770,0.0902808812,5.26215247],
[775,0.0913798779,5.30073734],
[780,0.0924863029,5.33928145],
[785,0.0936001473,5.37778585],
[790,0.0947214024,5.41625157],
[795,0.0958500597,5.45467957],
[800,0.0969861114,5.49307082],
[805,0.0981295497,5.53142621],
[810,0.0992803674,5.56974663],
[815,0.100438557,5.60803292],
[820,0.101604112,5.64628591],
[825,0.102777027,5.6845064],
[830,0.103957293,5.72269514],
[835,0.105144906,5.76085287],
[840,0.10633986,5.79898032],
[845,0.107542148,5.83707818],
[850,0.108751766,5.87514711],
[855,0.109968709,5.91318776],
[860,0.11119297,5.95120077],
[865,0.112424546,5.98918674],
[870,0.113663431,6.02714626],
[875,0.114909621,6.0650799],
[880,0.116163111,6.10298821],
[885,0.117423898,6.14087174],
[890,0.118691977,6.178731],
[895,0.119967344,6.2165665],
[900,0.121249995,6.25437873],
[905,0.122539926,6.29216816],
[910,0.123837134,6.32993526],
[915,0.125141615,6.36768048],
[920,0.126453366,6.40540424],
[925,0.127772384,6.44310697],
[930,0.129098664,6.48078909],
[935,0.130432204,6.518451],
[940,0.131773002,6.55609308],
[945,0.133121053,6.5937157],
[950,0.134476355,6.63131925],
[955,0.135838906,6.66890407],
[960,0.137208702,6.70647052],
[965,0.13858574,6.74401892],
[970,0.139970019,6.78154962],
[975,0.141361536,6.81906293],
[980,0.142760288,6.85655917],
[985,0.144166273,6.89403863],
[990,0.145579488,6.93150162],
[995,0.146999931,6.96894843],
[1000,0.1484276,7.00637932],
[1005,0.149862493,7.04379459],
[1010,0.151304607,7.0811945],
[1015,0.152753941,7.1185793],
[1020,0.154210492,7.15594926],
[1025,0.155674259,7.19330461],
[1030,0.157145239,7.23064561],
[1035,0.158623431,7.26797248],
[1040,0.160108832,7.30528547],
[1045,0.161601442,7.34258478],
[1050,0.163101257,7.37987065],
[1055,0.164608277,7.41714328],
[1060,0.166122499,7.45440289],
[1065,0.167643922,7.49164968],
[1070,0.169172544,7.52888384],
[1075,0.170708364,7.56610559],
[1080,0.17225138,7.60331509],
[1085,0.173801591,7.64051255],
[1090,0.175358994,7.67769813],
[1095,0.176923589,7.71487203],
[1100,0.178495373,7.75203441],
[1105,0.180074346,7.78918544],
[1110,0.181660506,7.82632529],
[1115,0.183253852,7.86345412],
[1120,0.184854381,7.9005721],
[1125,0.186462094,7.93767936],
[1130,0.188076988,7.97477608],
[1135,0.189699062,8.01186239],
[1140,0.191328315,8.04893844],
[1145,0.192964745,8.08600438],
[1150,0.194608351,8.12306033],
[1155,0.196259133,8.16010645],
[1160,0.197917088,8.19714286],
[1165,0.199582216,8.23416969],
[1170,0.201254515,8.27118708],
[1175,0.202933984,8.30819514],
[1180,0.204620622,8.345194],
[1185,0.206314427,8.38218378],
[1190,0.2080154,8.41916459],
[1195,0.209723538,8.45613656],
[1200,0.21143884,8.4930998],
[1205,0.213161305,8.53005441],
[1210,0.214890933,8.5670005],
[1215,0.216627722,8.60393819],
[1220,0.21837167,8.64086757],
[1225,0.220122778,8.67778875],
[1230,0.221881043,8.71470183],
[1235,0.223646465,8.75160691],
[1240,0.225419043,8.78850408],
[1245,0.227198776,8.82539343],
[1250,0.228985663,8.86227506],
[1255,0.230779702,8.89914907],
[1260,0.232580893,8.93601553],
[1265,0.234389235,8.97287454],
[1270,0.236204726,9.00972619],
[1275,0.238027366,9.04657055],
[1280,0.239857155,9.08340771],
[1285,0.24169409,9.12023775],
[1290,0.24353817,9.15706076],
[1295,0.245389396,9.1938768],
[1300,0.247247766,9.23068595]
])
res[:,0] = res[:,0]*1e-3
Ag = DiscreteIndex(
         res[:,:2],
         res[:,::2],
         name = "Ag LP300"
     )

# This is the same as Ag2, but with some extra lorentzians
res = array([
[300,1.51294414,0.713309223],
[305,1.41012588,0.61448513],
[310,1.27520573,0.519335392],
[315,1.09913986,0.434934378],
[320,0.866201301,0.377236798],
[325,0.552941829,0.401128302],
[330,0.273661263,0.637690802],
[335,0.209272373,0.890371495],
[340,0.199091875,1.04958484],
[345,0.189471549,1.16454994],
[350,0.172681553,1.25910109],
[355,0.150475518,1.34699325],
[360,0.128213079,1.43392835],
[365,0.10995663,1.51991113],
[370,0.096750997,1.60301629],
[375,0.0876381117,1.68189445],
[380,0.0811703743,1.75632111],
[385,0.0762231015,1.82677084],
[390,0.0721403917,1.89391937],
[395,0.0686018183,1.95837791],
[400,0.0654657545,2.02062037],
[405,0.0626692724,2.08099721],
[410,0.0601794979,2.139769],
[415,0.0579738606,2.19713488],
[420,0.0560330927,2.25325193],
[425,0.0543390349,2.30824756],
[430,0.0528740656,2.36222746],
[435,0.0516210385,2.41528102],
[440,0.0505633665,2.46748497],
[445,0.0496851329,2.5189061],
[450,0.0489711937,2.56960317],
[455,0.0484072564,2.6196284],
[460,0.0479799331,2.66902856],
[465,0.0476767685,2.71784586],
[470,0.0474862465,2.76611857],
[475,0.0473977774,2.8138816],
[480,0.0474016714,2.86116691],
[485,0.047489101,2.90800385],
[490,0.0476520558,2.95441946],
[495,0.0478832914,3.0004387],
[500,0.0481762769,3.04608464],
[505,0.0485251402,3.0913787],
[510,0.0489246141,3.13634071],
[515,0.0493699842,3.18098912],
[520,0.0498570386,3.22534106],
[525,0.05038202,3.26941252],
[530,0.0509415817,3.31321835],
[535,0.051532746,3.35677242],
[540,0.052152866,3.40008769],
[545,0.0527995914,3.44317624],
[550,0.053470836,3.48604936],
[555,0.0541647492,3.52871762],
[560,0.0548796904,3.57119092],
[565,0.0556142051,3.61347852],
[570,0.0563670042,3.65558911],
[575,0.0571369453,3.69753083],
[580,0.0579230157,3.73931133],
[585,0.0587243179,3.78093781],
[590,0.0595400561,3.822417],
[595,0.0603695244,3.86375529],
[600,0.0612120964,3.90495864],
[605,0.0620672161,3.9460327],
[610,0.0629343898,3.98698279],
[615,0.0638131783,4.02781394],
[620,0.0647031915,4.06853088],
[625,0.0656040816,4.1091381],
[630,0.0665155393,4.14963983],
[635,0.0674372884,4.1900401],
[640,0.0683690827,4.2303427],
[645,0.0693107022,4.27055123],
[650,0.0702619499,4.31066913],
[655,0.0712226499,4.35069962],
[660,0.0721926441,4.39064581],
[665,0.0731717908,4.43051062],
[670,0.0741599627,4.47029684],
[675,0.0751570453,4.51000712],
[680,0.0761629353,4.54964401],
[685,0.0771775397,4.58920989],
[690,0.0782007743,4.62870709],
[695,0.079232563,4.66813779],
[700,0.0802728369,4.70750409],
[705,0.0813215334,4.746808],
[710,0.0823785956,4.78605142],
[715,0.0834439718,4.8252362],
[720,0.0845176147,4.86436408],
[725,0.0855994812,4.90343675],
[730,0.086689532,4.94245581],
[735,0.0877877308,4.98142281],
[740,0.0888940446,5.02033924],
[745,0.0900084428,5.0592065],
[750,0.0911308973,5.09802598],
[755,0.0922613822,5.13679897],
[760,0.0933998737,5.17552674],
[765,0.0945463496,5.21421049],
[770,0.0957007893,5.2528514],
[775,0.0968631737,5.29145059],
[780,0.0980334852,5.33000914],
[785,0.0992117072,5.36852808],
[790,0.100397824,5.40700842],
[795,0.101591822,5.44545112],
[800,0.102793687,5.48385713],
[805,0.104003406,5.52222733],
[810,0.105220968,5.5605626],
[815,0.106446362,5.59886377],
[820,0.107679576,5.63713165],
[825,0.108920602,5.67536703],
[830,0.110169429,5.71357066],
[835,0.111426048,5.75174327],
[840,0.112690452,5.78988557],
[845,0.113962633,5.82799824],
[850,0.115242582,5.86608195],
[855,0.116530293,5.90413733],
[860,0.117825759,5.942165],
[865,0.119128973,5.98016557],
[870,0.120439929,6.01813961],
[875,0.121758621,6.0560877],
[880,0.123085043,6.09401038],
[885,0.124419191,6.13190817],
[890,0.125761058,6.1697816],
[895,0.12711064,6.20763117],
[900,0.128467931,6.24545735],
[905,0.129832928,6.28326062],
[910,0.131205626,6.32104144],
[915,0.132586021,6.35880025],
[920,0.133974108,6.39653747],
[925,0.135369883,6.43425354],
[930,0.136773343,6.47194885],
[935,0.138184484,6.5096238],
[940,0.139603303,6.54727878],
[945,0.141029795,6.58491417],
[950,0.142463958,6.62253031],
[955,0.143905788,6.66012758],
[960,0.145355282,6.69770631],
[965,0.146812438,6.73526684],
[970,0.148277252,6.7728095],
[975,0.149749721,6.8103346],
[980,0.151229842,6.84784245],
[985,0.152717613,6.88533336],
[990,0.154213032,6.92280763],
[995,0.155716095,6.96026553],
[1000,0.1572268,6.99770734],
[1005,0.158745144,7.03513335],
[1010,0.160271126,7.07254381],
[1015,0.161804743,7.10993898],
[1020,0.163345992,7.14731912],
[1025,0.164894872,7.18468447],
[1030,0.16645138,7.22203528],
[1035,0.168015515,7.25937177],
[1040,0.169587273,7.29669418],
[1045,0.171166653,7.33400273],
[1050,0.172753654,7.37129764],
[1055,0.174348272,7.40857911],
[1060,0.175950507,7.44584737],
[1065,0.177560356,7.48310261],
[1070,0.179177818,7.52034504],
[1075,0.18080289,7.55757483],
[1080,0.182435571,7.5947922],
[1085,0.18407586,7.63199731],
[1090,0.185723754,7.66919036],
[1095,0.187379252,7.70637151],
[1100,0.189042353,7.74354095],
[1105,0.190713054,7.78069883],
[1110,0.192391354,7.81784533],
[1115,0.194077251,7.85498062],
[1120,0.195770745,7.89210483],
[1125,0.197471833,7.92921814],
[1130,0.199180514,7.96632069],
[1135,0.200896786,8.00341264],
[1140,0.202620649,8.04049412],
[1145,0.2043521,8.07756528],
[1150,0.206091138,8.11462626],
[1155,0.207837763,8.15167719],
[1160,0.209591971,8.18871821],
[1165,0.211353763,8.22574945],
[1170,0.213123137,8.26277103],
[1175,0.214900091,8.29978308],
[1180,0.216684624,8.33678573],
[1185,0.218476736,8.3737791],
[1190,0.220276423,8.41076329],
[1195,0.222083687,8.44773843],
[1200,0.223898524,8.48470463],
[1205,0.225720934,8.52166201],
[1210,0.227550916,8.55861066],
[1215,0.229388468,8.5955507],
[1220,0.23123359,8.63248223],
[1225,0.23308628,8.66940535],
[1230,0.234946536,8.70632016],
[1235,0.236814359,8.74322677],
[1240,0.238689746,8.78012526],
[1245,0.240572696,8.81701574],
[1250,0.242463209,8.85389829],
[1255,0.244361283,8.89077301],
[1260,0.246266917,8.92763998],
[1265,0.24818011,8.9644993],
[1270,0.250100861,9.00135105],
[1275,0.252029169,9.03819531],
[1280,0.253965032,9.07503216],
[1285,0.25590845,9.1118617],
[1290,0.257859422,9.14868399],
[1295,0.259817946,9.18549912],
[1300,0.261784022,9.22230716]
])
res[:,0] = res[:,0]*1e-3
Ag2 = DiscreteIndex(
         res[:,:2],
         res[:,::2],
         name = "Ag LP300 complex"
     )
