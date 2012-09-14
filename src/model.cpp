#include "model.h"
#include "eigen.h"

const char *modelNames[numModels]={
	"JC69",
	"K80",
	"K81",
	"K81ne",
	"F81",
	"F84",
	"HKY",
	"T92",
	"TN93",
	"TN93eq",
	"TIM",
	"TIMeq",
	"TVM",
	"TVMeq",
	"SYM",
	"GTR",
	"JTT",
	"WAG",
	"PAM",
	"BLOSUM",
	"MTREV",
	"CPREV",
	"GENERAL"
};

const char *modelTitles[numModels]={
	"JC69:   Jukes TH and Cantor CR (1969). Evolution of Protein Molecules. New York: Academic Press. p. 21Ð132",
	"K80:    Kimura M (1980). \"A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences\". Journal of Molecular Evolution 16: 111Ð120",
	"K81:    ",
	"K81ne   ",
	"F81:    Felsenstein J (1981). \"Evolutionary trees from DNA sequences: a maximum likelihood approach\". Journal of Molecular Evolution 17: 368Ð376",
	"F84:    Felsenstein (1984)",
	"HKY:    Hasegawa, Kishino & Yano (1985). \"Dating of human-ape splitting by a molecular clock of mitochondrial DNA\". Journal of Molecular Evolution 22: 160Ð174",
	"T92:    Tamura K (1992). \"Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C content biases\". Molecular Biology and Evolution 9: 678Ð687",
	"TN93:   Tamura K, Nei M (1993). \"Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees\". Molecular Biology and Evolution 10: 512Ð526",
	"TN93eq: ",
	"TIM:    ",
	"TIMeq:  ",
	"TVM:    ",
	"TVMeq:  ",
	"SYM:    Zharkikh (1994)",
	"GTR:    General time reversible (nucleotides)",
	"JTT:    Jones, Taylor & Thornton (1992) CABIOS  8:275-282\n             DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
	"WAG:    Whelan & Goldman (2001) Mol Biol Evol 18:691Ð699",
	"PAM:    Dayhoff, Schwartz & Orcutt (1978)\n             DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
	"BLOSUM62: Henikoff & Henikoff (1992) PNAS USA 89:10915-10919",
	"MTREV24: Adachi & Hasegawa (1996) J Mol Evol 42:459-468",
	"CPREV45: Adachi et al. (2000) J Mol Evol 50:348-358",
	"GENERAL: General time reversible (amino acids)"
};

// JTT model for amino acid evolution -- DCMut
// D.T. Jones, W.R. Taylor, and J.M. Thornton 
// The rapid generation of mutation data matrices from protein sequences 
// CABIOS  vol. 8 no. 3 1992 pp. 275-282 
double RateMatrix::jttRelativeRates[NUM_AA_REL_RATES] = {
		0.531678, 0.557967, 0.827445, 0.574478, 0.556725, 1.066681, 1.740159, 0.219970, 0.361684, 0.310007, 0.369437, 0.469395, 0.138293, 1.959599, 3.887095, 4.582565, 0.084329, 0.139492, 2.924161,
		0.451095, 0.154899, 1.019843, 3.021995, 0.318483, 1.359652, 3.210671, 0.239195, 0.372261, 6.529255, 0.431045, 0.065314, 0.710489, 1.001551, 0.650282, 1.257961, 0.235601, 0.171995,
		5.549530, 0.313311, 0.768834, 0.578115, 0.773313, 4.025778, 0.491003, 0.137289, 2.529517, 0.330720, 0.073481, 0.121804, 5.057964, 2.351311, 0.027700, 0.700693, 0.164525,
		0.105625, 0.521646, 7.766557, 1.272434, 1.032342, 0.115968, 0.061486, 0.282466, 0.190001, 0.032522, 0.127164, 0.589268, 0.425159, 0.057466, 0.453952, 0.315261,
		0.091304, 0.053907, 0.546389, 0.724998, 0.150559, 0.164593, 0.049009, 0.409202, 0.678335, 0.123653, 2.155331, 0.469823, 1.104181, 2.114852, 0.621323,
		3.417706, 0.231294, 5.684080, 0.078270, 0.709004, 2.966732, 0.456901, 0.045683, 1.608126, 0.548807, 0.523825, 0.172206, 0.254745, 0.179771,
		1.115632, 0.243768, 0.111773, 0.097485, 1.731684, 0.175084, 0.043829, 0.191994, 0.312449, 0.331584, 0.114381, 0.063452, 0.465271,
		0.201696, 0.053769, 0.069492, 0.269840, 0.130379, 0.050212, 0.208081, 1.874296, 0.316862, 0.544180, 0.052500, 0.470140,
		0.181788, 0.540571, 0.525096, 0.329660, 0.453428, 1.141961, 0.743458, 0.477355, 0.128193, 5.848400, 0.121827,
		2.335139, 0.202562, 4.831666, 0.777090, 0.098580, 0.405119, 2.553806, 0.134510, 0.303445, 9.533943,
		0.146481, 3.856906, 2.500294, 1.060504, 0.592511, 0.272514, 0.530324, 0.241094, 1.761439,
		0.624581, 0.024521, 0.216345, 0.474478, 0.965641, 0.089134, 0.087904, 0.124066,
		0.436181, 0.164215, 0.285564, 2.114728, 0.201334, 0.189870, 3.038533,
		0.148483, 0.943971, 0.138904, 0.537922, 5.484236, 0.593478,
		2.788406, 1.176961, 0.069965, 0.113850, 0.211561,
		4.777647, 0.310927, 0.628608, 0.408532,
		0.080556, 0.201094, 1.14398,
		0.747889, 0.239697,
		0.165473
	};
double RateMatrix::jttFrequencies[NUM_AA] = {
		0.076862, 0.051057, 0.042546, 0.051269, 0.020279, 0.041061, 0.061820, 0.074714, 0.022983, 0.052569, 0.091111, 0.059498, 0.023414, 0.040530, 0.050532, 0.068225, 0.058518, 0.014336, 0.032303, 0.066374
	};
	
// WAG model of amino acid evolution 
// Whelan, S. and Goldman, N. (2001) A general empirical model of protein 
// evolution derived from multiple protein families using a maximum-likelihood 
// approach. Mol. Biol. Evol. 18, 691Ð699. 
double RateMatrix::wagRelativeRates[NUM_AA_REL_RATES] = {
		0.610810, 0.569079, 0.821500, 1.141050, 1.011980, 1.756410, 1.572160, 0.354813, 0.219023, 0.443935, 1.005440, 0.989475, 0.233492, 1.594890, 3.733380, 2.349220, 0.125227, 0.268987, 2.221870, 
		0.711690, 0.165074, 0.585809, 3.360330, 0.488649, 0.650469, 2.362040, 0.206722, 0.551450, 5.925170, 0.758446, 0.116821, 0.753467, 1.357640, 0.613776, 1.294610, 0.423612, 0.280336, 
		6.013660, 0.296524, 1.716740, 1.056790, 1.253910, 4.378930, 0.615636, 0.147156, 3.334390, 0.224747, 0.110793, 0.217538, 4.394450, 2.257930, 0.078463, 1.208560, 0.221176, 
		0.033379, 0.691268, 6.833400, 0.961142, 1.032910, 0.043523, 0.093930, 0.533362, 0.116813, 0.052004, 0.472601, 1.192810, 0.417372, 0.146348, 0.363243, 0.169417, 
		0.109261, 0.023920, 0.341086, 0.275403, 0.189890, 0.428414, 0.083649, 0.437393, 0.441300, 0.122303, 1.560590, 0.570186, 0.795736, 0.604634, 1.114570, 
		6.048790, 0.366510, 4.749460, 0.131046, 0.964886, 4.308310, 1.705070, 0.110744, 1.036370, 1.141210, 0.954144, 0.243615, 0.252457, 0.333890, 
		0.630832, 0.635025, 0.141320, 0.172579, 2.867580, 0.353912, 0.092310, 0.755791, 0.782467, 0.914814, 0.172682, 0.217549, 0.655045, 
		0.276379, 0.034151, 0.068651, 0.415992, 0.194220, 0.055288, 0.273149, 1.486700, 0.251477, 0.374321, 0.114187, 0.209108, 
		0.152215, 0.555096, 0.992083, 0.450867, 0.756080, 0.771387, 0.822459, 0.525511, 0.289998, 4.290350, 0.131869, 
		3.517820, 0.360574, 4.714220, 1.177640, 0.111502, 0.353443, 1.615050, 0.234326, 0.468951, 8.659740, 
		0.287583, 5.375250, 2.348200, 0.462018, 0.382421, 0.364222, 0.740259, 0.443205, 1.997370, 
		1.032220, 0.098843, 0.619503, 1.073780, 1.537920, 0.152232, 0.147411, 0.342012, 
		1.320870, 0.194864, 0.556353, 1.681970, 0.570369, 0.473810, 2.282020, 
		0.179896, 0.606814, 0.191467, 1.699780, 7.154480, 0.725096, 
		1.786490, 0.885349, 0.156619, 0.239607, 0.351250, 
		4.847130, 0.578784, 0.872519, 0.258861, 
		0.126678, 0.325490, 1.547670, 
		2.763540, 0.409817, 
		0.347826
	};
double RateMatrix::wagFrequencies[NUM_AA] = {
		0.0866, 0.0440, 0.0391, 0.0570, 0.0193, 0.0367, 0.0581, 0.0833, 0.0244, 0.0485, 0.0862, 0.0620, 0.0195, 0.0384, 0.0458, 0.0695, 0.0610, 0.0144, 0.0353, 0.0709
	};

// Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1978)
// A model of evolutionary change in proteins.
// Dayhoff, M.O. (ed.) Atlas of Protein Sequence Structur., Vol5, Suppl. 3,
// National Biomedical Research Foundation, Washington DC, pp. 345-352.
double RateMatrix::dayhoffRelativeRates[NUM_AA_REL_RATES] = {
		0.267828, 0.984474, 1.199805, 0.360016, 0.887753, 1.961167, 2.386111, 0.228116, 0.653416, 0.406431, 0.258635, 0.717840, 0.183641, 2.485920, 4.051870, 3.680365, 0.000000, 0.244139, 2.059564,
		0.327059, 0.000000, 0.232374, 2.439939, 0.000000, 0.087791, 2.383148, 0.632629, 0.154924, 4.610124, 0.896321, 0.136906, 1.028313, 1.531590, 0.265745, 2.001375, 0.078012, 0.240368,
		8.931515, 0.000000, 1.028509, 1.493409, 1.385352, 5.290024, 0.768024, 0.341113, 3.148371, 0.000000, 0.138503, 0.419244, 4.885892, 2.271697, 0.224968, 0.946940, 0.158067,
		0.00000, 1.348551, 11.388659, 1.240981, 0.868241, 0.239248, 0.000000, 0.716913, 0.000000, 0.000000, 0.133940, 0.956097, 0.660930, 0.000000, 0.000000, 0.178316,
		0.000000, 0.000000, 0.107278, 0.282729, 0.438074, 0.000000, 0.000000, 0.000000, 0.000000, 0.187550, 1.598356, 0.162366, 0.000000, 0.953164, 0.484678,
		7.086022, 0.281581, 6.011613, 0.180393, 0.730772, 1.519078, 1.127499, 0.000000, 1.526188, 0.561828, 0.525651, 0.000000, 0.000000, 0.346983,
		0.811907, 0.439469, 0.609526, 0.112880, 0.830078, 0.304803, 0.000000, 0.507003, 0.793999, 0.340156, 0.000000, 0.214717, 0.36725,
		0.106802, 0.000000, 0.071514, 0.267683, 0.170372, 0.153478, 0.347153, 2.322243, 0.306662, 0.000000, 0.000000, 0.538165,
		0.076981, 0.443504, 0.270475, 0.000000, 0.475927, 0.933709, 0.353643, 0.226333, 0.270564, 1.265400, 0.438715,
		2.556685, 0.460857, 3.332732, 1.951951, 0.119152, 0.247955, 1.900739, 0.000000, 0.374834, 8.810038,
		0.180629, 5.230115, 1.565160, 0.316258, 0.171432, 0.331090, 0.461776, 0.286572, 1.745156,
		2.411739, 0.000000, 0.335419, 0.954557, 1.350599, 0.000000, 0.132142, 0.10385,
		0.921860, 0.170205, 0.619951, 1.031534, 0.000000, 0.000000, 2.565955,
		0.110506, 0.459901, 0.136655, 0.762354, 6.952629, 0.123606,
		2.427202, 0.782857, 0.000000, 0.000000, 0.485026,
		5.436674, 0.740819, 0.336289, 0.303836,
		0.000000, 0.417839, 1.561997,
		0.608070, 0.000000,
		0.279379
	};
double RateMatrix::dayhoffFrequencies[NUM_AA] = {
		0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255, 0.049530, 0.088612, 0.033619, 0.036886, 0.085357, 0.080481, 0.014753, 0.039772, 0.050680, 0.069577, 0.058542, 0.010494, 0.029916, 0.064718
	};
	
// BLOSUM62 model of amino acid evolution 
// Henikoff, S., and J. G. Henikoff. 1992. PNAS USA 89:10915-10919.
double RateMatrix::blosumRelativeRates[NUM_AA_REL_RATES] = {
		7.3579038969751e-01, 4.8539105546575e-01, 5.4316182089867e-01, 1.4599953104700e+00, 1.1997057046020e+00, 1.1709490427999e+00, 1.9558835749595e+00, 7.1624144499779e-01, 6.0589900368677e-01, 8.0001653051838e-01, 1.2952012667833e+00, 1.2537582666635e+00, 4.9296467974759e-01, 1.1732759009239e+00, 4.3250926870566e+00, 1.7291780194850e+00, 4.6583936772479e-01, 7.1820669758623e-01, 2.1877745220045e+00, 
		1.2974467051337e+00, 5.0096440855513e-01, 2.2782657420895e-01, 3.0208336100636e+00, 1.3605741904203e+00, 4.1876330851753e-01, 1.4561411663360e+00, 2.3203644514174e-01, 6.2271166969249e-01, 5.4111151414889e+00, 9.8369298745695e-01, 3.7164469320875e-01, 4.4813366171831e-01, 1.1227831042096e+00, 9.1466595456337e-01, 4.2638231012175e-01, 7.2051744121611e-01, 4.3838834377202e-01, 
		3.1801000482161e+00, 3.9735894989702e-01, 1.8392161469920e+00, 1.2404885086396e+00, 1.3558723444845e+00, 2.4145014342081e+00, 2.8301732627800e-01, 2.1188815961519e-01, 1.5931370434574e+00, 6.4844127878707e-01, 3.5486124922252e-01, 4.9488704370192e-01, 2.9041016564560e+00, 1.8981736345332e+00, 1.9148204624678e-01, 5.3822251903674e-01, 3.1285879799342e-01, 
		2.4083661480204e-01, 1.1909457033960e+00, 3.7616252083685e+00, 7.9847324896839e-01, 7.7814266402188e-01, 4.1855573246161e-01, 2.1813157759360e-01, 1.0324479249521e+00, 2.2262189795786e-01, 2.8173069420651e-01, 7.3062827299842e-01, 1.5827541420653e+00, 9.3418750943056e-01, 1.4534504627853e-01, 2.6142220896504e-01, 2.5812928941763e-01, 
		3.2980150463028e-01, 1.4074889181440e-01, 4.1820319228376e-01, 3.5405810983129e-01, 7.7489402279418e-01, 8.3184264014158e-01, 2.8507880090648e-01, 7.6768882347954e-01, 4.4133747118660e-01, 3.5600849876863e-01, 1.1971884150942e+00, 1.1198313585160e+00, 5.2766441887169e-01, 4.7023773369610e-01, 1.1163524786062e+00, 
		5.5289191779282e+00, 6.0984630538281e-01, 2.4353411311401e+00, 2.3620245120365e-01, 5.8073709318144e-01, 3.9452776745146e+00, 2.4948960771127e+00, 1.4435695975031e-01, 8.5857057567418e-01, 1.9348709245965e+00, 1.2774802945956e+00, 7.5865380864172e-01, 9.5898974285014e-01, 5.3078579012486e-01, 
		4.2357999217628e-01, 1.6268910569817e+00, 1.8684804693170e-01, 3.7262517508685e-01, 2.8024271516787e+00, 5.5541539747043e-01, 2.9140908416530e-01, 9.2656393484598e-01, 1.7698932389373e+00, 1.0710972360073e+00, 4.0763564893830e-01, 5.9671930034577e-01, 5.2425384633796e-01, 
		5.3985912495418e-01, 1.8929629237636e-01, 2.1772115923623e-01, 7.5204244030271e-01, 4.5943617357855e-01, 3.6816646445253e-01, 5.0408659952683e-01, 1.5093262532236e+00, 6.4143601140497e-01, 5.0835892463812e-01, 3.0805573703500e-01, 2.5334079019018e-01, 
		2.5271844788492e-01, 3.4807220979697e-01, 1.0225070358890e+00, 9.8431152535870e-01, 7.1453370392764e-01, 5.2700733915060e-01, 1.1170297629105e+00, 5.8540709022472e-01, 3.0124860078016e-01, 4.2189539693890e+00, 2.0155597175031e-01, 
		3.8909637733035e+00, 4.0619358664202e-01, 3.3647977631042e+00, 1.5173593259539e+00, 3.8835540920564e-01, 3.5754441245967e-01, 1.1790911972601e+00, 3.4198578754023e-01, 6.7461709322842e-01, 8.3118394054582e+00, 
		4.4557027426059e-01, 6.0305593795716e+00, 2.0648397032375e+00, 3.7455568747097e-01, 3.5296918452729e-01, 9.1525985769421e-01, 6.9147463459998e-01, 8.1124585632307e-01, 2.2314056889131e+00, 
		1.0730611843319e+00, 2.6692475051102e-01, 1.0473834507215e+00, 1.7521659178195e+00, 1.3038752007987e+00, 3.3224304063396e-01, 7.1799348690032e-01, 4.9813847530407e-01, 
		1.7738551688305e+00, 4.5412362510273e-01, 9.1872341574605e-01, 1.4885480537218e+00, 8.8810109815193e-01, 9.5168216224591e-01, 2.5758507553153e+00, 
		2.3359790962888e-01, 5.4002764482413e-01, 4.8820611879305e-01, 2.0743248934965e+00, 6.7472604308008e+00, 8.3811961017754e-01, 
		1.1691295777157e+00, 1.0054516831488e+00, 2.5221483002727e-01, 3.6940531935451e-01, 4.9690841067567e-01, 
		5.1515562922704e+00, 3.8792562209837e-01, 7.9675152076106e-01, 5.6192545744165e-01, 
		5.1312812689059e-01, 8.0101024319939e-01, 2.2530740511763e+00, 
		4.0544190065580e+00, 2.6650873142646e-01, 
		1.0000000000000e+00
	};
double RateMatrix::blosumFrequencies[NUM_AA] = {
		0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.032, 0.073
	};

// mtRev model - complete sequence data of mtDNA from 24 vertebrate species 
// Adachi, J., and Hasegawa, M. 1996. J. Mol. Evol. 42:459-468. 
double RateMatrix::mtrevRelativeRates[NUM_AA_REL_RATES] = {
		1.2199217606346e+01, 1.4182139942122e+01, 9.2985091873208e+00, 3.1542792981957e+01, 1.0025852846688e+00, 5.1418866803338e+00, 6.3531246495131e+01, 7.3137132861715e+00, 5.0782382656186e+01, 1.3399741808481e+01, 4.4021672780560e+00, 7.4673480520104e+01, 3.3513021631978e+00, 2.8582502221773e+01, 2.0413623195312e+02, 2.5301305153906e+02, 1.0000000000000e+00, 3.4084158197615e+00, 1.0266468401249e+02, 
		6.9661274444534e+00, 1.0000000000000e+00, 5.4384584796568e+01, 1.1631134513343e+02, 1.0000000000000e+00, 1.2122831341194e+01, 8.6961067087353e+01, 1.0000000000000e+00, 8.1976829394538e+00, 7.4423215395318e+01, 1.0000000000000e+00, 2.4659158338099e+00, 1.2439947713615e+01, 3.1791814866372e+00, 1.0935327216119e+00, 1.1550775790126e+01, 1.0000000000000e+00, 4.0211417480338e+00, 
		4.1809325468160e+02, 3.1020979842967e+01, 9.1349622725361e+01, 3.3185663516310e+01, 2.8052324651124e+01, 2.6112087577885e+02, 1.4261453863336e+01, 7.9775653461977e+00, 3.2036829276162e+02, 3.4424354918739e+01, 7.9996445145608e+00, 3.8586541461044e+01, 2.6020426225852e+02, 1.2550758780474e+02, 5.6207759736659e+00, 1.0071406219571e+02, 1.0000000000000e+00, 
		1.0000000000000e+00, 2.9097352675564e+01, 3.0713149855302e+02, 2.9877072751897e+01, 5.9995408885817e+01, 2.2827096245105e+00, 1.0000000000000e+00, 1.2183938185384e+00, 1.0000000000000e+00, 2.6221929413096e+00, 7.0708004204733e+00, 3.6327934317139e+01, 1.4743408713748e+01, 1.0453246057102e+01, 1.1165627147496e+01, 1.0000000000000e+00, 
		3.9599394038972e+01, 1.0000000000000e+00, 1.6163581056674e+01, 7.4467985406234e+01, 3.3018175376623e+01, 1.3500725995091e+01, 1.0000000000000e+00, 3.2504095376923e+00, 3.7264767083096e+01, 1.6454136037822e+01, 1.4581783243113e+02, 9.4720031458442e+01, 1.7684087896962e+01, 1.3409157685926e+02, 1.0000000000000e+00, 
		1.6503249008836e+02, 3.5530760735494e+00, 3.0652523140859e+02, 4.3905393139325e+00, 2.0895470525345e+01, 2.4504076430724e+02, 2.4931300477797e+01, 1.0059428264289e+01, 7.2256314165467e+01, 2.8480937892158e+01, 4.9962974409828e+01, 1.0000000000000e+00, 2.0430790980529e+01, 9.9986289000676e+00, 
		1.4884496769963e+01, 2.5853576435567e+01, 1.7418201388328e+00, 1.0000000000000e+00, 1.6519126809071e+02, 1.0000000000000e+00, 1.4067850525292e+00, 6.7547121641947e+00, 2.8794794140840e+01, 7.8001372062558e+00, 1.0000000000000e+00, 6.9067239183061e+00, 1.1127702362585e+01, 
		1.0000000000000e+00, 3.1466649021550e+00, 1.2699794194865e+00, 1.1962111069278e+01, 1.0000000000000e+00, 1.0000000000000e+00, 1.0000000000000e+00, 6.6277950574411e+01, 5.8800079133028e+00, 5.7494182626674e+00, 1.6887657206208e+00, 1.3320553471351e+00, 
		6.4536986087271e+00, 6.0472584534958e+00, 6.7197196398961e+01, 6.2977633277779e+00, 2.5347805183364e+01, 3.2089868698728e+01, 4.0766987134407e+01, 2.3570850628539e+01, 3.7286635325194e+00, 3.5270764890474e+02, 1.0000000000000e+00, 
		1.7320653206333e+02, 1.0298655619743e+01, 2.7262244199514e+02, 4.4561065036310e+01, 1.0856482766156e+01, 2.5107659603898e+01, 1.9391167162525e+02, 1.0000000000000e+00, 1.3161329199391e+01, 6.4365086389428e+02, 
		7.8314019154706e+00, 2.8290920517725e+02, 1.1371735519833e+02, 2.1105885757279e+01, 3.8741359395934e+01, 6.6524559321657e+01, 1.7071378554833e+01, 2.3234516108847e+01, 4.8247261078055e+01, 
		4.8092094826036e+01, 3.3887559483420e+00, 2.6368577564199e+01, 5.5679895711418e+01, 7.1750284708933e+01, 1.2631893872825e+01, 2.6932728996777e+01, 1.0000000000000e+00, 
		4.7798798034572e+01, 9.9165053447429e+00, 5.8505442466161e+01, 2.7798190504760e+02, 1.1427000119701e+01, 2.1029990530586e+01, 2.0397078683768e+02, 
		9.1089574817139e+00, 3.3835737720574e+01, 1.7815549567056e+01, 4.1272404968214e+00, 2.4504156395152e+02, 3.3435675442163e+00, 
		8.9421193040709e+01, 6.7485067008375e+01, 2.2161693733113e+00, 8.5338209390745e+00, 4.3342126659660e+00, 
		3.1432036618746e+02, 2.0305343047059e+01, 3.4167877957799e+01, 1.0000000000000e+00, 
		5.2559565123081e+00, 2.0382362288681e+01, 1.0765527137500e+02, 
		1.3814733274637e+01, 2.8259139240676e+00, 
		1.0000000000000e+00
	};
double RateMatrix::mtrevFrequencies[NUM_AA] = {
		0.072, 0.019, 0.039, 0.019, 0.006, 0.025, 0.024, 0.056, 0.028, 0.088, 0.168, 0.023, 0.054, 0.061, 0.054, 0.072, 0.086, 0.029, 0.033, 0.043
	};

// CPREV 45 model of amino acid evolution 
// Adachi, J., P.J. Waddell, W. Martin, and M. Hasegawa. 2000. JME 50:348-358 
double RateMatrix::cprevRelativeRates[NUM_AA_REL_RATES] = {
		-105, -227, -175, -669, -157, -499, -665, -66, -145, -197, -236, -185, -68, -490, -2440, -1340, -14, -56, -968, 
		-357, -43, -823, -1745, -152, -243, -715, -136, -203, -4482, -125, -53, -87, -385, -314, -230, -323, -92, 
		-4435, -538, -768, -1055, -653, -1405, -168, -113, -2430, -61, -97, -173, -2085, -1393, -40, -754, -83, 
		-10, -400, -3691, -431, -331, -10, -10, -412, -47, -22, -170, -590, -266, -18, -281, -75, 
		-10, -10, -303, -441, -280, -396, -48, -159, -726, -285, -2331, -576, -435, -1466, -592, 
		-3122, -133, -1269, -92, -286, -3313, -202, -10, -323, -396, -241, -53, -391, -54, 
		-379, -162, -148, -82, -2629, -113, -145, -185, -568, -369, -63, -142, -200, 
		-19, -40, -20, -263, -21, -25, -28, -691, -92, -82, -10, -91, 
		-29, -66, -305, -10, -127, -152, -303, -32, -69, -1971, -25, 
		-1745, -345, -1772, -454, -117, -216, -1040, -42, -89, -4797, 
		-218, -1351, -1268, -219, -516, -156, -159, -189, -865, 
		-193, -72, -302, -868, -918, -10, -247, -249, 
		-327, -100, -93, -645, -86, -215, -475, 
		-43, -487, -148, -468, -2370, -317, 
		-1202, -260, -49, -97, -122, 
		-2151, -73, -522, -167, 
		-29, -71, -760, 
		-346, -10, 
		-119
	};
double RateMatrix::cprevFrequencies[NUM_AA] = {
		0.076, 0.062, 0.041, 0.037, 0.009, 0.038, 0.049, 0.084, 0.025, 0.081, 0.101, 0.050, 0.022, 0.051, 0.043, 0.062, 0.054, 0.018, 0.031, 0.066
	};

int model, numStates, isNucModel;

std::string stateCharacters;

void RateMatrix::InitializeSubstitutionVectors()
{
	int size;
		
	if (isNucModel) size = 4;
	else size = 20;

	pi.assign(size, 1.0/size);
	Qij.assign(size*size, 1);
	Sij.assign(size*size-size, 1.0);
	Cijk.assign(size*size*size, 1.0);
	Pij.assign(MAX_RATE_CATS, vector<double> (size*size));
	catRate.assign(MAX_RATE_CATS, 0.0);
	catRate.at(0) = 1.0;		/// If there are no categories.
	freqRate.assign(MAX_RATE_CATS, 0.0);
}

void RateMatrix::setModelConditions(
									string& rate_matrix, 
									string& rates_in, 
									string& freqs_in
								   )
{
	size_t j = 0;
	list<string> arg_split;
	vector<double>::iterator R_it;

	substitution_model = FindModel(rate_matrix);
	InitializeSubstitutionVectors();

	//////////
	/// Each entry in the input of frequencies is separated by the ',' token.
	//////////
	arg_split = split(freqs_in, ",");
	if (arg_split.size() == numStates) {
		vector<double>::iterator pi_t = pi.begin();
		for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it, ++pi_t)
			(*pi_t) = atof((*it).c_str());
	} else if ( 
			   substitution_model == SYM 
			   || substitution_model == JC69 
			   || substitution_model == K80 
			   || substitution_model == K81 
			   || substitution_model == TN93eq 
			   || substitution_model == TIMeq
			   || substitution_model == TVMeq
			   || numStates > 4 
			  ) {
		// Do not need to set frequencies from command-line. Already have preset  values.
	} else cerr << "setModelConditions: bad freqs in." << endl;

	arg_split = split(rates_in, ",");
	switch (substitution_model) {
		case SYM:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case GTR:
			if(arg_split.size() != NUM_NUC_REL_RATES) {
				cout << "GTR and SYM rate matrices require " << NUM_NUC_REL_RATES << " values in option --rel_rates (-r), only " << arg_split.size() << " supplied:" << endl;
				for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it, j++)
					cout << j << ": " << (*it) << endl;
				cout << endl;
				cout << " <list> = a,b,c,d,e,f--------------------------\n";
				cout << "     model GTR/SYM\n";
				cout << "           A C G T   \n";
				cout << "         | * a b c | \n";
				cout << "     Q = | a * d e | \n";
			    cout << "         | b d * f | \n";
			    cout << "         | c e f * | \n";
				cout << "\n";
				cout << " ----------------------------------------------\n";
				cout << "\n";
				exit(EXIT_FAILURE);
			}
			R_it = Sij.begin();
			for(list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it, ++R_it) 
				(*R_it) = atof((*it).c_str());
			if(Sij.back() != 1.0) {
				for(R_it = Sij.begin(); R_it != Sij.end(); ++R_it) 
					(*R_it) /= Sij.back();
			}
			break;
		case K80:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case F84:
		case HKY:
			if (arg_split.size() != 1) {
				cout << "K80 and HKY rate matrices require 1 values (kappa) in option --rel_rates (-r): " << rates_in << endl << endl;
				cout << " <list> = a -----------------------------------\n";
				cout << "     model K80/HKY/F84\n";
				cout << "           A C G T  \n";
				cout << "         | *   a   |\n";
				cout << "     Q = |   *   a |\n";
				cout << "         | a   *   |\n";
				cout << "         |   a   * |\n";
				cout << "\n";
				cout << " ----------------------------------------------\n";
				cout << "\n";
				exit(EXIT_FAILURE);
			}
			for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it) 
				kappa.push_back(atof((*it).c_str()));
		break;
		case TN93eq:
		case K81:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case TN93:
		case K81ne:
			if (arg_split.size() != 2) {
				cout << "K81 and TN93 rate matrices require 2 values in option --rel_rates (-r): " << rates_in << endl << endl;
				cout << " <list> = a,b ---------------------------------\n";
				cout << "     model K81/K81ne    model TN93/TN93eq\n";
				cout << "           A C G T            A C G T\n";
				cout << "         | *   a b |        | *   a   |\n";
				cout << "     Q = |   * b a |    Q = |   *   b |\n";
			    cout << "         | a b *   |        | a   *   |\n";
			    cout << "         | b a   * |        |   b   * |\n";
				cout << "\n";
				cout << " ----------------------------------------------\n";
				cout << "\n";
				exit(EXIT_FAILURE);
			}
			for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it) 
				kappa.push_back(atof((*it).c_str()));
			break;
		case TIMeq:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case TIM:
			if (arg_split.size() != 3) {
				cout << "TIM rate matrix require 3 values in option --rel_rates (-r): " << rates_in << endl << endl;
				cout << " <list> = a,b,c -------------------------------\n";
				cout << "     model TIM/TIMeq\n";
				cout << "           A C G T   \n";
				cout << "         | *   a b | \n";
				cout << "     Q = |   * b c | \n";
			    cout << "         | a b *   | \n";
			    cout << "         | b c   * | \n";
				cout << "\n";
				cout << " ----------------------------------------------\n";
				cout << "\n";
				exit(EXIT_FAILURE);
			}
			for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it) 
				kappa.push_back(atof((*it).c_str()));
			break;
		case TVMeq:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case TVM:
			if (arg_split.size() != 4) {
				cout << "TVM rate matrix require 4 values in option --rel_rates (-r): " << rates_in << endl << endl;
				cout << " <list> = a,b,c,d -----------------------------\n";
				cout << "     model TVM/TVMeq\n";
				cout << "           A C G T   \n";
				cout << "         | *   a b | \n";
				cout << "     Q = |   * c b | \n";
			    cout << "         | a c * d | \n";
			    cout << "         | b b d * | \n";
				cout << "\n";
				cout << " ----------------------------------------------\n";
				cout << "\n";
				exit(EXIT_FAILURE);
			}
			for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it) 
				kappa.push_back(atof((*it).c_str()));
			break;
		case JC69:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case F81:
			if ( arg_split.empty() ) {
				cerr << "Input model does not require relative rates: " << modelNames[substitution_model] << endl;
				exit(EXIT_FAILURE);
			}
			break;
		case JTT:
		case PAM:
		case WAG:
		case BLOSUM:
		case MTREV:
		case CPREV:
			if ( arg_split.empty() ) {
				cerr << "Protein models do not require relative rates: " << modelNames[substitution_model] << endl;
				exit(EXIT_FAILURE);
			}
			break;
		default:
			cerr << "Unknown input model: " << rate_matrix << endl;
			exit(EXIT_FAILURE);
		break;
	}
}

int RateMatrix::FindModel(
						  string& theModel
						 ) 
{
	int model = -1;

	trim(theModel);

	for(size_t i = 0; i < 3; i++) 
		theModel.at(i) = toupper(theModel.at(i));
	if (theModel.compare("TRN") == 0) theModel.assign("TN93");
	if (theModel.compare("TRNeq") == 0) theModel.assign("TN93eq");
	if (theModel.compare("K2P") == 0) theModel.assign("K80");
	if (theModel.compare("K3P") == 0) theModel.assign("K81");
	if (theModel.compare("K3Pne") == 0) theModel.assign("K81ne");
	for(int i = JC69; i<numModels; i++) {
		if(theModel.compare(modelNames[i]) == 0) {
			model = i;
			if (model <= GTR) {
				if(isNucModel == -1) {
					isNucModel = 1;
					numStates = NUM_NUC;
				} else {
					if(!isNucModel) {
						cerr << "Substitution model " << modelNames[model] << " is a nucleotide model, but simulation run is for amino acids." << endl;
						exit(EXIT_FAILURE);							
					}
				}
			} else {
				if(isNucModel == -1) {
					isNucModel = 0;
					numStates = NUM_AA;
				} else {
					if(isNucModel) {
						cerr << "Substitution model " << modelNames[model] << " is an amino acid model, but simulation run is for nucleotides." << endl;
						exit(EXIT_FAILURE);
					}
				}
			}
		} else if (theModel.compare(0,3,"REV",0,3)==0) {
			if(isNucModel == -1) {
				isNucModel = 1;
				numStates = NUM_NUC;
			} else {
				if(!isNucModel) {
					
				}
			}
			model = GTR;
		}
	}
	if(model == -1) {
		cerr << "Unknown model: " << theModel << endl << endl;
		exit(EXIT_FAILURE);
	}

	numStates_squared = numStates * numStates;
	numStates_cubed = numStates_squared * numStates;

	return model;
}

void RateMatrix::SetModel(
						  inClade *branch
						 )
{	
	int i;
	bool freqSet = true;
	double a, b, c, d, e, f, a_, b_, c_, d_, e_, f_, kappa1, kappa2, alpha, beta, gamma, delta;
	double sum_freq = 0;	
	
	cerr << "Point-> RateMatrix::SetModel(inClade *branch) IN" << endl;
	
	if (isNucModel) stateCharacters = nucleotides;
	else stateCharacters = aminoAcids;

	//////////
	/// freqSet does not matter for nucleotides. It should always be set for nuc models that do
	/// not require equal frequencies. For amino acid models, though, need to know if the freqs
	/// are set, otherwise they can be overwritten by model-specific pi values.
	//////////
	if (branch != NULL) {
		if ( !(branch->values2Export2Freq.empty()) ) {
			for (int j = 0; j < numStates; j++)
				pi.at(j) = branch->values2Export2Freq.at(j);
		} else freqSet = false;
	} else freqSet = false;

	for (vector<double>::iterator it = pi.begin(); it != pi.end(); ++it)
		sum_freq += (*it);
	for (vector<double>::iterator it = pi.begin(); it != pi.end(); ++it)
		(*it) /= sum_freq;

	//////////
	/// Putting model parameters explicitly, easier to understand.
	//////////
	switch (substitution_model)
	{
		case NUC_JC69:
			////////// OK
			///     | * 1 1 1 |
			/// Q = | 1 * 1 1 |
			///     | 1 1 * 1 |
			///     | 1 1 1 * |
			//////////
			a=b=c=d=e=f=a_=b_=c_=d_=e_=f_=1;
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		break;
		case NUC_K80:
			////////// OK
			///     | * 1 K 1 |
			/// Q = | 1 * 1 K |
			///     | K 1 * 1 |
			///     | 1 K 1 * |
			//////////
			kappa1 = kappa.front();
			b = b_ = e = e_ = kappa1;
			a = a_ = c = c_ = d = d_ = f = f_ = 1.0;
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		break;
		case NUC_K81: 
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case NUC_K81ne:
			////////// OK
			///     | *  1  K A |
			/// Q = | 1  *  A K |
			///     | K A  *  1 |
			///     | A K  1  * |
			//////////
			kappa1 = kappa.at(0);
			alpha = kappa.at(1);
			a = 1; b = kappa1; c = alpha;	// *fA
			a_ = 1; d = alpha; e = kappa1;	// *fC
			b_ = kappa1; d_ = alpha; f = 1;	// *fG
			c_ = alpha; e_ = kappa1; f_ = 1;	// *fT
		break;
		case NUC_F81:
			////////// OK
			///     | * fA  fA fA |
			/// Q = | fC *  fC fC |
			///     | fG fG *  fG |
			///     | fT fT fT *  |
			//////////
			a=b=c=d=e=f=a_=b_=c_=d_=e_=f_=1;
		break;
		case NUC_HKY:
			////////// OK
			///     | *   fA  KfA  fA  |
			/// Q = | fC  *   fC   KfC |
			///     | KfG fG  *    fG  |
			///     | fT  KfT fT   *   |
			//////////
			kappa1 = kappa.at(0);
			a = c = 1; b = kappa1;		// *fA
			a_ = d = 1; e = kappa1;		// *fC
			b_ = kappa1; d_ = f = 1;		// *fG
			c_ = f_ = 1; e_ = kappa1;	// *fT
		break;
		case NUC_F84:
			////////// OK
			/// INPUT: kappa
			///     | *               fA              (1+K/(fA+fG))fA fA              |
			/// Q = | fC              *               fC              (1+K/(fC+fT))fC |
			///     | (1+K/(fA+fG))fG fG              *               fG              |
			///     | fT              (1+K/(fC+fT))fT fT              *               |
			//////////
			kappa1 = kappa.at(0);
			a = c = 1; b = (1.0+kappa1/(pi.at(A)+pi.at(G)));		// *fA
			a_ = d = 1; e = (1.0+kappa1/(pi.at(C)+pi.at(T)));		// *fC
			b_ = (1.0+kappa1/(pi.at(A)+pi.at(G))); d_ = f = 1;		// *fG
			c_ = f_ = 1; e_ = (1.0+kappa1/(pi.at(C)+pi.at(T)));		// *fT
		break;
		case NUC_T92:
			//////////
			/// INPUT: kappa
			///     | *             (1-(fG+fC))/2  K(1-(fG+fC))/2 (1-(fG+fC))/2 |
			/// Q = | (fG+fC)/2     *              (fG+fC)/2      K(fG+fC)/2    |
			///     | K(fG+fC)/2    (fG+fC)/2      *              (fG+fC)/2     |
			///     | (1-(fG+fC))/2 K(1-(fG+fC))/2 (1-(fG+fC))/2  *             |
			//////////
		break;
		case NUC_TN93eq:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case NUC_TN93:
			//////////
			/// INPUT: kappa, kappa2;
			///     | *    fA   K1fA fA   |
			/// Q = | fC   *    fC   K2fC |
			///     | K1fG fG   *    fG   |
			///     | fT   K2fT fT   *    |
			//////////
			kappa1 = kappa.at(0);
			kappa2 = kappa.at(1);
			a = c = 1; b = kappa1;		// *fA
			a_ = d = 1; e = kappa2;		// *fC
			b_ = kappa1; d_ = f = 1;		// *fG
			c_ = f_ = 1; e_ = kappa2;	// *fT
		break;
		case NUC_TIMeq:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case NUC_TIM:
			//////////
			/// INPUT: alpha, beta, gamma
			///     | *   fA  AfA  BfA |
			/// Q = | fC  *   BfC  GfC |
			///     | AfG BfG *    fG  |
			///     | BfT GfT fT   *   |
			///////////
			alpha = kappa.at(0);
			beta = kappa.at(1);
			gamma = kappa.at(2);
			a = 1; b = alpha; c = beta;		// *fA
			a_ = 1; d = beta; e = gamma;	// *fC
			b_ = alpha; d_ = beta; f = 1;	// *fG
			c_ = beta; e_ = gamma; f_ = 1;	// *fT
		break;
		case NUC_TVMeq:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
		case NUC_TVM:
			//////////
			/// INPUT: alpha, beta, gamma, delta
			///     | *   fA  AfA  BfA |
			/// Q = | fC  *   GfC  AfC |
			///     | AfG GfG *    DfG  |
			///     | BfT AfT DfT  *   |
			///////////
			alpha = kappa.at(0);
			beta = kappa.at(1);
			gamma = kappa.at(2);
			delta = kappa.at(3);
			a = 1; b = alpha; c = beta;			// *fA
			a_ = 1; c = gamma; d = alpha;		// *fC
			b_ = alpha; d_ = gamma; f = delta;	// *fG
			c_ = beta; e_ = alpha; f_ = delta;	// *fT
		break;
		case NUC_SYM:
			pi.at(A) = pi.at(C) = pi.at(G) = pi.at(T) = 0.25;
			/// SYM and GTR are same, except SYM has equal nucleotide frequencies. No break.
		case NUC_GTR:
			//////////
			/// INPUT: relative rates abcdef
			///     | *     rACfC rAGfA rATfA |
			/// Q = | rACfC *     rCGfG rCTfC |
			///     | rAGfG rCGfG *     rGTfG |
			///     | rATfT rCTfT rGTfT *     |
			//////////
			a  = Sij.at(0); b  = Sij.at(1); c  = Sij.at(2);		// *fA
			a_ = Sij.at(0); d  = Sij.at(3); e  = Sij.at(4);		// *fC
			b_ = Sij.at(1); d_ = Sij.at(3); f  = Sij.at(5);		// *fG
			c_ = Sij.at(2); e_ = Sij.at(4); d_ = Sij.at(5);		// *fT
		break;
		case AA_JTT:
			SetRelativeRates(jttRelativeRates); 
			if (!freqSet) SetFrequencies(jttFrequencies);
			break;
		case AA_WAG: 
			SetRelativeRates(wagRelativeRates); 
			if (!freqSet) SetFrequencies(wagFrequencies);
			break;
		case AA_DAYHOFF78: 
			SetRelativeRates(dayhoffRelativeRates); 
			if (!freqSet) SetFrequencies(dayhoffFrequencies);
			break;
		case AA_BLOSUM62: 
			SetRelativeRates(blosumRelativeRates); 
			if (!freqSet) SetFrequencies(blosumFrequencies);
			break;
		case AA_MTREV24: 
			SetRelativeRates(mtrevRelativeRates); 
			if (!freqSet) SetFrequencies(mtrevFrequencies);
			break;
		case AA_CPREV45: 
			SetRelativeRates(cprevRelativeRates); 
			if (!freqSet) SetFrequencies(cprevFrequencies);
			break;
	default:
		cerr << "Do not know what to do with your model." << endl;
		exit(EXIT_FAILURE);
	}

	if (numStates == NUM_NUC) {
		Sij.at(0) = a;
		Sij.at(1) = b;
		Sij.at(2) = c;
		Sij.at(3) = d;
		Sij.at(4) = e;
		Sij.at(5) = f;
	}

	SetupMatrix(true);
}

void RateMatrix::SetRelativeRates(
								  double *inRelativeRate
								 ) 
{
	int i;
	for (i=0; i<NUM_AA_REL_RATES; i++) 
		Sij.at(i) = inRelativeRate[i];
}

void RateMatrix::SetFrequencies(
								double *inFrequencies, 
								inClade *branch
							   )
{
	int i;
	if (branch == NULL) {
		for (i=0; i<numStates; i++)
			pi.at(i) = inFrequencies[i];
	} else {
		for (i=0; i<numStates; i++)
			branch->values2Export2Freq.at(i) = inFrequencies[i];
	}
}

//////////
// Ensures that frequencies are not smaller than MINFREQ and
// that two frequencies differ by at least 2*MINFDIFF.
// This avoids potential problems later when eigenvalues
// are computed.
//////////
void RateMatrix::CheckFrequencies()
{
	int i, j;
	double diff;
	
	// required frequency difference
	double MINFDIFF = 1.0E-10;

	// lower limit on frequency
	double MINFREQ = 1.0E-10;

	int maxi = 0;
	double sum = 0.0;
	double maxfreq = 0.0;
	for (i = 0; i < numStates; i++) {
		if (pi.at(i) < MINFREQ) pi.at(i) = MINFREQ;
		if (pi.at(i) > maxfreq) {
			maxfreq = pi.at(i);
			maxi = i;
		}
		sum += pi.at(i);
	}
	
	diff = 1.0 - sum;
	pi.at(maxi) += diff;

	for (i = 0; i < numStates - 1; i++) {
		for (j = i+1; j < numStates; j++) {
			if (pi.at(i) == pi.at(j)) {
				pi.at(i) += MINFDIFF;
				pi.at(j) -= MINFDIFF;
			}
		}
	}
}

void RateMatrix::CheckInputAAFrequencies(
										 vector<double>& inFreq
										)
{
	int i, j;
	double diff;
	
	// required frequency difference
	double MINFDIFF = 1.0E-10;

	// lower limit on frequency
	double MINFREQ = 1.0E-10;

	int maxi = 0;
	double sum = 0.0;
	double maxfreq = 0.0;
	for (i = 0; i < NUM_AA; i++) {
		if (inFreq.at(i) < MINFREQ) inFreq.at(i) = MINFREQ;
		if (inFreq.at(i) > maxfreq) {
			maxfreq = inFreq.at(i);
			maxi = i;
		}
		sum += inFreq.at(i);
	}
	
	diff = 1.0 - sum;
	inFreq.at(maxi) += diff;

	for (i = 0; i < NUM_AA - 1; i++) {
		for (j = i+1; j < NUM_AA; j++) {
			if (inFreq.at(i) == inFreq.at(j)) {
				inFreq.at(i) += MINFDIFF;
				inFreq.at(j) -= MINFDIFF;
			}
		}
	}
}

void RateMatrix::setPij(
						Site site, 
						double BL, 
						int rateHetero
					   )
{
    switch (rateHetero) {
	case GammaRates:
		BL *= site.returnGamma();
	case NoRates:
		SetMatrix(Pij.at(0), BL);
	    break;
	case DiscreteGammaRates:
	    for (int j=0; j < num_categories; j++)
	    	SetMatrix(Pij.at(j), catRate.at(j)*BL);
	    break;
	case CodonRates:
//	    for (j=0; j < 3; j++)
//	    	SetMatrix(matrix[j], des->nodeEnv->catRate[j]*inlen);
//		for (i = 0; i < sequence_size; i++) {
//			if (random_placement)
//				pos = (double)rndu() * (des->seq_evo.size()-1);
//			else pos = i;
//			vector<Site>::iterator posit = des->seq_evo.begin()+pos;
//			cat=(pos+iTree->codon_offset)%3;
//			calling_routine.assign("MutateSequence->GONE, CodonRates");
//			result=SetState(
//							matrix[cat]+((*posit).returnState() * numStates), 
//							calling_routine
//						   );
//			if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) 
//				do_sub=true;
//			else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) )
//				do_sub = true;
//			else rejected++;
//
//			if (do_sub) {
//				// Check if stop codon is result of substitution.
//				int codon[3];
//				if (pos - (3-iTree->codon_offset) < 0 || pos + (3-iTree->codon_offset) >= des->seq_evo.size()) {
//					; // Punt, since I do not keep track of other partitions (what about subsequences??)
//				} else if ( (pos+iTree->codon_offset) % 3 == 0 ) {
//					codon[0] = result;
//					codon[1] = (*(posit+1)).returnState();
//					codon[2] = (*(posit+2)).returnState();
//				} else if ( (pos+iTree->codon_offset) % 3 == 1 ) {
//					codon[0] = (*(posit-1)).returnState();
//					codon[1] = result;
//					codon[2] = (*(posit+1)).returnState();
//				} else if ( (pos+iTree->codon_offset) % 3 == 2 ) {
//					codon[0] = (*(posit-2)).returnState();
//					codon[1] = (*(posit-1)).returnState();
//					codon[2] = result;
//				} else {
//					cerr << "Huh?? how can there be more than 3 codon positions?" << endl;
//					exit(EXIT_FAILURE);
//				}
//				if ( !Stop_Codon(codon) ) {
//					(*posit).setState(result);
//				} else rejected++;
//			}
//	    }
//
		break;
	}
}

void RateMatrix::PrintTransitionMatrices()
{
	cerr << "Qij:";
	printQij();

	cerr << "Pij:" << endl;
	cerr << num_categories << endl;
	printPij();
}

void RateMatrix::printPij()
{
	int i = 0;
	
	for (int cat = 0; cat < num_categories; cat++) {
		cerr << "  Category " << cat;
		for (vector<double>::iterator it = Pij.at(cat).begin(); it != Pij.at(cat).end(); ++it, i++) {
			if (i % numStates == 0) cerr << endl << "    ";
			cerr << (*it) << " ";
		}
		cerr << endl;
	}
	cerr << endl;

}

void RateMatrix::printQij()
{
	int i = 0;
	for (vector<double>::iterator it = Qij.begin(); it != Qij.end(); ++it, i++) {
		if (i % numStates == 0) cerr << endl;
		cerr << (*it) << " ";
	}
	cerr << endl;
}

void RateMatrix::fullReport()
{
	PrintTransitionMatrices();
	printCategories();
}

void RateMatrix::printCategories()
{
	size_t i;

	cerr << "Gamma categories: " << num_categories << "  alpha: " << alphaGamma << endl;
	for (i = 0; i < num_categories; i++) {
		cerr << i+1 << ":" << catRate.at(i) << "  ";
	}
	cerr << endl;

}

void RateMatrix::SetupMatrix(bool transition_probability_independent_sites_setup)
{
	unsigned int i,j,k;
	double mr;
	double sum;
	double U[numStates_squared+1], V[numStates_squared+1], T1[numStates_squared+1], T2[numStates_squared+1];
	double Qij2pass[numStates_squared+1];

	//////////
	/// transition_probability_independent_sites_setup:
	/// If we are calculating the transition probabilities as independent sites models, then we need
	/// to run the code in the if statements. However, if we are doing the pseudo-dependent sites
	/// calculations, we want to skip those code segments. 
	/// Segment #1: 
	///   * Sets the Qij entries to be the input substitution model. 
	///     -> For pseudo-dependent sites, we excise the Q matrix from the 4^N X 4^N sequence model.
	/// Segment #2:
	///	  * Sets Qij so that branch length corresponds to expected number of changes.
	///     -> pseudo-dependent sites branch lengths correspond not to 1 change per site, but whatever
	///        the dependent sites model constrains it to.
	//////////
	CheckFrequencies();
	k=0;
	if (transition_probability_independent_sites_setup) {
		for (i=0; i<numStates-1; i++) {
			for (j=i+1; j<numStates; j++) {
				Qij.at(i*numStates+j) = Qij.at(j*numStates+i) = Sij.at(k++);
			}
		}

		for (i=0; i<numStates; i++) {
			for (j=0; j<numStates; j++) { 
				Qij.at(i*numStates+j) *= pi.at(j);
			}
		}

		mr=0;		
		for (i=0; i<numStates; i++) {
			sum = 0;
			Qij.at(i*numStates+i)=0; 
			for (j=0; j<numStates; j++) { 
				sum += Qij.at(i*numStates+j);
			}
			Qij.at(i*numStates+i) = -sum; 
			mr += pi.at(i) * sum;
		}

		//////////
		/// Multiplies every element of Qij by 1.0/mr. This is abyx.
		//////////
		for (i=0; i<numStates_squared; Qij.at(i)*=1.0/mr,i++) ; 
	}

	// Testcode: Print Qij.
	//if (!transition_probability_independent_sites_setup) printQij();
	
	//Testcode: does sum = 1?
	sum = 0;
	double partial = 0;
	for (i = 0; i < numStates; i++) {
		partial = 0;
		for (j = 0; j < numStates; j++) {
			if (i != j)
				partial += Qij.at(i*numStates+j);
		}
		sum += pi.at(i) * partial;
	}
	//if (!transition_probability_independent_sites_setup) cerr << "Sum: " << sum << endl;

	for (i = 0; i < numStates_squared; i++) Qij2pass[i]=Qij.at(i);
	if ((k=eigen(1, Qij2pass, numStates, Root, T1, U, V, T2))!=0) {
		fprintf(stderr, "\ncomplex roots in SetupMatrix");
		exit(EXIT_FAILURE);
	}
	xtoy (U, V, numStates_squared);
	matinv (V, numStates, numStates, T1);
	vector<unsigned int> squared (numStates, 0), power_1 (numStates, 0);
	for (i=0; i<numStates; i++) { squared[i]=i*numStates_squared; power_1[i]=i*numStates; }
	for (i=0; i<numStates; i++) {
   		for (j=0; j<numStates; j++) {
   			for (k=0; k<numStates; k++) {
				Cijk.at(squared[i]+power_1[j]+k) = U[power_1[i]+k]*V[power_1[k]+j];
  			}
  		}
   	}
}

void RateMatrix::SetMatrix(
						   vector<double>& matrix, 
						   double len
						  )
{
	int i,j,k;
	double expt[numStates];
	vector<double>::iterator P;
	
	// P(t)ij = SUM Cijk * exp{Root*t}

	P=matrix.begin();

	if (len<1e-6) { 
		for (i=0; i<numStates; i++) {
			for (j=0; j<numStates; j++) {
				if (i==j) 
					(*P)=1.0;
				else 
					(*P)=0.0;
				P++;
			}
		}
		return; 
	}
	
	for (k=1; k<numStates; k++) {
		expt[k]=exp(len*Root[k]);
	}

	vector<unsigned int> squared (numStates, 0), power_1 (numStates, 0);
	for (i=0; i<numStates; i++) { squared[i]=i*numStates_squared; power_1[i]=i*numStates; }
	for (i=0; i<numStates; i++) {
		for (j=0; j<numStates; j++) {
			(*P)=Cijk.at(squared[i]+power_1[j]+0);
			for (k=1; k<numStates; k++) {
				(*P)+=Cijk.at(squared[i]+power_1[j]+k)*expt[k];
			}
			P++;
		}
	}
}

void RateMatrix::SetAAFreqs(int theModel, inClade *branch) 
{
	switch (theModel) {
		case AA_JTT: SetFrequencies(jttFrequencies, branch); break;
		case AA_WAG: SetFrequencies(wagFrequencies, branch); break;
		case AA_DAYHOFF78: SetFrequencies(dayhoffFrequencies, branch); break;
		case AA_BLOSUM62: SetFrequencies(blosumFrequencies, branch); break;
		case AA_MTREV24: SetFrequencies(mtrevFrequencies, branch); break;
		case AA_CPREV45: SetFrequencies(cprevFrequencies, branch); break;
		default: cerr << "You cannot set the general matrix in lineages." << endl; exit(EXIT_FAILURE);
	}
}

void RateMatrix::printCijk()
{
	cerr << "Cijk: " << endl;
	for (vector<double>::iterator it = Cijk.begin(); it != Cijk.end(); it++)
		cerr << (*it) << " ";
	cerr << endl;
}

void RateMatrix::printRoot()
{
	cerr << "Root: " << endl;
	for (int i = 0; i < NUM_AA; i++) {
		cerr << Root[i] << " ";
	}
	cerr << endl;
}

