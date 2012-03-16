#include "model.h"

char *modelNames[numModels]={
	"HKY",
	"F84",
	"GTR",
	"JTT",
	"WAG",
	"PAM",
	"BLOSUM",
	"MTREV",
	"CPREV",
	"GENERAL"
};

char *modelTitles[numModels]={
	"HKY: Hasegawa, Kishino & Yano (1985)",
	"F84: Felsenstein (1984)",
	"GTR: General time reversible (nucleotides)",
	"JTT: Jones, Taylor & Thornton (1992) CABIOS  8:275-282\n             DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
	"WAG: Whelan & Goldman (2001) Mol Biol Evol 18:691Ð699",
	"PAM: Dayhoff, Schwartz & Orcutt (1978)\n             DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
	"BLOSUM62: Henikoff & Henikoff (1992) PNAS USA 89:10915-10919",
	"MTREV24: Adachi & Hasegawa (1996) J Mol Evol 42:459-468",
	"CPREV45: Adachi et al. (2000) J Mol Evol 50:348-358",
	"GENERAL: General time reversible (amino acids)"
};

int model, numStates, isNucModel, userFreqs, equalFreqs;

std::string stateCharacters;

double *freq, *addFreq;

int  FindModel(string theModel) 
{
	int model = -1;

	trim(theModel);
	for(size_t i = 0; i < theModel.size(); i++) {
		theModel.at(i) = toupper(theModel.at(i));
	}
	for(int i = F84; i<numModels; i++) {
		if(theModel.compare(0,3,modelNames[i], 0, 3) == 0) {
			model = i;
			if (model <= GTR) {
				if(isNucModel == -1) {
					isNucModel = 1;
					numStates = 4;
				} else {
					if(!isNucModel) {
						cerr << "Substitution model " << modelNames[model] << " is a nucleotide model, but simulation run is for amino acids." << endl;
						exit(EXIT_FAILURE);							
					}
				}
			} else {
				if(isNucModel == -1) {
					isNucModel = 0;
					numStates = 20;
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
				numStates = 4;
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

	return model;
}

void SetFrequencies(std::string frequencies) 
{
	list<string> arg_split;
	int j;
	double sumfreq = 0;

	if(isNucModel) {
		if(toupper(frequencies.at(0)) == 'E') ;
		else {
			equalFreqs=0;
			arg_split = split(frequencies, ",");
			j = 0;
			for(list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++) {
				sumfreq += atof((*it).c_str());
				nucFreq[j] = atof((*it).c_str());
				j++;
			}
			if(j != NUM_NUC) {
				cerr << "Bad Nucleotide Frequencies: " << frequencies << endl << endl;
				exit(EXIT_FAILURE);
			}
			j = 0;
			for(list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++,j++) {
				nucFreq[j] /= sumfreq;
				cerr << "NUCFREQ[" << j << "] = " << nucFreq[j] << endl;
			}
		}
	} else {
		aaFreqSet = 1;
		if(toupper(frequencies.at(0)) == 'E') {
			equalFreqs = 1;
			for(j=0; j<NUM_AA; j++) {
				aaFreq[j] = 0.05;
			}
		} else {
			equalFreqs = 0;
			arg_split = split(frequencies, ",");
			j = 0;
			for(list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++) {
				sumfreq += atof((*it).c_str());
				aaFreq[j] = atof((*it).c_str());
				j++;
			}
			if(j != NUM_AA) {
				cerr << "Bad amino acid frequencies: " << frequencies << endl << endl;
				exit(EXIT_FAILURE);
			}
			j = 0;
			for(list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++,j++) {
				aaFreq[j] /= sumfreq;
			}
		}
	}
}

void SetModelFreqs(int theModel, inClade *branch) 
{
	if (isNucModel) {
		SetNucFreqs(theModel, branch);
	} else {
		SetAAFreqs(theModel - numNucModels, branch);
	}
}

void SetModel(int theModel, inClade *branch)
{	
	int i;
	
	model=theModel;

	if (isNucModel) {
		numStates = NUM_NUC;

		if (branch != NULL) {
			if ( !(branch->values2Export2Freq.empty()) )
				for (int j = 0; j < numStates; j++)
					nucFreq[j] = branch->values2Export2Freq.at(j);
		}
		
		SetNucModel(theModel);

		freq = nucFreq;
		addFreq = nucAddFreq;
		stateCharacters = nucleotides;
	} else {
		numStates = NUM_AA;

		if (branch != NULL) {
			if ( !(branch->values2Export2Freq.empty()) ) {
				aaFreqSet = 1;
				for (int j = 0; j < numStates; j++)
					aaFreq[j] = branch->values2Export2Freq.at(j);		
			}
		}
		
		SetAAModel(theModel - numNucModels);

		freq = aaFreq;
		addFreq = aaAddFreq;
		stateCharacters = aminoAcids;
	}
	
	addFreq[0]=freq[0];
	for (i = 1; i < numStates; i++) {
		addFreq[i] = addFreq[i-1] + freq[i];
	}
}

void SetMatrix(double *matrix, double len)
{
	if (isNucModel) {
		SetNucMatrix(matrix, len);
	} else {
		SetAAMatrix(matrix, len);
	}
}

void SetVector(double *vector, short state, double len)
{
	if (isNucModel) {
		SetNucVector(vector, state, len);
	} else {
		SetAAVector(vector, state, len);
	}
}

