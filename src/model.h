/*  Header file for model.c                                */
#ifndef _MODEL_H_
#define _MODEL_H_

#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cstring>
#include <math.h>
#include "evolve.h"
#include "inTree.h"
#include "nucmodels.h"
#include "aamodels.h"

class inClade;

extern std::string stateCharacters;

enum { NONE=-1, F84, HKY, GTR, JTT, WAG, PAM, BLOSUM, MTREV, CPREV, GENERAL, numModels };

extern char *modelNames[numModels];
extern char *modelTitles[numModels];

extern int model, isNucModel, numStates, userFreqs, equalFreqs;

extern double *freq, *addFreq;

void SetFrequencies(std::string frequencies);
int  FindModel(std::string theModel);
void SetModel(int theModel, inClade *branch);
void SetMatrix(double *matrix, double len);
void SetVector(double *vector, short state, double len);
void SetModelFreqs(int theModel, inClade *branch);

#endif /* _MODEL_H_ */
