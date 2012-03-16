// Copyright (C) 2008 Cory Strope <corystrope@gmail.com>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

#ifndef _TREEFILE_H_
#define _TREEFILE_H_

#include <cstring>
#include <iostream>
#include <exception>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <ctype.h>
#include <string>
#include <list>
#include <iostream>
#include <iomanip>
#include "tree.h"
#include "inTree.h"
#include "twister.h"
#include "evolve.h"

using namespace std;

extern int invariableSites;
extern size_t TNode_InheritMotifSites;

void initializeBranchRun(inTree *iTree, TNode *des, double *scale, double *invar_scale, int inNumSites);
void initializeCladeParameters(TTree *tree, TNode *node);
void initializeTreeParameters(inTree *iTree, double *scale, double *invar_scale);
TNode *NewNode(TTree *tree);
char ReadToNextChar(string tree_str, int *pos);
void ReadUntil(string tree_str, int *pos, char stopChar, char *what);

#endif /* _TREEFILE_H_ */

