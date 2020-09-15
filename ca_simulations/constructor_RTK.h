template <typename Comp_Type> cellField<Comp_Type>::cellField(int* _size) : SubSize(7) {
	size[0]=_size[0];
	size[1]=_size[1];
	size[2]=_size[2];
	CellsSize=size[0]*size[1]*size[2];
	numStates=1;
    numComparts=1;
	SpeciesSize=CellsSize*numSpecies;
    maxRecip=new int[numSpecies];
	for (int i=0; i<DIM; i++) {
		anIso_pos[i]=new int[numStates];
		anIso_neg[i]=new int[numStates];
        anIso_pos[i][0]=((size[i]>1)?1:0);
        anIso_neg[i][0]=((size[i]>1)?1:0);
    }
	Cells = new Comp_Type[SpeciesSize];
	State = new int[CellsSize];
    Center = new vec3d[CellsSize];

#ifdef _MATLAB_MEX
    double stoFactor=*(static_cast<double*>(mxGetData(mxGetField(Params, 0, "StochVal"))));
#else
    double stoFactor=0.0;
#endif
    
    //srand(19);

    for (int i=0; i<CellsSize; i++) {
        vec3d tCenter(((size[0]>1)?(0.5-rand()*rand_norm):0.),\
                      ((size[1]>1)?(0.5-rand()*rand_norm):0.),\
                      ((size[2]>1)?(0.5-rand()*rand_norm):0.));

        Center[i]=vec3d(0.5,0.5,0.5)+stoFactor*tCenter;
    }
    
	initSubGrid();

    int Pos=0;
    PMArea=0;
    ERArea=0;
    REArea=0;
    cytArea=0;
    intArea=0;

#ifdef _MATLAB_MEX
    double *DiffRad=static_cast<double*>(mxGetData(mxGetField(Params, 0, "DiffRad")));
#else
    double DiffRad[1]={12.0};
#endif
    
        
	for (int z=0;z<size[2];z++)
		for (int y=0;y<size[1];y++)
			for (int x=0;x<size[0];x++) {
                State[Pos]=imageData[Pos];
                Cells[Pos]=initVals[Pos];

                if (State[Pos]==1) PMArea++;
                else if (State[Pos]==63) REArea++;
                else if (State[Pos]==95) ERArea++;
                else if (State[Pos]==127) ERArea++;
                else if (State[Pos]==255) cytArea++;
                
                if (State[Pos]>0) intArea++;
                Pos++;
			}

    for (;Pos<SpeciesSize;Pos++) {
        Cells[Pos]=initVals[Pos];
    }

    initNeighb(DiffRad);
};
