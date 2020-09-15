template <typename Comp_Type> void cellField<Comp_Type>::interact() {
    Comp_Type RTKa, RTKi, PTPac, PTPic, PTPam, PTPim;

    double* p=static_cast<double*>(mxGetData(mxGetField(Params, 0, "IntPar")));
//Numbering of parameters in C starts at 0, in matlab starts at 1 !!!
    
    for (int i=0; i<CellsSize; ++i) {
        Comp_Type &RTKa_old =Cells[i];             RTKa = RTKa_old;
        Comp_Type &RTKi_old =Cells[i+  CellsSize]; RTKi = RTKi_old;
        Comp_Type &PTPam_old=Cells[i+2*CellsSize]; PTPam=PTPam_old;
        Comp_Type &PTPim_old=Cells[i+3*CellsSize]; PTPim=PTPim_old;
        Comp_Type &PTPac_old=Cells[i+4*CellsSize]; PTPac=PTPac_old;
        Comp_Type &PTPic_old=Cells[i+5*CellsSize]; PTPic=PTPic_old;

        if (State[i]==127) {

            //initialize parameters sent by MATLAB
            
            double k0=p[0];
            double k1=p[1];
            double k2=p[2];
            double k3=p[3];
            double k4=p[4];
            double k5=p[5];
            double k6=p[6];
            double k7=p[7];
            double k2p=p[8];
            
            double Km0=p[10];
            double Km1=p[11];
            double Km2=p[12];
            double Km3=p[13];
            double Km4=p[14];
            double Km5=p[15];
            double Km6=p[16];
            double Km7=p[17];
            double Km2p=p[18];
            double Hill=p[19];
            double Rt=p[20];
            double not_used=p[22];
            
            double krand1=p[23];
            double krand2=p[24];
            double krand3=p[25];
            
            //interface with original implementation 
            //based on receptor tyrosine kinase simulations

            double Ga=RTKa;
            double Gi=RTKi;
            
            double Ma=PTPam;
            double Mi=PTPim;
            
            double Ra=PTPac;
            double Ri=PTPic;
            
            //normal distribute random number Box-Muller Method: https://en.wikipedia.org/wiki/Normal_distribution 
            
            double rand_u=(double)rand()/((double)RAND_MAX+1.);
            double rand_v=(double)rand()/((double)RAND_MAX+1.);
            double rand_x1=sqrt(-2*log(rand_u))*cos(2*M_PI*rand_v);
            
            rand_u=(double)rand()/((double)RAND_MAX+1.);
            rand_v=(double)rand()/((double)RAND_MAX+1.);
            double rand_x2=sqrt(-2*log(rand_u))*cos(2*M_PI*rand_v);
            
                        rand_u=(double)rand()/((double)RAND_MAX+1.);
            rand_v=(double)rand()/((double)RAND_MAX+1.);
            double rand_x3=sqrt(-2*log(rand_u))*cos(2*M_PI*rand_v);
            
            //calculate dG/dt

            exchangeConc(
                    k3*Gi*Ra
                    -k4*Ga*Ma
                    +krand1*rand_x1,
                    RTKi_old,RTKa_old);
           
            //calculate dM/dt

            exchangeConc(
                    ((k5*Ra*pow(Mi,Hill))/(Km5+pow(Mi,Hill)))
                    +((k7*Mi)/(Km7+Mi))
                    -((k6*Ma)/(Km6+Ma))
                    +krand2*rand_x2,
                    PTPim_old,PTPam_old);
            
            //calculate dR/dt
            
            exchangeConc(
                    ((k0*Ga*Ri)/(Km1+Ri))
                    -((k2*Ra)/(Km2+Ra))
                    +krand3*rand_x3,
                    PTPic_old,PTPac_old);

            //add constant increase in G (not used)
            RTKi_old=RTKi_old+p[21];
        }
        
        

        if (State[i]==255) {
        }
        

        if (State[i]==1) {
        }
        

        if (State[i]==127) {
        }
        

        if (State[i]==95) {
        }
        

        if (State[i]==63) {
        }
    }
}

template <typename Comp_Type> void cellField<Comp_Type>::exchangeConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere) {
  if (fromWhere<howMuch || toWhere<-howMuch) {
     if (howMuch>0){
          howMuch=fromWhere;
      }
      else 
      {
          howMuch=-toWhere;
      }
  }
  fromWhere-=howMuch;toWhere+=howMuch;
}

template <typename Comp_Type> void cellField<Comp_Type>::oneWayConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere) {
  if (howMuch>fromWhere) howMuch=fromWhere;
  toWhere+=howMuch;
  fromWhere-=howMuch; 
}

template <typename Comp_Type> void cellField<Comp_Type>::dissociateConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere, Comp_Type &andWhere) {
  if (howMuch>fromWhere) howMuch=fromWhere;
  toWhere+=howMuch;
  andWhere+=howMuch;
  fromWhere-=howMuch; 
}

template <typename Comp_Type> void cellField<Comp_Type>::associateConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andWhere, Comp_Type &toWhere) {
  if (howMuch>fromWhere) howMuch=fromWhere;
  if (howMuch>andWhere) howMuch=andWhere;
  toWhere+=howMuch;
  fromWhere-=howMuch; 
  andWhere-=howMuch;
}

template <typename Comp_Type> void cellField<Comp_Type>::asDissConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andWhere, Comp_Type &toWhere) {
// A+B <-> AB
  if (fromWhere<howMuch || andWhere<howMuch || toWhere<-howMuch) {
      if (howMuch>0) howMuch=MIN(fromWhere,andWhere);
      else howMuch=-toWhere;
  }
  toWhere+=howMuch;
  andWhere-=howMuch;
  fromWhere-=howMuch; 
}


template <typename Comp_Type> void cellField<Comp_Type>::asymAsDissConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andFromWhere, Comp_Type &toWhere, Comp_Type &andToWhere) {
// A+B -> AB -> A+B*
  if (fromWhere<howMuch || andFromWhere<howMuch || toWhere<-howMuch) {
      if (howMuch>0) howMuch=MIN(fromWhere,andFromWhere);
      else howMuch=-toWhere;
  }
  toWhere+=howMuch;
  if (howMuch<0) andToWhere-=howMuch;
  if (howMuch>0) andFromWhere-=howMuch;
  fromWhere-=howMuch; 
}

