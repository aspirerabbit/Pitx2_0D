//const double ach = 0.1; 
//const double iso=1;
//const double Flecainide=2;
//double Flecainide_kon= Flecainide*5.5;//Flecainide*5500;
//double Flecainide_kon= 0.011;
//double Flecainide_DO =0;
//double dFlecainide_DO =0;
// gkr = 0.0294*sqrt(ko/5.4)*(1.0/(1.0+Flecainide/1.5));
 
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define bcl 1000     /* Basic cycle length (ms) */
#define beats 100     /* Number of Beats */
#define S2 1000

//#define Vc				0.016404 //111111111
//#define Vsr				0.001094//11111111111
//#define Vss				0.00005468//11111111



/* List of variables and paramaters (this code uses all global variables) */

        /* Creation of Data File */
        FILE *ap3;//changed by Bai
		char * str3;//changed by Bai
		char * strfmax;//changed by Bai
        FILE *fmaxs;
       // FILE *fpara;
        void prttofile ();      
        int printdata;  
        int printval;   
		
		const double ach = 0.0; /* Acetylcholine concentration */

		const double iso=0;//ISO
		
		const double AF=0;

		const double Flecainide=0;//2;

		const double PITX2=1;
		
		double LA=1;
		double RA=0;

		double Vc,Vsr,Vss;

			//Calcium buffering dynamics
	double Bufc=0.2;//11111111111111111
	double Kbufc=0.001;//1111111111111111
	double Bufsr=10.;//111111111111111111
	double Kbufsr=0.3;//111111111111111111
	double Bufss=0.4;/////////////////11111111111111111
	double Kbufss=0.00025;////////////11111111111111111111111

	//Intracellular calcium flux dynamics
	//double Vmaxup=0.006375*1.3;//11111111111111 
	double Vmaxup=0.006375;//*2;//11111111111111 
	double Kup=0.00025*(2-1*iso);
	double Vrel=0.102;////////////////11111111111111111111  

	double k1_=0.15;/////////////////111111111111111111111111
	double k2_=0.045;///////////////1111111111111111111111
	double k3=0.060;////////////////1111111111111111111111
	double k4=0.005;//0.000015;///////11111111111111111111111
	double EC=1.5;//////////////////////1111111111111111111111
	double maxsr=2.5;/////////////////1111111111111111111111
	double minsr=1.;//////////////////111111111111111111
	double Vleak=0.00036;/////////////1111111111111111111111
	double Vxfer=0.0038;//////////////111111111111111111111111111
	
	double k1=0;////////////////////////1111111111
	double k2=0;///////////////////////11111111111111
	double kCaSR=0;////////////////////111111111111111

	//double dNai=0;//1111111111111
	//double dKi=0;//1111111111111
	//double dCai=0;//1111111111111
	double dCaSR=0;//111111111111
	double dCaSS=0;///////////////////////111111111
	double dRR=0;/////////////////////////111111111
	
	double CaCSQN=0;//1111111111111
	double bjsr=0;//11111111111111
	double cjsr=0;//111111111111111
	double CaSSBuf=0;////////////////////1111111111111
	double CaBuf=0;//11111111111
	double bc=0;//1111111111
	double cc=0;//111111111
	double CAPACITANCE=0.185;
	
	double bcss;///////////////////////1111111111
	double ccss;////////////////////////1111111111
	
	double Irel,Ileak,Irelease,Iup,Ixfer;
	
	
	double Flecainide_IC=2;
    double Flecainide_kon= Flecainide*0.0055;//5.5;//Flecainide*5500;
    double Flecainide_koff= Flecainide_IC*0.0055;
    double Flecainide_DO ;
    double dFlecainide_DO =0;
	


        /* Cell Geometry */
        const double l = 0.01;       /* Length of the cell (cm) */
        const double a = 0.0008;     /* Radius of the cell (cm) */
        const double pi = 3.141592;  /* Pi */
        double vcell;   /* Cell volume (uL) */
        double ageo;    /* Geometric membrane area (cm^2) */
        double acap;    /* Capacity */
        double vmyo;    /* Myoplasm volume (uL) */
        double vmito;   /* Mitochondria volume (uL) */
        double vsr;     /* SR volume (uL) */
        double vnsr;    /* NSR volume (uL) */
        double vjsr;    /* JSR volume (uL) */
        
        /* Voltage */
        double v;       /* Membrane voltage (mV) */
        double vnew;    /* New Voltage (mV) */
        double dvdt;    /* Change in Voltage / Change in Time (mV/ms) */
        double dvdtnew; /* New dv/dt (mV/ms) */
        double boolien; /* Boolien condition to test for dvdtmax */
        
        /* Time Step */
        double dt;      /* Time step (ms) */
        double t;       /* Time (ms) */
        double udt;     /* Universal Time Step */
        int utsc;       /* Universal Time Step Counter */
        int nxstep;     /* Interval Between Calculating Ion Currents */
        int steps;      /* Number of Steps */
        int increment;  /* Loop Control Variable */ 
        
        /* Action Potential Duration and Max. Info */
        double vmax[beats+1];           /* Max. Voltage (mV) */
        double dvdtmax[beats+1];        /* Max. dv/dt (mV/ms) */
        double apd[beats+1];            /* Action Potential Duration */
        double toneapd[beats+1];        /* Time of dv/dt Max. */
        double ttwoapd[beats+1];        /* Time of 90% Repolarization */
        double trep[beats+1];           /* Time of Full Repolarization */
        double di[beats+1];             /* Diastolic Interval */
        double rmbp[beats+1];           /* Resting Membrane Potential */
        double nair[beats+1];           /* Intracellular Na At Rest */
        double cair[beats+1];           /* Intracellular Ca At Rest */
        double caimax[beats+1];         /* Peak Intracellular Ca */
        int i;                        /* Stimulus Counter */
        
        /* Total Current and Stimulus */
        double st;       /* Constant Stimulus (uA/cm^2) */
        double tstim;    /* Time Stimulus is Applied (ms) */
        double stimtime; /* Time period during which stimulus is applied (ms) */
        double it;       /* Total current (uA/cm^2) */

        /* Terms for Solution of Conductance and Reversal Potential */
        const double R = 8314;      /* Universal Gas Constant (J/kmol*K) */
        const double frdy = 96485;  /* Faraday's Constant (C/mol) */
        const double temp = 310;    /* Temperature (K) */

        /* Ion Valences */
        const double zna = 1;  /* Na valence */
        const double zk = 1;   /* K valence */
        const double zca = 2;  /* Ca valence */

        /* Ion Concentrations */
        double nai;    /* Intracellular Na Concentration (mM) */
        double nao;    /* Extracellular Na Concentration (mM) */
        double ki;     /* Intracellular K Concentration (mM) */
        double ko;     /* Extracellular K Concentration (mM) */
        double cai;    /* Intracellular Ca Concentration (mM) */
        double cao;    /* Extracellular Ca Concentration (mM) */
        double cmdn;   /* Calmodulin Buffered Ca Concentration (mM) */
        double trpn;   /* Troponin Buffered Ca Concentration (mM) */
        double nsr;    /* NSR Ca Concentration (mM) */
        double jsr;    /* JSR Ca Concentration (mM) */
        double csqn;   /* Calsequestrin Buffered Ca Concentration (mM) */
        double CaSS;  //TP06
		double CaSR; //TP06
		
		
        /* Myoplasmic Na Ion Concentration Changes */
        double naiont;  /* Total Na Ion Flow (mM/ms) */
        double dnai;    /* Change in Intracellular Na Concentration (mM) */

        /* Myoplasmic K Ion Concentration Changes */
        double kiont; /* Total K Ion Flow (mM/ms) */
        double dki;   /* Change in Intracellular K Concentration (mM) */

        /* NSR Ca Ion Concentration Changes */
        double dnsr;   /* Change in [Ca] in the NSR (mM) */
        double iup;    /* Ca uptake from myo. to NSR (mM/ms) */
        double ileak;  /* Ca leakage from NSR to myo. (mM/ms) */
        double kleak;  /* Rate constant of Ca leakage from NSR (ms^-1) */
        const double kmup = 0.00092;    /* Half-saturation concentration of iup (mM) */
        const double iupbar = 0.005;  /* Max. current through iup channel (mM/ms) */
        const double nsrbar = 15;       /* Max. [Ca] in NSR (mM) */
        
        /* JSR Ca Ion Concentration Changes */
        double djsr;                   /* Change in [Ca] in the JSR (mM) */
        double urel;                   /* Activation gate u of Ca release from jsr*/
		double urelss;                 /* Steady state of activation gate u*/
		double tauurel;                /* Time constant of activation gate u*/
		double vrel;                   /* Activation gate v of Ca release from jsr*/
		double vrelss;                 /* Steady state of activation gate v*/
		double tauvrel;                /* Time constant of activation gate v*/
		double wrel;                   /* Inactivation gate w of Ca release from jsr*/
		double wrelss;                 /* Steady state of inactivation gate w*/
		double tauwrel;                /* Time constant of inactivation gate w*/
		double fn;
		const double grelbarjsrol = 30; /* Rate constant of Ca release from JSR due to overload (ms^-1)*/
        double greljsrol;               /* Rate constant of Ca release from JSR due to CICR (ms^-1)*/
        double ireljsrol;               /* Ca release from JSR to myo. due to JSR overload (mM/ms)*/
        const double csqnbar = 10;      /* Max. [Ca] buffered in CSQN (mM)*/
        const double kmcsqn = 0.8;      /* Equalibrium constant of buffering for CSQN (mM)*/
//        double bjsr;                    /* b Variable for analytical computation of [Ca] in JSR (mM)*/
  //      double cjsr;                    /* c Variable for analytical computation of [Ca] in JSR (mM)*/
        double on;                      /* Time constant of activation of Ca release from JSR (ms)*/
        double off;                     /* Time constant of deactivation of Ca release from JSR (ms)*/
        double magrel;                  /* Magnitude of Ca release*/

        /* Translocation of Ca Ions from NSR to JSR */
        double itr;                /* Translocation current of Ca ions from NSR to JSR (mM/ms)*/
        const double tautr = 180;  /* Time constant of Ca transfer from NSR to JSR (ms)*/
        
        /* Myoplasmic Ca Ion Concentration Changes */
        double caiont;  /* Total Ca Ion Flow (mM/ms) */
        double dcai;    /* Change in myoplasmic Ca concentration (mM) */
        double b1cai;
		double b2cai;
		const double cmdnbar = 0.050;   /* Max. [Ca] buffered in CMDN (mM) */
        const double trpnbar = 0.070;   /* Max. [Ca] buffered in TRPN (mM) */
        const double kmcmdn = 0.00238;  /* Equalibrium constant of buffering for CMDN (mM) */
        const double kmtrpn = 0.0005;   /* Equalibrium constant of buffering for TRPN (mM) */

        /* Fast Sodium Current (time dependant) */
        double ina;    /* Fast Na Current (uA/uF) */
        double gna;    /* Max. Conductance of the Na Channel (mS/uF) */
        double ena;    /* Reversal Potential of Na (mV) */
        double ah;     /* Na alpha-h rate constant (ms^-1) */
        double bh;     /* Na beta-h rate constant (ms^-1) */
        double aj;     /* Na alpha-j rate constant (ms^-1) */
        double bj;     /* Na beta-j rate constant (ms^-1) */
        double am;     /* Na alpha-m rate constant (ms^-1) */
        double bm;     /* Na beta-m rate constant (ms^-1) */
        double h;      /* Na activation */
        double j;      /* Na inactivation */
        double m;      /* Na inactivation */
        double gB;
	    double f1;//TP06
	
        /* Current through L-type Ca Channel */
        double ilca;    /* Ca current through L-type Ca channel (uA/uF) */
        double ilcatot; /* Total current through the L-type Ca channel (uA/uF) */
        double ibarca;  /* Max. Ca current through Ca channel (uA/uF) */
        double d;       /* Voltage dependant activation gate */
        double dss;     /* Steady-state value of activation gate d  */
        double taud;    /* Time constant of gate d (ms^-1) */
        double f;       /* Voltage dependant inactivation gate */
        double fss;     /* Steady-state value of inactivation gate f */
        double tauf;    /* Time constant of gate f (ms^-1) */
        double fca;     /* Ca dependant inactivation gate */
        double taufca;  /* Time constant of gate fca (ms^-1) */
		double fcass;   /* Steady-state value of activation gate fca  */

	   double fss1,tauf1;//TP06
	   double sRR,sOO;

		
		
		
		const double gcalbar = 0.1238;
        
    	/* Acetylcholine-Activated Potassium Current */
		/* modified from Matsuoka et al., Jap J Physiol 2003;53:105-123 */
		double ikach; /* Acetylcholine-activated K current (uA/uF) */
		double gkach; /* Channel conductance of acetylcholine-activated K current (mS/uF) */
        double ekach; /* Reversal potential of acetylcholine-activated K current (mV) */
    	double alphayach; /* Alpha rate constant (ms^-1) */
	    double betayach; /* Beta rate constant (ms^-1) */
	    double tauyach; /* Time constant (ms) */
	    double yachss; /* Steady-state value */
	    double yach;
	    
		
	

        /* Ultra-Rapidly Activating Potassium Current */
        double ikur;   /* Ultra-rapidly activating K current (uA/uF) */
        double gkur;   /* Channel conductance of ultra-rapidly activating K current (mS/uF) */
        double ekur;   /* Reversal potential of ultra-rapidly activating K current (mV) */
        double uakur;    /* Ultra-rapidly activating K activation gate ua */
        double uakurss;  /* Steady-state value of activation gate ua */
        double tauuakur; /* Time constant of gate ua (ms^-1) */
        double alphauakur; /* Alpha rate constant of activation gate ua (ms^-1) */
		double betauakur;  /* Beta rate constant of activation gate ua (ms^-1) */
		double uikur;    /* Ultra-rapidly activating K activation gate ui*/
        double uikurss;  /* Steady-state value of activation gate ui */
        double tauuikur; /* Time constant of gate ui (ms) */
		double alphauikur; /* Alpha rate constant of activation gate ui (ms^-1) */
		double betauikur; /* Beta rate constant of activation gate ui (ms^-1) */

		double uskur;    /* Ultra-rapidly activating K activation gate ua *///BAI
        double uskurss;  /* Steady-state value of activation gate ua */     //BAI
        double tauuskur; /* Time constant of gate ua (ms^-1) */             //BAI

		/* Rapidly Activating Potassium Current */
        double ikr;   /* Rapidly activating K current (uA/uF) */
        double gkr;   /* Channel conductance of rapidly activating K current (mS/uF) */
        double ekr;   /* Reversal potential of rapidly activating K current (mV) */
        double xr;    /* Rapidly activating K time-dependant activation */
        double xrss;  /* Steady-state value of inactivation gate xr */
        double tauxr; /* Time constant of gate xr (ms^-1) */
        double r;     /* K time-independant inactivation */
        
        /* Slowly Activating Potassium Current */
        double iks;   /* Slowly activating K current (uA/uF) */
        double gks;   /* Channel conductance of slowly activating K current (mS/uF) */
        double eks;   /* Reversal potential of slowly activating K current (mV) */
        double xs;    /* Slowly activating potassium current activation gate*/
        double xsss;  /* Steady-state value of activation gate xs */
        double tauxs; /* Time constant of gate xs (ms^-1) */
		const double prnak = 0.01833;  /* Na/K Permiability Ratio */
        
        /* Time-Independent Potassium Current */
		/*Partly modified from Matsuoka, et al, Jap J Physiol,2003:53:105-123*/
        double iki;    /* Time-independant K current (uA/uF) */
        double gki;    /* Channel conductance of time independant K current (mS/uF) */
        double eki;    /* Reversal potential of time independant K current (mV) */
        double kin;    /* K inactivation */
		double iku ;	/*Attaching rate constant of Magnesium to iki*/
		double ikl ;    /*Detaching rate constant of Magnesium to iki*/
		double ikay ;   /*Attaching rate constant of spermine to iki*/
		double ikby ;   /*Detaching rate constant of spermine to iki*/
		double tauiky ; /*Time constant of spermine attachment*/
		double ikyss ;  /*Steady state of spermine attachment*/
		double iky ;    /*Spermine attachment*/
		double foiki ;  /*Fraction of channel free from attachment of Magnesium*/
		double fbiki ;  /*Fraction of channel with attachment of Magnesium*/
        
        /* Transient Outward Potassium Current */
        double ito;       /* Transient outward current */
        double gito;      /* Maximum conductance of Ito */
        double erevto;    /* Reversal potential of Ito */
        double ato;       /* Ito activation */
        double alphaato;  /* Ito alpha-a rate constant */
        double betaato;   /* Ito beta-a rate constant */
        double tauato;    /* Time constant of a gate */
        double atoss;     /* Steady-state value of a gate */
        double iito;      /* Ito inactivation */
        double alphaiito; /* Ito alpha-i rate constant */
        double betaiito;  /* Ito beta-i rate constant */
        double tauiito;   /* Time constant of i gate */
        double iitoss;    /* Steady-state value of i gate */
                
        /* Sodium-Calcium Exchanger */
        double inaca;               /* NaCa exchanger current (uA/uF) */
		const double kmnancx = 87.5;  /* Na saturation constant for NaCa exchanger */
		const double ksatncx = 0.1;   /* Saturation factor for NaCa exchanger */
		const double kmcancx = 1.38;  /* Ca saturation factor for NaCa exchanger */
        const double gammas = 0.35;  /* Position of energy barrier controlling voltage dependance of inaca */

        /* Sodium-Potassium Pump */
        double inak;    /* NaK pump current (uA/uF) */
        double fnak;    /* Voltage-dependance parameter of inak */
        double sigma;   /* [Na]o dependance factor of fnak */
        const double ibarnak = 1.0933;   /* Max. current through Na-K pump (uA/uF) */
        const double kmnai = 10*(1-0.25*iso);    /* Half-saturation concentration of NaK pump (mM) */
        const double kmko = 1.5;    /* Half-saturation concentration of NaK pump (mM) */
        
        /* Sarcolemmal Ca Pump */
        double ipca;                 /* Sarcolemmal Ca pump current (uA/uF) */
        const double ibarpca = 0.275; /* Max. Ca current through sarcolemmal Ca pump (uA/uF) */
        const double kmpca = 0.0005; /* Half-saturation concentration of sarcolemmal Ca pump (mM) */
        
        /* Ca Background Current */
        double icab;  /* Ca background current (uA/uF) */
        double gcab;  /* Max. conductance of Ca background (mS/uF) */
        double ecan;  /* Nernst potential for Ca (mV) */

        /* Na Background Current */
        double inab;  /* Na background current (uA/uF) */
        double gnab;  /* Max. conductance of Na background (mS/uF) */
        double enan;  /* Nernst potential for Na (mV) */
      
        /* Total Ca current */
        double icatot;

        /* Ion Current Functions */     
        void comp_ina ();    /* Calculates Fast Na Current */
        void comp_ical ();   /* Calculates Currents through L-Type Ca Channel */
        void comp_ikr ();    /* Calculates Rapidly Activating K Current */
        void comp_iks ();    /* Calculates Slowly Activating K Current */
        void comp_iki ();    /* Calculates Time-Independant K Current */
		void comp_ikach();   /* Calculates Acetylcholine-sensitive potassium*/
		void comp_ikur ();   /* Calculates Ultra-Rapidly activation K Current*/
        void comp_ito ();    /* Calculates Transient Outward Current */
        void comp_inaca ();  /* Calculates Na-Ca Exchanger Current */
        void comp_inak ();   /* Calculates Na-K Pump Current */
        void comp_ipca ();   /* Calculates Sarcolemmal Ca Pump Current */
        void comp_icab ();   /* Calculates Ca Background Current */
        void comp_inab ();   /* Calculates Na Background Current */
        void comp_it ();     /* Calculates Total Current */

        /* Ion Concentration Functions */
        void conc_nai ();    /* Calculates new myoplasmic Na ion concentration */
        void conc_ki ();     /* Calculates new myoplasmic K ion concentration */
        void conc_nsr ();    /* Calculates new NSR Ca ion concentration */
        void conc_jsr ();    /* Calculates new JSR Ca ion concentration */
        void calc_itr ();    /* Calculates Translocation of Ca from NSR to JSR */
        void conc_cai ();    /* Calculates new myoplasmic Ca ion concentration */
		void TPcai ();

void main ()
{       
        /* Opening of Datafiles */



   str3 =(char *) malloc(32*sizeof(char));
   if (str3 == NULL){
		printf("\n out of memory ");
		exit(1);
	}
	sprintf(str3,"Pitx2_3_LA_%d.txt", bcl);

	ap3 = fopen(str3,"w");
	if (ap3 == NULL)
	{
		printf("\n cannot open  file ");
		exit(1);
	}
	free(str3);



	   strfmax =(char *) malloc(32*sizeof(char));
   if (strfmax == NULL){
		printf("\n out of memory ");
		exit(1);
	}
     sprintf(strfmax,"Pitx2_3_APD_LA_%d.txt", bcl);

	fmaxs = fopen(strfmax,"w");
	if (fmaxs == NULL)
	{
		printf("\n cannot open  file ");
		exit(1);
	}
	free(strfmax);



      
     //   fpara = fopen("fpara.txt", "w");  
    //    fmaxs = fopen("fmaxs.txt", "w");


        /* Cell Geometry */
        vcell = 0.020106;
	    ageo = 0.00005831;
		acap = 0.00011661;
	   // vcell = 1000*pi*a*a*l;     /*   3.801e-5 uL */
       // ageo = 2*pi*a*a+2*pi*a*l;
	  //  acap = ageo*2;             /*   1.534e-4 cm^2 */
        vmyo = vcell*0.68;
       // vmito = vcell*0.26;
       // vsr = vcell*0.06;
       // vnsr = vcell*0.0552;
       // vjsr = vcell*0.0048;

		//Vc=vcell*0.815;
		Vc=vcell;
		Vsr=Vc*0.066;
		Vss=Vc*0.0033;
        
        /* Time Loop Conditions */
        t = 0;           /* Time (ms) */
        udt =0.005;// 0.01;     /* Time step (ms) */
        steps = (S2 + bcl*beats)/udt; /* Number of ms */
        st =-80; // -200;        /* Stimulus */ changed by 
        tstim = 10;       /* Time to begin stimulus */
        stimtime = 10;   /* Initial Condition for Stimulus */
        v = -81.2;       /* Initial Voltage (mv) */

        /* Beginning Ion Concentrations */
        nai = 11.2;       /* Initial Intracellular Na (mM) */
        nao = 140;      /* Initial Extracellular Na (mM) */
        ki = 139;       /* Initial Intracellular K (mM) */
        ko = 4.5;       /* Initial Extracellular K (mM) */
        cai = 0.000102;  /* Initial Intracellular Ca (mM) */
        cao = 1.8;      /* Initial Extracellular Ca (mM) */
        CaSS= 0.000102;  //TP06
        CaSR= 1.3;  //TP06
	   /* Initial Gate Conditions */
        m = 0.00291;
        h = 0.965;
        j = 0.978;
        d = 0.000137;
        f = 0.999837;
		xs = 0.0187;
        xr = 0.0000329;
		ato = 0.0304;
		iito = 0.999;
		uakur = 0.00496;
		uikur = 0.999;
		uskur = 0.999;//BAI 
		f1=1.0;//BAI
		sOO=0.0;//BAI
		sRR=1.0;//BAI
		Flecainide_DO=0.0;
		
		
		fca = 0.775;
		ireljsrol=0;

        
        /* Initial Conditions */
        jsr = 1.49;
        nsr = 1.49;
		trpn = 0.0118;
        cmdn = 0.00205;
        csqn = 6.51;
        boolien = 1;
        dt = udt;
        utsc = 50;
		urel = 0.00;
		vrel = 1.00;
		wrel = 0.999;
		yach = 2.54e-2;
        iky = 0.6;
		i=-1;

        /* Beginning of Time Loop */
        for (increment = 0; increment < steps; increment++)
        {               
               /* List of functions called for each timestep, currents commented out are only used when modeling pathological conditions */
                comp_ina ();
                comp_ical ();
                comp_ikr ();
                comp_iks ();
                comp_iki ();
				comp_ikach ();
                comp_ikur ();
				comp_ito (); 
                comp_inaca ();
                comp_inak ();
                comp_ipca ();
                comp_icab ();
                comp_inab ();
                comp_it ();
                
                conc_nai ();
                conc_ki ();
                TPcai ();
				// conc_nsr ();
               // conc_jsr ();
               // calc_itr ();
               //  conc_cai ();
                
                utsc = 0;
                dt = udt;
		
	prttofile();      

        vnew = v-it*udt;
        dvdtnew = (vnew-v)/udt;

                if (vnew>vmax[i])
                        vmax[i] = vnew;
                if (cai>caimax[i])
                        caimax[i] = cai;
                if (dvdtnew>dvdtmax[i])
                        {dvdtmax[i] = dvdtnew;
                        toneapd[i] = t;}
                if (vnew>=(vmax[i]-0.9*(vmax[i]-rmbp[i])))
                        ttwoapd[i] = t;
                 if (vnew>=(vmax[i]-0.98*(vmax[i]-rmbp[i])))
                        trep[i] = t;


        v = vnew;
        utsc = utsc+1;      
        t = t+udt;     
        }

        printf ("Main loop passed...\n"); 

        if (beats > 0) {       
          apd[beats] = ttwoapd[beats]-toneapd[beats];
          di[i] = toneapd[beats]-ttwoapd[beats-1]; 
          printf("DI = %g, APD = %g\n", di[beats], apd[beats]);
        }

    //    fprintf(fpara,"%.3f\n", t);

        for(i=1;i<beats;i++)
        {apd[i] = ttwoapd[i]-toneapd[i];
         di[i] = toneapd[i]-ttwoapd[i-1]; 
         fprintf(fmaxs, "%i\t%g\t%g\t%g\t%g\n", i, di[i], apd[i], toneapd[i], trep[i]);} 


	fclose(ap3);
	//fclose(fpara);
	fclose(fmaxs);

}

/********************************************************/

/* Functions that describe the currents begin here */


void comp_ina ()
{
        gna = 7.8*(1-AF*0.1)*(1.0/(1.0+Flecainide/84));//PITX2
        ena = ((R*temp)/frdy)*log(nao/nai);

        am = 0.32*(v+47.13)/(1-exp(-0.1*(v+47.13)));
        bm = 0.08*exp(-v/11);
        
        if (v < -40) 
        {ah = 0.135*exp((80+v)/-6.8);  
        bh = 3.56*exp(0.079*v)+310000*exp(0.35*v);
        aj = (-127140*exp(0.2444*v)-0.00003474*exp(-0.04391*v))*((v+37.78)/(1+exp(0.311*(v+79.23))));
        bj = (0.1212*exp(-0.01052*v))/(1+exp(-0.1378*(v+40.14)));}

        else  
        {ah = 0.0;
        bh = 1/(0.13*(1+exp((v+10.66)/-11.1)));
        aj = 0.0;
        bj = (0.3*exp(-0.0000002535*v))/(1+exp(-0.1*(v+32)));}

        h = ah/(ah+bh)-((ah/(ah+bh))-h)*exp(-dt/(1/(ah+bh)));
        j = aj/(aj+bj)-((aj/(aj+bj))-j)*exp(-dt/(1/(aj+bj)));
        m = am/(am+bm)-((am/(am+bm))-m)*exp(-dt/(1/(am+bm)));
        
        ina = gna*m*m*m*h*j*(v-ena);
}

void comp_ical ()
	
{

      dss = 1./(1.+exp((-8-v-3*iso)/7.5));
	  taud = 1.4/(1.+exp((-35-v-3*iso)/13))+0.25+1.4/(1.+exp((v+5+3*iso)/5))+1./(1.+exp((50-v-3*iso)/20));
   	   
	  fss = 1./(1.+exp((v+20+3*iso)/7));
	  tauf= 1102.5*exp(-(v+27+3*iso)*(v+27+3*iso)/225)+200./(1+exp((13-v-3*iso)/10.))+(180./(1+exp((v+30+3*iso)/10)))+20;
	     
	  fss1= 0.67/(1.+exp((v+35+3*iso)/7))+0.33;
	  tauf1=600*exp(-(v+25+3*iso)*(v+25+3*iso)/170)+31/(1.+exp((25-v-3*iso)/10))+16/(1.+exp((v+30+3*iso)/10));
	   
	  fcass = 0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
	  taufca= (80./(1+(CaSS/0.05)*(CaSS/0.05))+2.);
	   d = dss-(dss-d)*exp(-dt/taud);
       f = fss-(fss-f)*exp(-dt/tauf);
	   f1=fss1-(fss1-f1)*exp(-dt/tauf1);
	   fca = fcass-(fcass-fca)*exp(-dt/taufca);
        
		 ibarca=0.00003980*(1+0.5*iso)*(1-0.5*AF)*(1-PITX2*0.5);
		 ilca=ibarca*d*f*f1*fca*4*(v-15)*(frdy*frdy/(R*temp))*(0.25*exp(2*(v-15)*frdy/(R*temp))*CaSS-cao)/(exp(2*(v-15)*frdy/(R*temp))-1.);
	  
	   /*
	   dss = 1/(1+exp(-(v+10)/8));
       taud = (1-exp((v+10)/-6.24))/(0.035*(v+10)*(1+exp((v+10)/-6.24)));

        fss = 1/(1+exp((v+28)/6.9));
        tauf = 9/(0.0197*exp(-pow((0.0337*(v+10)),2))+0.02);

		fcass = 1/(1+CaSS/0.00035);
		taufca = 2;
		ibarca = gcalbar*(v-65);    
		 d = dss-(dss-d)*exp(-dt/taud);
        f = fss-(fss-f)*exp(-dt/tauf);
	    fca = fcass-(fcass-fca)*exp(-dt/taufca);
        
        ilca = d*f*fca*ibarca;
		 */        
        ilcatot = ilca;
		
		
}

void comp_ikr ()
{
        gkr = 0.0294*(1+LA*0.6)*(1+RA*0.0)*sqrt(ko/5.4)*(1.0/(1.0+Flecainide/1.5));
       
		
		ekr = ((R*temp)/frdy)*log(ko/ki);

        xrss = 1/(1+exp(-(v+14.1)/6.5));
        tauxr = 1/(0.0003*(v+14.1)/(1-exp(-(v+14.1)/5))+0.000073898*(v-3.3328)/(exp((v-3.3328)/5.1237)-1));
                
        xr = xrss-(xrss-xr)*exp(-dt/tauxr);
        
        r = 1/(1+exp((v+15)/22.4));
                
        ikr = gkr*xr*r*(v-ekr);
}

void comp_iks ()
{
        gks = 0.129*(1+2.0*iso+1*AF);
        eks = ((R*temp)/frdy)*log(ko/ki);
		tauxs = 0.5/(0.00004*(v-19.9+40*iso)/(1-exp(-(v-19.9+40*iso)/17))+0.000035*(v-19.9+40*iso)/(exp((v-19.9+40*iso)/9)-1));
		xsss = 1/pow((1+exp(-(v-19.9+40*iso)/12.7)),0.5);
		xs = xsss-(xsss-xs)*exp(-dt/tauxs);

		iks = gks*xs*xs*(v-eks);
}

void comp_iki ()
{
        gki = (1+1*AF)*0.09*pow(ko/5.4,0.4);
        eki = ((R*temp)/frdy)*log(ko/ki);

        kin = 1/(1+exp(0.07*(v+80)));
        
        iki = gki*kin*(v-eki);
/*

  // modified from Matsuoka, et al Jap J Physiol 2003:53:105-123
	iku = 0.75*exp(0.035*(v-eki-10))/(1+exp(0.015*(v-eki-140)));
	ikl = 3*exp(-0.048*(v-eki-10))*(1+exp(0.064*(v-eki-38)))/(1+exp(0.03*(v-eki-70)));
	ikay =1/(8000*exp((v-eki-97)/8.5)+7*exp((v-eki-97)/300));
	ikby =1/(0.00014*exp(-(v-eki-97)/9.1)+0.2*exp(-(v-eki-97)/500));
	tauiky = 1/(ikay+ikby);
	ikyss = ikay/(ikay+ikby);
	iky = ikyss - (ikyss-iky)*exp(-dt/tauiky);
	foiki = ikl/(iku+ikl);
	fbiki = iku/(iku+ikl);

	
	iki = gki*(pow(foiki,4)+8*pow(foiki,3)*fbiki/3+2*foiki*foiki*fbiki*fbiki)*iky*(v-eki); 
*/
  }

void comp_ikach ()
{
		gkach = 0.135;
		ekach = ((R*temp)/frdy)*log(ko/ki);
		alphayach= 1.232e-2/(1+0.0042/ach)+0.0002475;
		betayach = 0.01*exp(0.0133*(v+40));
		tauyach = 1/(alphayach+betayach);
		yachss = alphayach/(alphayach+betayach);

		yach = yachss-(yachss-yach)*exp(-dt/tauyach);
		ikach = gkach*yach*(v-ekach)/(1+exp((v+20)/20));
}

void comp_ikur ()
{
		gkur = (1+2.0*iso)*(1-0.5*AF)*(0.005+0.05/(1+exp(-(v-15)/13)));
        ekur = ((R*temp)/frdy)*log(ko/ki);
		
		alphauakur = 0.65/(exp(-(v+10)/8.5)+exp(-(v-30)/59.0));
		betauakur = 0.65/(2.5+exp((v+82)/17.0));
		tauuakur = 1/(3*(alphauakur+betauakur));
		uakurss = 1/(1+exp(-(v+30.3)/9.6));
		
		alphauikur = 1/(21+exp(-(v-185)/28));
		betauikur = exp((v-158)/16);
		tauuikur = 1/(3*(alphauikur+betauikur));
		uikurss = 1/(1+exp((v-99.45)/27.48));
		
		//tauuikur = 800*(2-(v/40.0));            //BAI
		//uikurss = 1+exp((v-35.0)/20.0);         //BAI

		//tauuskur = 5800/(1+exp(-(v+80.0)/11.0));//BAI
		//uskurss = 1+exp((v+5.0)/5.0);           //BAI

		uakur = uakurss-(uakurss-uakur)*exp(-dt/tauuakur);//BAI
		uikur = uikurss-(uikurss-uikur)*exp(-dt/tauuikur);//BAI
		
		//uskur = uskurss-(uskurss-uskur)*exp(-dt/tauuskur);//BAI
		//ikur = gkur*uakur*uakur*uakur*uikur*uskur*(v-ekur);//BAI

		ikur = gkur*uakur*uakur*uakur*uikur*(v-ekur);
		
		//uakurss=1./(1.+exp((v+6.)/-8.6)); //Atrial_1
        //tauuakur=9./(1.+exp((v+5.)/12.0))+0.5;//Atrial_1
	   // uikurss=1./(1.+exp((v+7.5)/10.));//Atrial_1
	   // tauuikur=(590./(1.+exp((v+60.)/10.0))+3050.);//Atrial_1
		//uakur = uakurss-(uakurss-uakur)*exp(-dt/tauuakur);//Atrial_1
		//uikur = uikurss-(uikurss-uikur)*exp(-dt/tauuikur);//Atrial_1
		//gkur =0.045; 
		//ikur = gkur*uakur*uikur*(v-ekur);
}

void comp_ito ()
{
        gito = 0.1652*(1-0.7*AF);
        erevto = ((R*temp)/frdy)*log(ko/ki);
        
        alphaato = 0.65/(exp(-(v+10)/8.5)+exp(-(v-30)/59));
        betaato = 0.65/(2.5+exp((v+82)/17));
        tauato = 1/(3*(alphaato+betaato));
        atoss = 1/(1+exp(-(v+20.47)/17.54));
        ato = atoss-(atoss-ato)*exp(-dt/tauato);

        alphaiito = 1/(18.53+exp((v+113.7)/10.95));
        betaiito = 1/(35.56+exp(-(v+1.26)/7.44));
        tauiito = 1/(3*(alphaiito+betaiito));
        iitoss = 1/(1+exp((v+43.1)/5.3));
        iito = iitoss-(iitoss-iito)*exp(-dt/tauiito);
        
        ito = gito*ato*ato*ato*iito*(v-erevto);
}

void comp_inaca ()
{
	inaca = (1+0.4*AF)*1750*(exp(gammas*frdy*v/(R*temp))*nai*nai*nai*cao-exp((gammas-1)*frdy*v/(R*temp))*nao*nao*nao*cai)/((pow(kmnancx,3)+pow(nao,3))*(kmcancx+cao)*(1+ksatncx*exp((gammas-1)*frdy*v/(R*temp))));
}

void comp_inak ()
{
        sigma = (exp(nao/67.3)-1)/7;

        fnak = 1/(1+0.1245*exp((-0.1*v*frdy)/(R*temp))+0.0365*sigma*exp((-v*frdy)/(R*temp)));
		//fnak=(v+150)/(v+200);
        inak = ibarnak*fnak*(1/(1+pow((kmnai/nai),1.5)))*(ko/(ko+kmko));
}

void comp_ipca ()
{
        ipca = (ibarpca*cai)/(kmpca+cai);       
}

void comp_icab ()
{
        gcab = 0.00113;
        ecan = ((R*temp)/frdy)*log(cao/cai);
                
        icab = gcab*(v-ecan);
}

void comp_inab ()
{
        gnab = 0.000674;
        enan = ((R*temp)/frdy)*log(nao/nai);
        
        inab = gnab*(v-enan);
}

/* Total sum of currents is calculated here, if the time is between stimtime = 0 and stimtime = 0.5, a stimulus is applied */
void comp_it ()
{
        naiont = ina+inab+3*inak+3*inaca+1.5e-2;
        kiont = ikr+iks+iki-2*inak+ito+ikur+ikach+1.5e-2;  //changed by Bai
		//kiont = ikr+iks+iki-2*inak+ito+ikur+1.5e-2; 
        caiont = ilca+icab+ipca-2*inaca;
        
        if (t>=tstim && t<(tstim+dt))     
        {stimtime = 0;
        i = i+1;
        if (i < beats-1) 
            tstim = tstim+bcl; 
        else tstim = tstim+S2; 
        printf ("Stimulus %d applied, ", i+1);
        printf ("Time = %f\n", t);
        boolien = 0;
        rmbp[i]=v;
        nair[i] = nai;
        cair[i] = cai;}
        
        if(stimtime>=0 && stimtime<0.5&&i<(beats-1)) 
        {it = st+naiont+kiont+caiont;}
        else
        {it = naiont+kiont+caiont;}

        stimtime = stimtime+dt;
}

/* Functions that calculate intracellular ion concentrations begins here */

void conc_nai ()
{
        dnai = -dt*naiont*acap/(vmyo*zna*frdy);
        nai = dnai + nai;
}

void conc_ki ()
{
        dki = -dt*kiont*acap/(vmyo*zk*frdy);
        ki = dki + ki;
}

void TPcai ()
{


	kCaSR=(maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))))*(1+2*AF+iso*(1-AF));//////////////
	k1=k1_/kCaSR;/////////////
	k2=k2_*kCaSR;/////////////
	dRR=k4*(1-sRR)-k2*CaSS*sRR;/////////////

	sRR+=dt*dRR;///////////
	sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);
	//sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS)-dt*(Flecainide_kon*sOO-Flecainide_koff*Flecainide_DO);//////////////
	dFlecainide_DO=dt*(Flecainide_kon*sOO-Flecainide_koff*Flecainide_DO);
	Flecainide_DO=dFlecainide_DO+dFlecainide_DO;
	
	//Irel=((1+PITX2*0.3)*Vrel*(sOO)*(CaSR-CaSS)+(1+0.25*AF)*Vleak*(sRR)*(CaSR-CaSS))*(1.0/(1.0+Flecainide/55.0));

		if(Flecainide!=0.0)
	{
		Irel=((1+PITX2*0.3)*Vrel*(1.0/(1.0+Flecainide/55.0))*(sOO)*(CaSR-CaSS)+(1+0.25*AF)*Vleak*(sRR)*(1-1.0/(1.0+Flecainide/55.0))*(CaSR-CaSS));
	}
	else
	{
		Irel=((1+PITX2*0.3)*Vrel*(sOO)*(CaSR-CaSS)+(1+0.25*AF)*Vleak*(sRR)*(CaSR-CaSS));
	}
			
	//Irel=Vrel*(sOO+0.2*Flecainide_DO)*(CaSR-CaSS)+(1+0.25*AF)*Vleak*(sRR+0.2*Flecainide_DO)*(CaSR-CaSS);/////////// Atrial
	Ileak=Vleak*sRR*(CaSR-CaSS);
	Irelease=Vrel*sOO*(CaSR-CaSS);
//	Ileak=Vleak*(CaSR-Cai);////////////////
	//Iup=(1+PITX2*2.0)*Vmaxup/(1.+((Kup*Kup)/(cai*cai)));///////////////
	Iup=(1+PITX2*2.0)*Vmaxup/(1.+((Kup*Kup)/(cai*cai)));///////////////
	Ixfer=Vxfer*(CaSS-cai);////////////

	CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);/////////////
	dCaSR=dt*(Iup-Irel);///////////////Atrial
	bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;//////////////
	cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);/////////////
	CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;/////////////
	
	//double Vc，Vss，Vsr；
    //double Vc,Vss,Vsr;

	CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);//////////////
	dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ilca*(1./(2*Vss*frdy))*CAPACITANCE));/////////////
	bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;///////////
	ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);///////////
	CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;//////////////////

	CaBuf=Bufc*cai/(cai+Kbufc);///////////
	dcai=dt*((-(icab+ipca-2*inaca)*(1./(2*Vc*frdy))*CAPACITANCE)-(Iup)*(Vsr/Vc)+Ixfer);///////////////Atrial
	bc=Bufc-CaBuf-dcai-cai+Kbufc;////////////////////
	cc=Kbufc*(CaBuf+dcai+cai);///////////////
	cai=(sqrt(bc*bc+4*cc)-bc)/2;////////////////////
}

/*
void conc_nsr ()
{
       
        kleak = iupbar/nsrbar;
        ileak = kleak*nsr;

        iup = iupbar*cai/(cai+kmup);
		csqn = csqnbar*(jsr/(jsr+kmcsqn));
		
		dnsr = dt*(iup-ileak-itr*vjsr/vnsr);
        nsr = dnsr+nsr;
}

void conc_jsr ()
{

		fn = vjsr*(1e-12)*ireljsrol-(5e-13)*(ilca/2+inaca/5)*acap/frdy; 
//		fn = vjsr*(1e-12)*ireljsrol-(1e-12)*caiont*acap/(2*frdy); 

		tauurel = 8.0;
		urelss = 1/(1+exp(-(fn-3.4175e-13)/13.67e-16));
		tauvrel = 1.91+2.09/(1+exp(-(fn-3.4175e-13)/13.67e-16));
		vrelss = 1-1/(1+exp(-(fn-6.835e-14)/13.67e-16));
		tauwrel = 6.0*(1-exp(-(v-7.9)/5))/((1+0.3*exp(-(v-7.9)/5))*(v-7.9));
		wrelss = 1-1/(1+exp(-(v-40)/17));

		urel = urelss-(urelss-urel)*exp(-dt/tauurel);
		vrel = vrelss-(vrelss-vrel)*exp(-dt/tauvrel);
		wrel = wrelss-(wrelss-wrel)*exp(-dt/tauwrel);
		
        greljsrol = grelbarjsrol*urel*urel*vrel*wrel;
        ireljsrol = greljsrol*(jsr-cai); 

        
		djsr = dt*(itr-0.5*ireljsrol)/(1+csqnbar*kmcsqn/pow((jsr+kmcsqn),2)); //LAI

        jsr = djsr+jsr;
}

void calc_itr ()
{
        itr = (nsr-jsr)/tautr;
}

void conc_cai ()
{
        trpn = trpnbar*(cai/(cai+kmtrpn));
        cmdn = cmdnbar*(cai/(cai+kmcmdn));
		
		b1cai = -caiont*acap/(2*frdy*vmyo)+(vnsr*(ileak-iup)+0.5*ireljsrol*vjsr)/vmyo; //LAI
		b2cai = 1+trpnbar*kmtrpn/pow((cai+kmtrpn),2)+cmdn*kmcmdn/pow((cai+kmcmdn),2);
		dcai = dt*b1cai/b2cai;


        cai = dcai+cai;
}
*/

/* Values are printed to a file called ap. The voltage and currents can be plotted versus time using graphing software. */

void prttofile ()
{
//if((increment%100 == 0)&&t>(beats-2)*bcl)
	if((increment%50 == 0)&&t>(beats-2)*bcl)
{
	fprintf (ap3, "%.3f\t%g\t%g\t%g\n", t-(beats-2)*bcl, v, cai,CaSR);

	
  // fprintf (ap3, "%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t-(beats-2)*bcl, v, ina, inab, inak, inaca, ilca, icab, ipca, ikr, iks, iki, ito, ikur,  ikach, ki, cai, nai, CaSS, CaSR, Irel, Iup, Ixfer);

}	
}
