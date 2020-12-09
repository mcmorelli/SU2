/*!
 * \file postprocessing_structure.hpp
 * \brief Headers of the post processing structure.
 * \author T. Albring, Beckett Y. Zhou
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once


#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

//#include "fluid_model.hpp"
//#include "numerics_structure.hpp"
//#include "variable_structure.hpp"

//#include "../../Common/include/geometry_structure.hpp"
//#include "../../Common/include/config_structure.hpp"
//#include "../../Common/include/matrix_structure.hpp"
//#include "../../Common/include/vector_structure.hpp"
//#include "../../Common/include/linear_solvers_structure.hpp"
//#include "../../Common/include/grid_movement_structure.hpp"
//#include "../../SU2_CFD/include/solver_structure.hpp"

//#include "numerics_machine_learning.hpp"

//MCM
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/CConfig.hpp"

#include "../../Common/include/grid_movement/CGridMovement.hpp"
#include "../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../Common/include/grid_movement/CVolumetricMovement.hpp"

#include "../include/linear_algebra/CSysMatrix.hpp"
#include "../include/linear_algebra/CSysVector.hpp"
#include "../include/linear_algebra/CSysSolve.hpp"
#include "../include/solvers/CSolver.hpp"

//#include "../../Common/include/math_op_structure.hpp"
//#include "../../Common/include/toolboxes/printing_toolbox.hpp"

using namespace std;


class FWHSolver {
private:
  //SU2_MPI::Comm SU2_Communicator; /*!< \brief MPI communicator of SU2.*/
  int rank, size;                 /*!< \brief MPI rank and size.*/

public:
  su2double  CFD_PressureFluctuation;
  su2double  CAA_PressureFluctuation;
  su2double U1, U2, U3, M, a_inf, AOA, beta_sq, FreeStreamDensity, FreeStreamPressure;
  su2double ***dJdU, ***dJdX;
  su2double *surface_geo,  *n, *RetSurf;
  su2double  *Q;
  su2double *F, *RHO;
  complex <su2double>  *pp_fft;
  unsigned long nObserver, nPanel, nSample, nDim, SamplingFreq,nFace, nqSample;
  unsigned long *PointID;
  su2double **Observer_Locations;
  su2double SPL;
  bool TimeDomain3D;
  su2double T_a, Freq_a, Amp_a;
  unsigned long **Conn;
  unsigned short *nNodeOnFace2;
  unsigned long *LocNode , *FaceNodeList, iSample, nFaceGlobal,nLocNode;
  su2double *AREA,*StartTime,*EndTime,*min_dtt,*max_dtt, *dt;
  su2double rho,rho_ux,rho_uy,rho_uz,rho_E,TKE;
  unsigned long nEdgeNodes, nGlobalEdgeNodes;
  su2double **edgeNodeNorms, **edgeNodeNormsGlobal;
  su2double **edgeNodeCoords, **edgeNodeCoordsGlobal;
  unsigned long *edgeNodeID, *edgeNodeIDGlobal;
  unsigned long *edgeNodeLocalID, *edgeNodeLocalIDGlobal;
  unsigned long *FaceNodeList2;
  su2double *RadVec;
  unsigned long nHaloFaces;
  unsigned long nSample2;
  unsigned short iZone, nZone;
  unsigned long **HaloLayer, **HaloLayer2;
  su2double pp_t, T1,T2,T3,T4,T5, *pp_out;
  su2double *x00,*y00,*z00;
  unsigned long Max_nPanel, Tot_nPanel, nFaceTot, nHaloFacesTot, nFaceMax;
  unsigned long *nConnSizeMaster, nConnSizeMax, nConnSizeTot,*Conn_single_master, *nFaceSizeMaster, *CellTypeMaster, *nPanelMaster;
  su2double *surfaceCoordsMaster, *surfaceCoords, *SingleMaster_pp, *Single_pp;
  su2double *Single_dJdU, *Single_dJdX, *SingleMaster_dJdU, *SingleMaster_dJdX;
	/*!
	 * \brief Constructor of the  class.
	 */
	FWHSolver(CConfig *config, CGeometry *geometry);

	/*!
	 * \brief Destructor of the class.
	 */
	~FWHSolver(void);



        void SetAeroacoustic_Analysis(CSolver *solver, CConfig *config, CGeometry *geometry, ofstream &CFD_pressure_file);
        void SetCFD_PressureFluctuation(CSolver *solver, CConfig *config, CGeometry *geometry, unsigned long iObserver);
        void F1A_SourceTimeDominant ( CConfig* config, CGeometry *geometry);
	void Initialize( CConfig *config, CGeometry *geometry);
	void F1A_Formulation ( CConfig *config,unsigned long iObserver, unsigned long iPanel,  unsigned long iSample, unsigned long iLocSample, unsigned long i);
	void ComputeVelocities ( CConfig *config,  unsigned long iSample, unsigned long iLocSample);
	void ComputeNorms( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long iObserver);
	void Read_surface_csv_files( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long Pre_Mot );
	void PrescribedMotion (CConfig *config, unsigned long iSample, unsigned long iLocSample);
	void ComputeMinMaxInc_Time( CConfig *config, CGeometry *geometry );
	void ComputeObserverTime( CConfig *config, unsigned long iObserver, unsigned long iPanel,unsigned long iSample, unsigned long iLocSample, unsigned long i);
	void Connectivity( CConfig *config, CGeometry *geometry);
	void Extract_NoiseSources( CConfig* config, CGeometry *geometry);
 	void Extract_NoiseSources2(CConfig* config, CGeometry *geometry);
        void Compute_FarfieldNoise(CSolver *solver, CConfig* config, CGeometry *geometry);
        void Compute_TimeDomainPanelSignal(CConfig* config);
        void Compute_GreensFunction2D (CConfig* config);
        void Compute_GreensFunction3D (CConfig* config);
        void Window_SourceTerms ();
        void FFT_SourceTermsR2 ();
	void FFT_AcousticPressureSignal(su2double pp_mean, su2double *pp_TimeDomain, unsigned long res_factor);
        void Integrate_Sources (CConfig* config);
        void iFFT_SignalR2 ();
        void Compute_SPL ();
        void Write_Sensitivities(CSolver *solver, CConfig *config, CGeometry *geometry);
        void Paraview_Output(CConfig *config, unsigned long iObserver, unsigned long iSample, unsigned long SURFACE_TYPE, unsigned long DERIVATIVES, unsigned long FIA_TERMS);
	void CombineZones( CConfig *config);
	su2double GetCFD_PressureFluctuation();
//        void iFFT_Signal (CSolver *solver, CConfig* config, CGeometry *geometry);
//        void FFT_SourceTerms (CSolver *solver, CConfig* config, CGeometry *geometry);
//        void SetCAA_PressureFluctuation(CSolver *solver, CConfig* config, CGeometry *geometry, ofstream &SRC_p_file, ofstream &SRC_ux_file, ofstream &SRC_uy_file, ofstream &SRC_rho_file, ofstream &time_file);
//        su2double GetCAA_PressureFluctuation();

};


class SNG {
  public:
  unsigned long  nDim, nSNGPts, NF, NT, N_Tij_Out;
  long Type_JBBN;
  su2double f_min, f_max, dt, U1, U2, U3, a_inf, TKE_ReDimFac, omega_ReDimFac;
  unsigned long *T_ij_OutputIdx;
  su2double **NoiseSourceZone;
  su2double **SNG_Coords, *TKE, *omega, *SNG_CellVol;
  su2double **k_n, *Psi_n, **sigma_n;
  su2double ***u_turb, ***T_tilda, **T_tilda_mean;
  bool GenNewRand;
  su2double J_BBN;
  su2double **dJBBN_dU;
  unsigned long *PointID;
    /*!
     * \brief Constructor of the  class.
     */
    SNG(CConfig *config, CGeometry *geometry);

    /*!
     * \brief Destructor of the class.
     */
    ~SNG(void);


    void SetSNG_Analysis(CSolver *solver, CConfig *config, CGeometry *geometry );
    void Perform_SNG_Analysis();
    void Extract_RANS(CSolver *solver, CConfig* config, CGeometry *geometry);
    void Write_ExtractedRANS();
    void SetRandomFourierModes();
    void Compute_TurbVelocity();
    void Compute_BroadBandNoiseSource();
    void Compute_BBN_ObjFunc();
    void Write_BroadBandNoiseSource();
    void Write_SNGSensitivities();

};

class F1A {
public:
    unsigned short nDim;
    unsigned short iZone, nZone;

    unsigned long nObserver, nPanel, nSample, SamplingFreq, nqSample;
    unsigned long nSurfaceNodes;

    unsigned long *PointID;
    unsigned long *globalIndexContainer;

    su2double FreeStreamPressure;
    su2double FreeStreamDensity;
    su2double U1, U2, U3, a_inf;
    su2double pp_t;
    su2double T1, T2, T3, T4, T5;
    su2double SPL;

    su2double *surface_geo;
    su2double *GridVel;
    su2double *Momentum;
    su2double *Normal;
    su2double *UnitaryNormal;
    su2double *Area;
    su2double *Q;
    su2double *F;
    su2double *RHO;
    su2double *StartTime;
    su2double *EndTime;
    su2double *dt;
    su2double *RadVec;
    su2double *Observer_Locations;
    su2double *pp_out;

    complex <su2double>  *pp_fft;

    F1A(CConfig *config, CGeometry *geometry);
    ~F1A(void);

    //new implementation
    void Initialize( CConfig *config, CGeometry *geometry);
    static void UpdateDualGrid(CGeometry *geometry, CConfig *config);
    void ComputeNormal( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long iObserver);
    void SetSurfaceGeom(CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long iObserver);
    void Read_TECPLOT_ASCII( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample );

    //legacy implementation
    void ComputeMinMaxInc_Time( CConfig *config, CGeometry *geometry );
    void ComputeObserverTime( CConfig *config, unsigned long iObserver, unsigned long iPanel,unsigned long iSample, unsigned long iLocSample, unsigned long i);
    void F1A_SourceTimeDominant ( CConfig* config, CGeometry *geometry);
    void F1A_Formulation ( CConfig *config,unsigned long iObserver, unsigned long iPanel,  unsigned long iSample, unsigned long iLocSample, unsigned long i);
    void ComputeVelocities ( CConfig *config,  unsigned long iSample, unsigned long iLocSample);
    void FFT_AcousticPressureSignal(su2double pp_mean, su2double *pp_TimeDomain, unsigned long res_factor);

private:
    int rank, size;

};

#include "../include/postprocessing_structure.inl"
