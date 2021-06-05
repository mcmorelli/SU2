/*!
 * \file postprocessing_structure.cpp
 * \brief Source file of the post-processing structure.
 * \author R. Omur Icke
 * \version 6.2.0 
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

#include "../include/postprocessing_structure.hpp"

#include <fstream>
#include <iostream>
//#include "../include/gsl_sf_bessel.h"
//#include <cmath>
#include <complex>
//#include "../include/gsl_errno.h"
//#include "../include/gsl_fft_real.h"
//#include "../include/gsl_fft_halfcomplex.h"
#include <valarray>
//#include "fft_r2.cpp"
//#include "bessel_new.c"
#include <string.h>
typedef std::complex<su2double> Complex;
typedef std::valarray<Complex> CArray;

//#include <vector>
#include "../../Common/include/math_op_structure.hpp"
#include "../../Common/include/toolboxes/printing_toolbox.hpp"
using namespace PrintingToolbox;

FWHSolver::FWHSolver(CConfig *config,CGeometry *geometry) {

    	unsigned long  i, nMarker,iMarker,panelCount, end_iter, start_iter, iVertex,iPoint,  iSample,iPanel,iObserver, iDim,jDim;
   	su2double FreeStreamTemperature;
    	su2double pi=3.141592653589793;
    	su2double *Coord,   *Normal;
    	su2double  x, y, z, nx, ny, nz,  Area;
    	su2double R = 287.058;
    	nDim = geometry->GetnDim();

    /*--- Store MPI rank and size ---*/

    rank = SU2_MPI::GetRank();
    size = SU2_MPI::GetSize();

	#ifdef HAVE_MPI
  		int rank, nProcessor;
    		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     		MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
	#endif

    	//TimeDomain3D = config->GetTimeDomain3D();
    	TimeDomain3D= true;
    	if (rank==MASTER_NODE){
    		if (TimeDomain3D){
        		if(nDim==2){
        			cout<<endl<<"***************  WARNING!!! Time domain FWH implementation NOT available in 2D!!! Use frequency domain!!!  ***************"<<endl;
          		}else{
         			cout<<endl<<"-------------- Initiating 3D FWH Solver in Time Domain --------------"<<endl;
          		}

      		}else{
        		if(nDim==2){
        			cout<<endl<<"-------------- Initiating 2D FWH Solver in Frequency Domain --------------"<<endl;
          		}else{
         			cout<<endl<<"-------------- Initiating 3D FWH Solver in Frequency Domain --------------"<<endl;
          		}
      		}
          	cout<<endl<<"-----------------------------------------VARIABLE SAMPLING FREQUENCY TEST --------------------------------------"<<endl;
	}


    	SPL = 0.0;
    	//end_iter  = config->GetnExtIter();
    	//start_iter  =  config->GetUnst_RestartIter();

    	start_iter = config->GetRestart_Iter();
    	end_iter = config->GetnTime_Iter();

    	SamplingFreq = config->GetWrt_Sol_Freq_DualTime();    //Sampling Freq: defined as the number of dt's between each snapsho (restart file)
    	nSample  =  config->GetIter_Avg_Objective()/SamplingFreq;
    	if (rank==MASTER_NODE) cout<<"Number of sample:"<<nSample<<endl;


    	/* Setting Observer locations -- change to config option later! */
    	string  text_line;
    	ifstream Observer_LocationFile;
    	string filename = "Observer_Locations.dat";
    	Observer_LocationFile.open(filename.data() , ios::in);
    	if (Observer_LocationFile.fail()) {
        	cout << "There is no file!!! " <<  filename.data()  << "."<< endl;
      		exit(EXIT_FAILURE);
    	}
    	getline (Observer_LocationFile, text_line);
    	istringstream point_line(text_line);
     	point_line >> nObserver ;
    	Observer_Locations = new su2double* [nObserver];
    	for(iObserver = 0;  iObserver <nObserver  ;  iObserver++){
       		Observer_Locations[iObserver] = new su2double[nDim];
       		for (iDim=0; iDim < nDim; iDim++){
         		Observer_Locations[ iObserver][iDim]= 0.0;
       		}
    	}

    	iObserver=0;
  	while (getline (Observer_LocationFile, text_line)) {
        	istringstream point_line(text_line);
        	if (nDim==2){
        		point_line >> Observer_Locations[iObserver][0]>> Observer_Locations[iObserver][1];
          	}
        	if (nDim==3){
        		point_line >> Observer_Locations[iObserver][0]>> Observer_Locations[iObserver][1]>> Observer_Locations[iObserver][2];
          	}
        	iObserver++;
    	}

    	M = config->GetMach();
    	beta_sq = 1-M*M;

    	FreeStreamPressure=config->GetPressure_FreeStream();
    	FreeStreamTemperature = config->GetTemperature_FreeStream();
    	FreeStreamDensity = FreeStreamPressure/R/FreeStreamTemperature;

    	a_inf = sqrt(config->GetGamma()*FreeStreamPressure / FreeStreamDensity);
    	AOA = config->GetAoA();

      	U1 = M*a_inf*cos(AOA*pi/180) ;
      	U2 = M*a_inf*sin(AOA*pi/180) ;
      	U3 = 0.0;    //FIX THIS LATER!!!


  	if (rank==MASTER_NODE) cout<<std::setprecision(15)<<"U1:"<<U1<<" U2:"<<U2<<" M:"<<M<<" a_0:"<<a_inf<<" p_0:"<<FreeStreamPressure<<" rho_0:"<<FreeStreamDensity<<" T_0:"<<FreeStreamTemperature<<endl;


	pp_fft          = new complex <su2double> [nSample];
        for (iSample=0; iSample < nSample; iSample++){
                pp_fft[iSample]= 0.0;
        }
       	StartTime= new su2double  [nObserver];
        EndTime  = new su2double  [nObserver];
        dt  = new su2double  [nObserver];
        for(iObserver = 0; iObserver<nObserver ; iObserver++){
               	StartTime[iObserver]= 0.0;
               	EndTime[iObserver]= 9999999999.0;
               	dt[iObserver]= 0.0;
	}

}//End of function






FWHSolver::~FWHSolver(void) {
  unsigned long iSample, iObserver, iPanel;
  /*
  if (Q != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
    if (Q[iPanel] != NULL)  delete [] Q[iPanel];
    delete [] Q;
  }

  if (n != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
    if (n[iPanel] != NULL)  delete [] n[iPanel];
    delete [] n;
  }
  if (surface_geo != NULL) {
    for (iPanel = 0; iPanel < nPanel; iPanel ++)
     if (surface_geo[iPanel] != NULL) delete []surface_geo[iPanel];
    delete [] surface_geo;
  }*/
  
  if (Observer_Locations != NULL) {
    for (iObserver = 0; iObserver < nPanel; iObserver++)
     if (Observer_Locations[iObserver] != NULL) delete []Observer_Locations[iObserver];
    delete []Observer_Locations ;
  }
}


void FWHSolver::SetAeroacoustic_Analysis(CSolver *solver, CConfig *config,CGeometry *geometry,
                              ofstream &CFD_pressure_file ){


	#ifdef HAVE_MPI
  		int rank;
    		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif
	Initialize(config, geometry);
	Connectivity(config, geometry);
	ComputeMinMaxInc_Time(config, geometry);
	F1A_SourceTimeDominant (config, geometry);
	
}



void FWHSolver::Initialize( CConfig *config, CGeometry *geometry){
	
	#ifdef HAVE_MPI
                int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        #endif
	unsigned long iPoint,iMarker,iPanel,iVertex,iDim,nMarker,iObserver;
	unsigned long nPoint = geometry->GetnPoint();
	
        cout<<"Processor-"<<rank<<" contains "<<nPoint<<" points."<<endl;

        FaceNodeList = new unsigned long[nPoint];
        FaceNodeList2  = new unsigned long[nPoint];
        for(iPoint = 0;  iPoint< nPoint; iPoint++){
                FaceNodeList[iPoint] = -5;
                FaceNodeList2[iPoint] = -5;
        }


        unsigned long iFace;
        unsigned long iFaceGlobal=0;
        unsigned long iMarker_nFace;
        unsigned short nNodeOnFace;
        unsigned long checkNode;
        unsigned long iLocNode=0;
        unsigned long CellNode[4]= {0,0,0,0};
        unsigned long Node_0,Node_1,Node_2,Node_3;
        unsigned long iEdgeNode,iGlobalEdgeNode,iHaloFace;
        iEdgeNode=0;
        iHaloFace=0;
        nMarker      = config->GetnMarker_All();

        nPanel=0;
        nFace=0;
        for (iMarker = 0; iMarker < nMarker; iMarker++){
                if (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {
                        nFace=nFace+geometry->GetnElem_Bound(iMarker);
                        nPanel=nPanel+ (geometry->GetnVertex(iMarker));
                }
        }

        LocNode = new unsigned long[nPanel];
        for(iVertex = 0; iVertex < nPanel ; iVertex++){
                LocNode[iVertex] = 0;
        }
	nNodeOnFace2= new unsigned short[nFace];
        for (int i=0; i<nFace;i++){
                nNodeOnFace2[i]=-5;
        }

        for (iMarker = 0; iMarker < nMarker; iMarker++){
                 /* --- Loop over boundary markers to select those on the FWH surface --- */
                if (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {
                        iMarker_nFace=geometry->GetnElem_Bound(iMarker);
                        for(iFace = 0; iFace < iMarker_nFace ; iFace++){
                                nNodeOnFace = geometry->bound[iMarker][iFace]->GetnNodes();
                                checkNode=0;
                                for (int i=0; i<nNodeOnFace; i++){
                                        CellNode[i]=geometry->bound[iMarker][iFace]->GetNode(i);
                                        //if (geometry->node[CellNode[i]]->GetDomain()) checkNode=checkNode+1;

                                        if ( geometry->nodes->GetDomain(CellNode[i]) ){
                                            checkNode=checkNode+1;
                                        }

					//cout<<"First reading CellNode[i]: "<<CellNode[i]<<" checkNode: "<<checkNode<<" GlobalID: "<<geometry->node[CellNode[i]]->GetGlobalIndex()<<" ";
                                }
				//cout<<endl;
                                if (checkNode==nNodeOnFace){//All nodes in the domain
                                        nNodeOnFace2[iFaceGlobal]=nNodeOnFace;
                                        for (int i=0; i<nNodeOnFace; i++){
                                                if (FaceNodeList[CellNode[i]]==-5){
                                                        LocNode[iLocNode]=CellNode[i];//Local node to Global node list
                                                        FaceNodeList[CellNode[i]]=iLocNode+1;// Global to local+1 node list
                                                        iLocNode++;
							//cout<<"CellNode[i]: "<<CellNode[i]<<" FaceNodeList: "<< FaceNodeList[CellNode[i]]<<" ";
                                                }
                                        }
					//cout<<endl;
                                        iFaceGlobal++;
                                }else if( checkNode>0 && checkNode<nNodeOnFace  ){//One or more nodes in the domain
                                        /*Here we evaluate halo layers.*/
                                        for (int i=0; i<nNodeOnFace; i++){

                                                //if (geometry->node[CellNode[i]]->GetDomain()){
                                                if ( geometry->nodes->GetDomain(CellNode[i]) ){
                                                        if (FaceNodeList[CellNode[i]]==-5){
                                                                LocNode[iLocNode]=CellNode[i];//Local node to Global node list
                                                                FaceNodeList[CellNode[i]]=iLocNode+1;// Global to local+1 node list
                                                                iLocNode++;
                                                        }
                                                        if (FaceNodeList2[CellNode[i]]==-5){
                                                                //FaceNodeList2[CellNode[i]]=geometry->node[CellNode[i]]->GetGlobalIndex();
                                                                FaceNodeList2[CellNode[i]]=geometry->nodes->GetGlobalIndex(CellNode[i]);

                                                                iEdgeNode++;
                                                        }
                                                }
                                        }
                                        iHaloFace++;
                                }
                        }//iFace Loop
                }//if heat flux
        }//Marker Loop
        nPanel=iLocNode;
        nEdgeNodes=iEdgeNode;
         nFace=iFaceGlobal;
        nHaloFaces=iHaloFace;
	cout<<"rank-"<<rank<<" contains "<<nPanel<<" panels."<<endl;
        cout<<"rank-"<<rank<<" contains "<<nFace<<" faces."<<endl;
        cout<<"rank-"<<rank<<" contains "<<nHaloFaces<<" halo layers."<<endl;
        cout<<"rank-"<<rank<<" contains "<<nEdgeNodes<<" edge nodes."<<endl;

        #ifdef HAVE_MPI
                HaloLayer = new unsigned long*[nHaloFaces];
                HaloLayer2= new unsigned long*[nHaloFaces];
                for (int i=0; i<nHaloFaces;i++){
                        HaloLayer[i] = new unsigned long[5];
                        HaloLayer2[i]= new unsigned long[9];
                        for(int j = 0; j < 9 ; j++){
                                HaloLayer2[i][j]=-5;
                                if(j<5){
                                        HaloLayer[i][j] =-5;
                                }
                        }
                }

                SU2_MPI::Allreduce(&nEdgeNodes, &nEdgeNodes, 1 , MPI_UNSIGNED_LONG, MPI_MAX , MPI_COMM_WORLD);
                SU2_MPI::Allreduce(&nEdgeNodes, &nGlobalEdgeNodes, 1 , MPI_UNSIGNED_LONG, MPI_SUM , MPI_COMM_WORLD);
                cout<<"rank-"<<rank<<" has "<<nEdgeNodes<<" nEdgeNodes and "<<nGlobalEdgeNodes<<" nGlobalEdgeNodes."<<endl;
                edgeNodeIDGlobal = new unsigned long [nGlobalEdgeNodes];
                edgeNodeLocalIDGlobal= new unsigned long [nGlobalEdgeNodes];
                edgeNodeID = new unsigned long [nEdgeNodes];
                edgeNodeLocalID = new unsigned long [nEdgeNodes];
                edgeNodeCoords = new su2double *[3];
                edgeNodeCoordsGlobal = new su2double *[3];

                for (int j=0; j<3; j++){
                        edgeNodeCoords[j] = new su2double [nEdgeNodes];
                        edgeNodeCoordsGlobal[j]= new su2double [nGlobalEdgeNodes];
                        for (int i=0; i<nEdgeNodes; i++){
                                if (j==0){
                                        edgeNodeID[i]=-5;
                                        edgeNodeLocalID[i]=-5;
                                }
                                edgeNodeCoords[j][i]=0.0;
                        }

                        for (int ii=0; ii<nGlobalEdgeNodes; ii++){
                                if (j==0){
                                        edgeNodeIDGlobal[ii]=-5;
                                        edgeNodeLocalIDGlobal[ii]=-5;
                                }
                                edgeNodeCoordsGlobal[j][ii]=0.0;
                        }
                }

        #endif
	
	Conn = new unsigned long*[nFace];
        for (int i=0; i<nFace;i++){
                if (nNodeOnFace2[i]==4){
                        Conn[i] = new unsigned long[5];
                        Conn[i][0] =4;
                        for(int iDim = 0; iDim < 4 ; iDim++){
                                Conn[i][iDim+1] =0;
                        }
                }else{
                        Conn[i] = new unsigned long[4];
                        Conn[i][0] =3;
                        for(int iDim = 0; iDim < 3 ; iDim++){
                                Conn[i][iDim+1] =0;
                        }
                }
        }

        /*Connectivit Array: Conn
 *
 *                 1st Col         2nd Col     3rd Col     4th Col    5th Col
 *                                 nNode           1st node    2nd node    3rd node   4thNode if exists)
 *
 *                                             End........*/

        PointID = new unsigned long [nPanel];
        for(iPanel = 0;  iPanel< nPanel; iPanel++){
                PointID [iPanel] = 0;
        }
	
	/*Allocate normals, Velocities, Pressure, surface information (radiation distance, vectors, areas and coords) */
        nqSample=20; //nqSample is quasi or local number of samples.
	surface_geo = new su2double [nPanel*nqSample*nDim];
	n = new su2double [nPanel*nqSample*nDim];
        Q = new su2double [nPanel*nqSample*nDim];
        //if (config->GetOutput_FileFormat() == PARAVIEW) RetSurf= new su2double [nPanel*nDim];

        RetSurf= new su2double [nPanel*nDim];

        AREA = new su2double [nPanel];
        F = new su2double [nPanel*nqSample];
        RadVec= new su2double [15];
        for (int i=0; i<15; i++) RadVec[i]=0.0;
        if (config->GetAD_Mode()) RHO= new su2double [nPanel*nqSample];
	for(iPanel = 0;  iPanel< nPanel; iPanel++){
                AREA[iPanel] = 0.0;
                //if (config->GetOutput_FileFormat() == PARAVIEW) for (iDim=0; iDim<nDim; iDim++)	RetSurf[iPanel*nDim+iDim] = 0.0;

                for (iDim=0; iDim<nDim; iDim++)	{
                    RetSurf[iPanel*nDim+iDim] = 0.0;
                }

                for (iSample=0; iSample < nqSample; iSample++){
                        F[iSample*nPanel+iPanel]= 0.0;
                        if (config->GetAD_Mode())RHO[iSample*nPanel+iPanel]=0.0;
			for (iDim=0; iDim < nDim; iDim++){
                                surface_geo[iSample*nPanel*nDim + iPanel*nDim + iDim ]=0.0;
				n[iSample*nPanel*nDim + iPanel*nDim + iDim ]=0.0;
				Q[iSample*nPanel*nDim + iPanel*nDim + iDim ]=0.0;
                        }
                }
        }

	 /*Allocate partial derivatives*/
        if (config->GetAD_Mode()){
                dJdU = new su2double** [nDim+3];
                for(int iDim = 0; iDim < nDim+3 ; iDim++){
                        dJdU[iDim] = new su2double*[nPanel];
                        for(int iPanel = 0;  iPanel< nPanel; iPanel++){
                                dJdU[iDim][iPanel] = new su2double [nSample];
                                for (int iSample=0; iSample <nSample; iSample++){
                                        dJdU[iDim][iPanel][iSample]= 0.0;
                                }
                        }
                }

                dJdX = new su2double** [nDim];
                for(int iDim = 0; iDim < nDim ; iDim++){
                        dJdX[iDim] = new su2double*[nPanel];
                        for(int iPanel = 0;  iPanel< nPanel; iPanel++){
                                dJdX[iDim][iPanel] = new su2double [nSample];
                                for (int iSample=0; iSample <nSample; iSample++){
                                        dJdX[iDim][iPanel][iSample]= 0.0;
                                }
                        }
                }
        }

}//End of subroutine Initialize

void FWHSolver::Connectivity( CConfig *config, CGeometry *geometry){
	#ifdef HAVE_MPI
                int rank, nProcessor;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
        #endif

	
        unsigned long Gnode[4]= {0,0,0,0},GlobID[4]= {0,0,0,0};
 	unsigned long CellNode[4]={0,0,0,0};

        unsigned long nMarker      = config->GetnMarker_All();
        unsigned long iEdgeNode,iGlobalEdgeNode, iPanel, nNodeOnFace, iFace, iMarker;
	unsigned long nNode;
	char buffer[32];

        unsigned long *Glob2Loc = new unsigned long [nPanel];
        for (int i=0; i <nPanel; i++){
                Glob2Loc[i]=0;
        }


        ofstream Connectivity_file;
        if ( config->GetWrt_Csv_Sol() ) { //writing connectivity data for output to other software
		snprintf(buffer, sizeof(char) * 32, "Connectivity_Zone_%i", (iZone));
                Connectivity_file.open(buffer);
        }

	unsigned long checkNode;
        unsigned long iFaceGlobal = 0;
        iEdgeNode=0;
        unsigned long iHaloFace=0;
	 for (iMarker = 0; iMarker < nMarker; iMarker++){
        /* --- Loop over boundary markers to select those on the FWH surface --- */
                if (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {
                        for (iFace = 0; iFace < geometry->GetnElem_Bound(iMarker) ; iFace++){
                                nNodeOnFace = geometry->bound[iMarker][iFace]->GetnNodes();
				checkNode=0;
				for (int i=0; i<nNodeOnFace; i++){
					CellNode[i]=geometry->bound[iMarker][iFace]->GetNode(i);
					//if (geometry->node[CellNode[i]]->GetDomain()) checkNode=checkNode+1;

					if (geometry->nodes->GetDomain(CellNode[i])){
                        checkNode=checkNode+1;
					}
				}
                                if (checkNode==nNodeOnFace){//All nodes in the domain
					if ( config->GetWrt_Csv_Sol() ) { //writing connectivity data for output to other software
						Connectivity_file<<nNodeOnFace<<' ';
					}
					for (int i=0; i<nNodeOnFace; i++){
						Conn[iFaceGlobal][i+1]=FaceNodeList[CellNode[i]];
						if ( config->GetWrt_Csv_Sol() ) { //writing connectivity data for output to other software
							Connectivity_file<<Conn[iFaceGlobal][i+1]<<' ';
						}				
					}
					if ( config->GetWrt_Csv_Sol() ) Connectivity_file<<endl;
                                	iFaceGlobal++;
                                }else if( checkNode>0 && checkNode<nNodeOnFace ){//One or more nodes in the domain
                                        /*Here we evaluate halo layers.*/
					HaloLayer[iHaloFace][0]=nNodeOnFace;
					for (int i=0; i<nNodeOnFace; i++){
						//HaloLayer[iHaloFace][i+1]=geometry->node[CellNode[i]]->GetGlobalIndex();
						HaloLayer[iHaloFace][i+1]=geometry->nodes->GetGlobalIndex(CellNode[i]);
					}
                                        iHaloFace++;
                                }
                        }//iFace Loop

                }//if statement for HeatFLUX
        }//iMarker loop
	iEdgeNode=0;
	for (iPanel=0 ; iPanel<nPanel; iPanel++){
		//PointID[iPanel] = geometry->node[ LocNode[iPanel] ]->GetGlobalIndex();
		PointID[iPanel] = geometry->nodes->GetGlobalIndex(LocNode[iPanel]);

		if (FaceNodeList2[LocNode[iPanel]]!=-5){
			edgeNodeID[iEdgeNode]=PointID[iPanel];//Global ID
                        edgeNodeLocalID[iEdgeNode]=iPanel+1;//Local ID
			iEdgeNode++;
		}	
	}
	#ifdef HAVE_MPI
		SU2_MPI::Allgather(edgeNodeID, nEdgeNodes, MPI_UNSIGNED_LONG, edgeNodeIDGlobal, nEdgeNodes, MPI_UNSIGNED_LONG,MPI_COMM_WORLD);
     		SU2_MPI::Allgather(edgeNodeLocalID, nEdgeNodes, MPI_UNSIGNED_LONG, edgeNodeLocalIDGlobal, nEdgeNodes, MPI_UNSIGNED_LONG,MPI_COMM_WORLD);		
	#endif

	int k;
	/*Connectivity for HaloLayer2*/
	for (iHaloFace=0; iHaloFace<nHaloFaces; iHaloFace++){
        	nNode=HaloLayer[iHaloFace][0];
		HaloLayer2[iHaloFace][0]=nNode;
                if (nNode!=-5){
                	for (int i=0; i<nNode; i++){
                		Gnode[i]=HaloLayer[iHaloFace][i+1];
             		}
			k=0;
			iGlobalEdgeNode=0;
			while (k<nNode){	
                		GlobID[k]=edgeNodeIDGlobal[iGlobalEdgeNode];
                        	if (GlobID[k]==Gnode[k]){
                                        HaloLayer2[iHaloFace][2*k+1]=floor(iGlobalEdgeNode/nEdgeNodes);//Rank ID
                                        HaloLayer2[iHaloFace][2*(k+1)]=edgeNodeLocalIDGlobal[iGlobalEdgeNode];
					k++;
					iGlobalEdgeNode=-1;
				}
				iGlobalEdgeNode++;
                	}
		}
	}


}//End of the subroutine connectivity

void FWHSolver::ComputeObserverTime( CConfig *config, unsigned long iObserver, unsigned long iPanel,unsigned long iSample, unsigned long iLocSample, unsigned long i){

/*This subroutine computes observer time and radiation vector for each panel */
 unsigned long iDim;
 su2double x,y,z;
 su2double xp[3]= {0.0,0.0,0.0};
 su2double xo[3]= {0.0,0.0,0.0};
 su2double U[3] = {U1, U2, U3};
 su2double dtt=0.1;
 su2double eps=1.0;
 su2double dtnew,diff,R,r1,r2,r3,r_mag,Time;

 //unsigned long start_iter  =  config->GetUnst_RestartIter();
 unsigned long start_iter = config->GetRestart_Iter();

 Time= config->GetDelta_UnstTime()*(start_iter+(iSample+iLocSample)*SamplingFreq);
 //x=surface_geo[iPanel][iLocSample][0];
 x=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 0 ];
 //y=surface_geo[iPanel][iLocSample][1];
 y=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 1 ];
 //z=surface_geo[iPanel][iLocSample][2];
 z=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 2 ];


 x=x-U1*Time;
 y=y-U2*Time;
 z=z-U3*Time;

 for (iDim=0; iDim<nDim; iDim++){ 
	xo[iDim]=Observer_Locations[iObserver][iDim]-U[iDim]*Time;
	xp[iDim]=xo[iDim]-U[iDim]*dtt;
 }
 while (eps>0.0000000001){
	R=sqrt((xp[0]-x)*(xp[0]-x)+(xp[1]-y)*(xp[1]-y)+(xp[2]-z)*(xp[2]-z));
	dtnew=R/a_inf;
        diff=dtnew-dtt;
        dtt=dtnew;
        eps=diff/dtnew*10000;
        eps=sqrt(eps*eps);
	for (iDim=0; iDim<nDim; iDim++) xp[iDim]=xo[iDim]-U[iDim]*dtt;
 }
 r1 = xp[0]-x;
 r2 = xp[1]-y;
 r3 = xp[2]-z;
 r_mag = sqrt(r1*r1+r2*r2+r3*r3);
 r1 = r1/r_mag; r2 = r2/r_mag;r3 = r3/r_mag;

 RadVec[5*i+0]=r1;
 RadVec[5*i+1]=r2;
 RadVec[5*i+2]=r3;
 RadVec[5*i+3]=r_mag;
 RadVec[5*i+4]=Time+dtt;//When the observer hear the signal

}//end of the subroutine

void FWHSolver::ComputeMinMaxInc_Time(CConfig *config, CGeometry *geometry){
	unsigned long iLocSample, iObserver, iPanel;

	#ifdef HAVE_MPI
                int rank, nProcessor;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
        #endif


	su2double tau;
	iLocSample=0;
	Read_surface_csv_files( config, geometry, 0, iLocSample, 0);
	Read_surface_csv_files( config, geometry, nSample-1, iLocSample+1, 0);
	
	for (iObserver=0; iObserver<nObserver; iObserver++){
		for (iLocSample=0; iLocSample<2; iLocSample++){
			for (iPanel=0; iPanel<nPanel; iPanel++){
				ComputeObserverTime( config, iObserver, iPanel, iLocSample*(nSample-1), iLocSample,0);
                                tau= RadVec[4];
				if (iLocSample==0){
					if (tau>StartTime[iObserver]){
						StartTime[iObserver]=tau;
					}
				}else{
					if (tau<EndTime[iObserver]){
                                                EndTime[iObserver]=tau;
                                        }
				}
			}
		}
		#ifdef HAVE_MPI
			if (EndTime[iObserver]==0.0){
                        	EndTime[iObserver]=999999999;
                	}
			SU2_MPI::Allreduce(&StartTime[iObserver],&StartTime[iObserver],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			SU2_MPI::Allreduce(&EndTime[iObserver],&EndTime[iObserver],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		#endif	
		dt[iObserver]=(EndTime[iObserver]-StartTime[iObserver])/nSample;
		StartTime[iObserver]=StartTime[iObserver]+5*dt[iObserver];
		EndTime[iObserver]=EndTime[iObserver]-5*dt[iObserver];
		dt[iObserver]=(EndTime[iObserver]-StartTime[iObserver])/nSample;
		if (rank==MASTER_NODE) cout<<"For Observer-"<<iObserver<<", Endtime="<<EndTime[iObserver]<<" and Startime="<<StartTime[iObserver]<<endl;
		if (rank==MASTER_NODE) cout<<"dt="<<dt[iObserver]<<endl;
	}
}//end of subroutine ComputeMinMaxInc_Time


void FWHSolver::Read_surface_csv_files( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long Pre_Mot ){

 //Pre_Mot is switch definition of prescribed motion. if it's 1, it's presribed motion, if 0, it's time dependent.
 #ifdef HAVE_MPI
        int rank, nProcessor;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
 #endif


 unsigned long iPanel,iDim;
 //unsigned long ExtIter = (iSample+iLocSample)*SamplingFreq + config->GetUnst_RestartIter();
 //unsigned long start_iter  =  config->GetUnst_RestartIter();

    unsigned long ExtIter = (iSample+iLocSample)*SamplingFreq + config->GetRestart_Iter();
    unsigned long start_iter = config->GetRestart_Iter();

 su2double dtt=0.1;
 su2double Time= config->GetDelta_UnstTime()*(start_iter+(iSample+iLocSample)*SamplingFreq);
 su2double DTT;
 DTT = config->GetDelta_UnstTime()*SamplingFreq;

 //cout<<"rank-"<<rank<<"ExtIter="<<ExtIter<<" iSample="<<iSample+iLocSample<<" Time="<<Time<<" dt="<<DTT<<endl;
 if (rank==MASTER_NODE){
 	cout<<"ExtIter="<<ExtIter<<" iSample="<<iSample+iLocSample<<" Time="<<Time<<" dt="<<DTT<<endl;
 }

 string data;
 string iddd,xxx,yyy,zzz,ppp,cpp,rho0;
 unsigned long iddd2;
 su2double ppp2,rho1,xxx2,yyy2,zzz2,x,y,z;
 ifstream csvFile;
 ofstream Outputs;
 if ( (iSample==0) && (iLocSample==0) ){
 	if ( config->GetWrt_Csv_Sol() ) {//Writing .csv file for output to other software
 		char buffer[32]; // The filename buffer.
     		snprintf(buffer, sizeof(char) * 32, "output_iZ%i_%05d.csv", (iZone),(ExtIter));		
   		Outputs.open(buffer);
 	}
 }
 char buffer[32]; // The filename buffer.
 if (nZone>1){
	snprintf(buffer, sizeof(char) * 32, "surface_flow_Z_%i_%05d.csv",iZone, ExtIter);
 }else{
 	snprintf(buffer, sizeof(char) * 32, "surface_flow_%05d.csv", ExtIter);
 }
 csvFile.open(buffer);
 //cout<<"rank-"<<rank<<" Reading -> "<<buffer<<endl;
 if (rank==MASTER_NODE){
       	cout<<"Reading -> "<<buffer<<endl;
 	if (csvFile.fail()) {
               	cout << "There is no file!!! " <<  buffer  << "."<< endl;
              	exit(EXIT_FAILURE);
        }
 }
 getline(csvFile,data); //to get rid of first line that is headerline
 csvFile.close();
 csvFile.open(buffer);
 getline(csvFile,data); //to get rid of first line that is headerline
 while(getline(csvFile,data)){
       	istringstream point_data(data);
        getline(point_data,iddd,',');
        getline(point_data,xxx,',');
        getline(point_data,yyy,',');
        getline(point_data,zzz,',');
        getline(point_data,ppp,',');
        if (config->GetAD_Mode()){
       		getline(point_data,cpp,',');
                getline(point_data,rho0,',');
                stringstream geek6(rho0);
                geek6 >> rho1;
    	}
        stringstream geek(iddd);
        geek >> iddd2;
        stringstream geek2(ppp);
        geek2 >> ppp2;
        if ((geometry->GetGlobal_to_Local_Point(iddd2)) > -1) {
		if ( FaceNodeList[geometry->GetGlobal_to_Local_Point(iddd2)]!=-5 ){
			//cout<<"Read-hey-6 iddd2-"<<iddd2<<" Global index of iddd2- "<<geometry->GetGlobal_to_Local_Point(iddd2)<<endl;
        		iPanel = FaceNodeList[geometry->GetGlobal_to_Local_Point(iddd2)]-1;
			//cout<<"Read-hey-7 iPanel= "<<iPanel<<" NodeID- "<<geometry->GetGlobal_to_Local_Point(iddd2)<<endl;
                	//F[iPanel][iLocSample]=ppp2-FreeStreamPressure;
			F[iLocSample*nPanel+iPanel]=ppp2-FreeStreamPressure;
      			if (config->GetAD_Mode()){
                      		//RHO[iPanel][iLocSample]=rho1;
				RHO[iLocSample*nPanel+iPanel]=rho1;
                	}
                	stringstream geek3(xxx);
                	stringstream geek4(yyy);
                	stringstream geek5(zzz);
			geek3 >> x;
			geek4 >> y;
			geek5 >> z;
			if ( (Pre_Mot==0) ||(iLocSample==0 && iSample==0) ){
				/* if the simulation requires sensitivity computation for dJ/dX::*/
           			if (config->GetAD_Mode()){
                			/*Registration for automatic differentiation*/
                        		AD::RegisterInput(x);
                        		AD::RegisterInput(y);
                        		AD::RegisterInput(z);
            			}
				surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 0 ] = x;
				surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 1 ] = y;
				surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 2 ] = z;
			}
		}
   	}
 }
 csvFile.close();

 if (Pre_Mot==1){
	PrescribedMotion (config,iSample, iLocSample);
 }	

 //cout<<"rank-"<<rank<<" finished reading"<<endl;
 /*Assign new coordinates on edge nodes of halo layers*/
 unsigned long iEdgeNode=0;
 for (iPanel=0; iPanel<nPanel; iPanel++){
	if (FaceNodeList2[LocNode[iPanel]]!=-5){
		for (iDim=0 ; iDim<nDim ;iDim++) edgeNodeCoords[iDim][iEdgeNode]=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + iDim ];
		iEdgeNode++;
	}
 }
 
 #ifdef HAVE_MPI
        for (int i=0; i<nDim ;i++){
                SU2_MPI::Allgather(edgeNodeCoords[i], nEdgeNodes, MPI_DOUBLE, edgeNodeCoordsGlobal[i], nEdgeNodes, MPI_DOUBLE,MPI_COMM_WORLD);
        }
 #endif
 
 if ( (iSample==0) && (iLocSample==0) ){
 	for (iPanel = 0; iPanel < nPanel; iPanel++){
		Outputs<<std::setprecision(15)<<iPanel+1<<","<<surface_geo[iLocSample*nPanel*nDim + iPanel*nDim+0]<<","<<surface_geo[iLocSample*nPanel*nDim + iPanel*nDim+1 ]<<","<<surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]<<endl;
	}
 }





 //if (rank==MASTER_NODE) cout<<".csv file allocated."<<endl;
}//end of subroutine (Read_surface_csv_files)



void FWHSolver::ComputeNorms( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long iObserver){

 #ifdef HAVE_MPI
  	int rank, nProcessor;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
 #endif

 su2double Area;
 unsigned long iFace, iDim, iPanel,iHaloFace,nNode, GlobID1,iEdgeNode, iGlobalEdgeNode, k;
 int i,j;
 su2double CG[3]= {0.0,0.0,0.0},ur[3] = {0.0,0.0,0.0}, vr[3] = {0.0,0.0,0.0};
 unsigned long Gnode[4]= {0,0,0,0},Lnode[4]= {0,0,0,0}, GlobID[4]= {0,0,0,0};


 for (iFace = 0; iFace<nFace; iFace++){
 	nNode=Conn[iFace][0];
        for (i=0; i<nNode; i++){
        	Lnode[i]=Conn[iFace][i+1]-1;
        }
        CG[0]=0.0;CG[1]=0.0;CG[2]=0.0;
        for ( i=0; i<nNode; i++){
        	for (iDim=0;iDim<nDim; iDim++){
			/*Center of gravity coordinates*/
			CG[iDim]=CG[iDim]+(surface_geo[iLocSample*nPanel*nDim + Lnode[i]*nDim + iDim ])/nNode;
      		}
    	}
        for (i=0; i<nNode; i++){
        	j=((i+1) % nNode);
                for (iDim=0; iDim<nDim; iDim++){
                	ur[iDim] = (surface_geo[iLocSample*nPanel*nDim + Lnode[j]*nDim + iDim ]-surface_geo[iLocSample*nPanel*nDim + Lnode[i]*nDim + iDim ])*0.5; //CG of first edge - X0
                        vr[iDim] = CG[iDim]-surface_geo[iLocSample*nPanel*nDim + Lnode[i]*nDim + iDim ]; //CG of face - X0
      		}
                iPanel=Lnode[i];
                /*Vectoral multiplication of vr and ur to find normal vectors*/
		n[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]=n[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]+ur[1]*vr[2]-ur[2]*vr[1];
		n[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]=n[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]+ur[2]*vr[0]-ur[0]*vr[2];
		n[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]=n[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]+ur[0]*vr[1]-ur[1]*vr[0];
                if ((n[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]==0) && (n[iLocSample*nPanel*nDim+iPanel*nDim+1]==0) && (n[iLocSample*nPanel*nDim+iPanel*nDim+2]==0)){
                	//cout<<"n eq 0!!!  pro-"<<rank<<" iPanel-"<<iPanel<<" NodeID:"<<LocNode[iPanel]<<" GlobalID:"<<geometry->node[LocNode[iPanel]]->GetGlobalIndex()<<endl;
                	cout<<"n eq 0!!!  pro-"<<rank<<" iPanel-"<<iPanel<<" NodeID:"<<LocNode[iPanel]<<" GlobalID:"<<geometry->nodes->GetGlobalIndex(LocNode[iPanel])<<endl;
       		}
   	}
 }

 /*Loop over Halo Layer*/
 for (iHaloFace=0; iHaloFace<nHaloFaces; iHaloFace++){
 	nNode=HaloLayer[iHaloFace][0];
	//cout<<"rank-"<<rank<<" iHaloFace="<<iHaloFace<<" HaloLayer- "<<HaloLayer[iHaloFace][1]<<" "<<HaloLayer[iHaloFace][2]<<" "<<HaloLayer[iHaloFace][3]<<" "<<endl;
        if (nNode!=-5){
		for (int i=0; i<nNode; i++){
                	Gnode[i]=HaloLayer[iHaloFace][i+1];
      		}
                k=0;
                iGlobalEdgeNode=0;
                while (k<nNode){
                	GlobID[k]=edgeNodeIDGlobal[iGlobalEdgeNode];
            		if (GlobID[k]==Gnode[k]){
				Lnode[k]=iGlobalEdgeNode;				
                        	k++;
				iGlobalEdgeNode=-1;
                   	}
                 	iGlobalEdgeNode++;
          	}
		//cout<<"rank-"<<rank<<" Lnode: "<<Lnode[0]<<" "<<Lnode[1]<<" "<<Lnode[2]<<" "<<GlobID[0]<<" "<<GlobID[1]<<" "<<GlobID[2]<<" "<<endl;	
		CG[0]=0.0;CG[1]=0.0;CG[2]=0.0;
		for (int i=0; i<nNode; i++){
			for (iDim=0;iDim<nDim; iDim++){
                        	CG[iDim]=CG[iDim]+(edgeNodeCoordsGlobal[iDim][Lnode[i]])/nNode;
                 	}
          	}
		//cout<<"rank-"<<rank<<" CG: "<<CG[0]<<" "<<CG[1]<<" "<<CG[2]<<" "<<endl;
                for (int i=0; i<nNode; i++){
                	j=((i+1) % nNode);
                        if ((geometry->GetGlobal_to_Local_Point(Gnode[i])) > -1) {//if it's in the domain
                        	for (iDim=0; iDim<nDim; iDim++){
                                	ur[iDim] = (edgeNodeCoordsGlobal[iDim][Lnode[j]]-edgeNodeCoordsGlobal[iDim][Lnode[i]])*0.5; //CG of first edge - X0
                                        vr[iDim] = CG[iDim]-edgeNodeCoordsGlobal[iDim][Lnode[i]]; //CG of face - X0
                        	}
                                iPanel = FaceNodeList[geometry->GetGlobal_to_Local_Point(Gnode[i])]-1;
				n[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]=n[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]+(ur[1]*vr[2]-ur[2]*vr[1]);
				n[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]=n[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]+(ur[2]*vr[0]-ur[0]*vr[2]);
				n[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]=n[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]+(ur[0]*vr[1]-ur[1]*vr[0]);
                	}
        	}
   	}
 }

 for (iPanel=0; iPanel<nPanel; iPanel++){
	Area  = 0.0; for ( iDim = 0; iDim < nDim; iDim++)   Area += n[iLocSample*nPanel*nDim + iPanel*nDim + iDim ]*n[iLocSample*nPanel*nDim + iPanel*nDim + iDim ];  Area  = sqrt( Area );
	if ((iSample+iLocSample+iObserver)==0){
		AREA[iPanel]=Area;
	}
	for (iDim=0 ; iDim<nDim ;iDim++){
		n[iLocSample*nPanel*nDim + iPanel*nDim + iDim ]= n[iLocSample*nPanel*nDim + iPanel*nDim + iDim ]/Area;
	}
 }


 //if (rank==MASTER_NODE) cout<<"Normals computed."<<endl;

}//end of subroutine ComputeNorms


void FWHSolver::ComputeVelocities ( CConfig *config,  unsigned long iSample, unsigned long iLocSample){

 unsigned long iPanel,iDim;
 su2double rho,rho_ux, rho_uy,rho_uz,rho_E,TKE,ux,uy,uz,p, StaticEnergy;
 su2double p1[3]= {0.0,0.0,0.0};
 su2double p2[3]= {0.0,0.0,0.0};
 su2double m1[3]= {0.0,0.0,0.0};
 su2double m2[3]= {0.0,0.0,0.0};
 su2double *U = new su2double [3];
 U[0]=U1;U[1]=U2; U[2]=U3;
 su2double DTT,Time_p1,Time_p2,Time_m1,Time_m2;
 DTT = config->GetDelta_UnstTime()*SamplingFreq;

 //unsigned long start_iter  =  config->GetUnst_RestartIter();
    unsigned long start_iter = config->GetRestart_Iter();

 Time_p1= config->GetDelta_UnstTime()*(start_iter+(iSample+iLocSample+1)*SamplingFreq);
 Time_p2= config->GetDelta_UnstTime()*(start_iter+(iSample+iLocSample+2)*SamplingFreq);
 Time_m1= config->GetDelta_UnstTime()*(start_iter+(iSample+iLocSample-1)*SamplingFreq);
 Time_m2= config->GetDelta_UnstTime()*(start_iter+(iSample+iLocSample-2)*SamplingFreq);
        /*Find velocities from node coordinates
 *        x_p1: x plus  1 =x(iSample+1)
 *        x_m1: x minus 1 =x(iSample-1)*/
 for (iPanel=0; iPanel<nPanel; iPanel++){
	for (int i=0; i<nDim; i++){
		p1[i]=surface_geo[(iLocSample+1)*nPanel*nDim + iPanel*nDim + i]-U[i]*Time_p1;
		//p2[i]=surface_geo[iPanel][iLocSample+2][i]-U[i]*Time_p2;
		m1[i]=surface_geo[(iLocSample-1)*nPanel*nDim + iPanel*nDim + i]-U[i]*Time_m1;
		//m2[i]=surface_geo[iPanel][iLocSample-2][i]-U[i]*Time_m2;
		/*4th order accurate central difference*/
		//Q[iPanel][iLocSample][i]=(-p2[i]+8*p1[i]-8*m1[i]+m2[i])/(12*DTT);
		/*2nd order accurate central difference*/
		Q[iLocSample*nPanel*nDim + iPanel*nDim + i ]=(p1[i]-m1[i])/(2*DTT);
	}
	if (config->GetAD_Mode()){
		if ( (iSample==0 && iLocSample==1) ){
			rho = RHO[(iLocSample-1)*nPanel+iPanel];
                        rho_ux = RHO[(iLocSample-1)*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]+U1);
			rho_uy = RHO[(iLocSample-1)*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]+U2);
                        if (nDim==3)  rho_uz = RHO[(iLocSample-1)*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]+U3);
			TKE = 0.0;
                        rho_E = rho*(((F[(iLocSample-1)*nPanel+iPanel]+FreeStreamPressure)/(rho*(config->GetGamma()-1)))+0.5*((ux)*(ux)+(uy)*(uy)+(uz)*(uz))+TKE);
			AD::RegisterInput(rho );
                        AD::RegisterInput(rho_ux );
                        AD::RegisterInput(rho_uy );
                        if (nDim==3) AD::RegisterInput(rho_uz );
                        AD::RegisterInput(rho_E );
                        AD::RegisterInput(TKE);
                }
		/*extract CONSERVATIVE flow data from a particular panel on the FWH surface*/
		rho = RHO[iLocSample*nPanel+iPanel];
                rho_ux = RHO[iLocSample*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]+U1);
                rho_uy = RHO[iLocSample*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]+U2);
		if (nDim==3)  rho_uz = RHO[iLocSample*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]+U3);
		TKE = 0.0;
		rho_E = rho*(((F[iLocSample*nPanel+iPanel]+FreeStreamPressure)/(rho*(config->GetGamma()-1)))+0.5*((ux)*(ux)+(uy)*(uy)+(uz)*(uz))+TKE);	
		AD::RegisterInput(rho );
                AD::RegisterInput(rho_ux );
                AD::RegisterInput(rho_uy );
                if (nDim==3) AD::RegisterInput(rho_uz );
                AD::RegisterInput(rho_E );
                AD::RegisterInput(TKE);
                ux = rho_ux/rho;
                uy = rho_uy/rho;
                uz = 0.0;
                if (nDim==3) uz= rho_uz/rho;
		StaticEnergy =  rho_E/rho-0.5*(ux*ux+uy*uy+uz*uz)-TKE;
                p = (config->GetGamma()-1)*rho*StaticEnergy;

                if ((rho!=0) && (rho_ux!=0) && (rho_uy!=0) && (rho_uz!=0) && (p!=0) ) {
                	 Q[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]= ux-U1;
                         Q[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]= ux-U2;
                         Q[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]= ux-U3;
			F[iLocSample*nPanel+iPanel]=p-FreeStreamPressure;
         	}
		if (iSample==nSample-nqSample-2 && iLocSample==nqSample-2){
			rho = RHO[(iLocSample+1)*nPanel+iPanel];
			rho_ux = RHO[(iLocSample+1)*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 0 ]+U1);	
			rho_uy = RHO[(iLocSample+1)*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 1 ]+U2);	
			if (nDim==3)  rho_uz = RHO[(iLocSample+1)*nPanel+iPanel]*(Q[iLocSample*nPanel*nDim + iPanel*nDim + 2 ]+U3);	
			TKE = 0.0;
			rho_E = rho*(((F[(iLocSample+1)*nPanel+iPanel]+FreeStreamPressure)/(rho*(config->GetGamma()-1)))+0.5*((ux)*(ux)+(uy)*(uy)+(uz)*(uz))+TKE);
			AD::RegisterInput(rho );
			AD::RegisterInput(rho_ux );
			AD::RegisterInput(rho_uy );
			if (nDim==3) AD::RegisterInput(rho_uz );
			AD::RegisterInput(rho_E );
			AD::RegisterInput(TKE);
		}
	}		
 }
 delete [] U;

}//End of the subroutine ComputeVelocities

void FWHSolver::PrescribedMotion (CConfig *config,  unsigned long iSample, unsigned long iLocSample){

 su2double w[3]= {0.0,0.0,0.0};
 su2double u[3]= {0.0,0.0,0.0};
 su2double x[3]= {0.0,0.0,0.0};
 su2double alpha[3]= {0.0,0.0,0.0};
 //su2double Time= config->GetDelta_UnstTime()*(start_iter+iSample*SamplingFreq);

 //unsigned long start_iter  =  config->GetUnst_RestartIter();
    unsigned long start_iter = config->GetRestart_Iter();

 su2double dt= config->GetDelta_UnstTime()*((iSample+iLocSample)*SamplingFreq);
 //su2double dt  = config->GetDelta_UnstTime()*(iSample*SamplingFreq);
 unsigned long iPanel,iDim;

 if (iSample==0 && iLocSample==0){
	x00= new su2double [nPanel]; y00= new su2double [nPanel]; z00= new su2double [nPanel];
	for (iPanel=0; iPanel<nPanel; iPanel++){
		x00[iPanel]=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 0 ];
		y00[iPanel]=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 1 ];
		z00[iPanel]=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 2 ];
	}
	//Read_surface_csv_files( config, geometry, iSample, iLocSample+1, 1);
 }	

if ( config->GetGrid_Movement() ) {
        if (config->GetKind_GridMovement() == RIGID_MOTION){
                //w[0]=config->GetRotation_Rate_X(iZone);
                //w[1]=config->GetRotation_Rate_Y(iZone);
                //w[2]=config->GetRotation_Rate_Z(iZone);

            w[0]=config->GetRotation_Rate(0);
            w[1]=config->GetRotation_Rate(1);
            w[2]=config->GetRotation_Rate(2);

		for (iDim=0;iDim<nDim;iDim++) alpha[iDim]=w[iDim]*dt;
		cout<<"alpha angles: "<<alpha[0]<<" "<<alpha[1]<<" "<<alpha[2]<<" dt: "<<dt<<endl;
	}
 }


for (iPanel=0; iPanel<nPanel; iPanel++){
 
 if ( config->GetGrid_Movement() ) {
 	//if (config->GetKind_GridMovement(iZone) == RIGID_MOTION){
 	if (config->GetKind_GridMovement() == RIGID_MOTION){
		x[0]=(z00[iPanel]*(sin(alpha[0])*sin(alpha[2]) + cos(alpha[0])*cos(alpha[2])*sin(alpha[1])) - y00[iPanel]*(cos(alpha[0])*sin(alpha[2]) - cos(alpha[2])*sin(alpha[0])*sin(alpha[1])) + x00[iPanel]*cos(alpha[1])*cos(alpha[2]));
		x[1]=(y00[iPanel]*(cos(alpha[0])*cos(alpha[2]) + sin(alpha[0])*sin(alpha[1])*sin(alpha[2])) - z00[iPanel]*(cos(alpha[2])*sin(alpha[0]) - cos(alpha[0])*sin(alpha[1])*sin(alpha[2])) + x00[iPanel]*cos(alpha[1])*sin(alpha[2]));
        	x[2]=(z00[iPanel]*cos(alpha[0])*cos(alpha[1]) - x00[iPanel]*sin(alpha[1]) + y00[iPanel]*cos(alpha[1])*sin(alpha[0]));	
		u[0]=-w[2]*x[1]+w[1]*x[2]-U1;
		u[1]= w[2]*x[0]-w[0]*x[2]-U2;
		u[2]=-w[1]*x[0]+w[0]*x[1]-U3;
	}else{
		x[0]=x00[iPanel]; x[1]=y00[iPanel]; x[2]=z00[iPanel];
		u[0]=-U1; u[1]=-U2; u[2]=-U3;
	}
 }else{
	x[0]=x00[iPanel]; x[1]=y00[iPanel];  x[2]=z00[iPanel];
        u[0]=-U1; u[1]=-U2; u[2]=-U3;
 }
 
 //cout<<"Coords for iPanel="<<iPanel<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;

 for (iDim=0; iDim<nDim; iDim++){
 	surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + iDim ]=x[iDim];
	          Q[iLocSample*nPanel*nDim + iPanel*nDim + iDim ]=u[iDim];
		//cout<<"u["<<iDim<<"]: "<<u[iDim];
 }
 //cout<<endl;
}//iPanel Loop

}//End of the subroutine PrescribedMotion






void FWHSolver::F1A_Formulation ( CConfig *config,unsigned long iObserver, unsigned long iPanel,  unsigned long iSample, unsigned long iLocSample, unsigned long i){

	#ifdef HAVE_MPI
        	int rank, nProcessor;
        	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        	MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
 	#endif
	unsigned long iDim;
	su2double Lnr,LnM, L_dotr, Qn, n_dot_Q, Qdotn, M, Mr,Mr1, M_dotr,n_dot,M_dot,L_dot;
 	su2double K,R,r,dS;
	su2double dt_src;
 	dt_src=config->GetDelta_UnstTime()*SamplingFreq;
	Lnr=0;LnM=0; L_dotr=0; Qn=0; n_dot_Q=0;  Qdotn=0; M=0; Mr=0; M_dotr=0;
	dS = AREA[iPanel];
  	for (iDim=0; iDim<nDim; iDim++){
		r = RadVec[5*i+iDim];
                M_dot=(Q[(iLocSample+1)*nPanel*nDim+iPanel*nDim+iDim]-Q[(iLocSample-1)*nPanel*nDim+iPanel*nDim+iDim])/2.0/dt_src/a_inf;
		M_dotr=M_dotr+M_dot*r;
                Mr=Mr+(Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*r)/a_inf;
		Mr1=1-Mr;
                M=M + Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim]/a_inf/a_inf;
		Qdotn=Qdotn+M_dot*a_inf*n[iLocSample*nPanel*nDim+iPanel*nDim+iDim];
		n_dot=(n[(iLocSample+1)*nPanel*nDim+iPanel*nDim+iDim]-n[(iLocSample-1)*nPanel*nDim+iPanel*nDim+iDim])/2.0/dt_src;
		n_dot_Q=n_dot_Q+n_dot*Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim];
		Qn=Qn+Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*n[iLocSample*nPanel*nDim+iPanel*nDim+iDim];
                L_dot=(F[(iLocSample-1)*nPanel+iPanel]*n[(iLocSample+1)*nPanel*nDim+iPanel*nDim+iDim]-F[(iLocSample-1)*nPanel+iPanel]*n[(iLocSample-1)*nPanel*nDim+iPanel*nDim+iDim])/2.0/dt_src;
		L_dotr= L_dotr+ L_dot*r;
                Lnr     =Lnr    + F[iLocSample*nPanel+iPanel]*n[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*r;
                LnM     =LnM    + F[iLocSample*nPanel+iPanel]*n[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim]/a_inf;
   	}
        M  = sqrt( M );
	R  = RadVec[5*i+nDim];
        K=M_dotr*R+Mr*a_inf-M*M*a_inf;
        T1=(FreeStreamDensity*(Qdotn+n_dot_Q))/R/Mr1/Mr1*(dS/4.0/M_PI);
        T2= (FreeStreamDensity*Qn*K)/R/R/Mr1/Mr1/Mr1*(dS/4.0/M_PI);
        T3= L_dotr/R/Mr1/Mr1/a_inf*(dS/4.0/M_PI);
        T4= (Lnr-LnM)/R/R/Mr1/Mr1*(dS/4.0/M_PI);
        T5= (Lnr*K)/R/R/Mr1/Mr1/Mr1/a_inf*(dS/4.0/M_PI);
	pp_t= (T1+T2+T3+T4+T5);
	
	if ((isnan(T1+T2+T3+T4+T5))||(isnan(T1))){
	//if ((iSample>nSample-2) && (iPanel>(nPanel-3))){
		cout<<"Here we have a NAN variable!!!!"<<endl;
		cout<<"iGlobal ID-"<<PointID[iPanel]<<endl;
        	cout<<"iSample- "<<iSample<<" iPanel-"<<iPanel<<" pp_t: "<<pp_t<<" T1: "<<T1<<" T2: "<<T2<<" T3: "<<T3<<" T4: "<<T4<<" T5: "<<T5<<endl;
 		cout<<"normals iLocSample-1: "<<n[(iLocSample-1)*nPanel*nDim+iPanel*nDim+0]<<" "<<n[(iLocSample-1)*nPanel*nDim+iPanel*nDim+1]<<" "<<n[(iLocSample-1)*nPanel*nDim+iPanel*nDim+2]<<" "<<endl;
		cout<<"normals iLocSample  : "<<n[(iLocSample)*nPanel*nDim+iPanel*nDim+0]<<" "<<n[(iLocSample)*nPanel*nDim+iPanel*nDim+1]<<" "<<n[(iLocSample)*nPanel*nDim+iPanel*nDim+2]<<" "<<endl;
		cout<<"normals iLocSample+1: "<<n[(iLocSample+1)*nPanel*nDim+iPanel*nDim+0]<<" "<<n[(iLocSample+1)*nPanel*nDim+iPanel*nDim+1]<<" "<<n[(iLocSample+1)*nPanel*nDim+iPanel*nDim+2]<<" "<<endl;
		cout<<"R: "<<R<<" M: "<<M<<" K: "<<K<<" Qn: "<<Qn<<" Mdotr: "<<M_dotr<<" Mdot: "<<M_dot<<endl;
		cout<<"Vels for iLocSample-1: "<<Q[(iLocSample-1)*nPanel*nDim+iPanel*nDim+0]<<" "<<Q[(iLocSample-1)*nPanel*nDim+iPanel*nDim+1]<<" "<<Q[(iLocSample-1)*nPanel*nDim+iPanel*nDim+2]<<" "<<endl;
		cout<<"Vels for iLocSample  : "<<Q[(iLocSample)*nPanel*nDim+iPanel*nDim+0]<<" "<<Q[(iLocSample)*nPanel*nDim+iPanel*nDim+1]<<" "<<Q[(iLocSample)*nPanel*nDim+iPanel*nDim+2]<<" "<<endl;
		cout<<"Vels for iLocSample+1: "<<Q[(iLocSample+1)*nPanel*nDim+iPanel*nDim+0]<<" "<<Q[(iLocSample+1)*nPanel*nDim+iPanel*nDim+1]<<" "<<Q[(iLocSample+1)*nPanel*nDim+iPanel*nDim+2]<<" "<<endl;
		cout<<"Pressure: "<<F[iLocSample*nPanel+iPanel]<<endl;
		cout<<"Area    : "<<dS<<endl;
		cout<<"r-vector: "<<RadVec[5*i+0]<<" "<<RadVec[5*i+1]<<" "<<RadVec[5*i+2]<<" "<<endl;
		cout<<"Qdotn: "<<Qdotn<<" n_dot_Q: "<<n_dot_Q<<" n_dot: "<<n_dot<<" Mr: "<<Mr<<endl;
		cout<<"iLocSample-1 x: "<<surface_geo[(iLocSample-1)*nPanel*nDim + iPanel*nDim + 0]<<" y: "<<surface_geo[(iLocSample-1)*nPanel*nDim + iPanel*nDim + 1]<<" z: "<<surface_geo[(iLocSample-1)*nPanel*nDim + iPanel*nDim + 2]<<endl;
		cout<<"iLocSample   x: "<<surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 0]<<" y: "<<surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 1]<<" z: "<<surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 2]<<endl;
		cout<<"iLocSample+1 x: "<<surface_geo[(iLocSample+1)*nPanel*nDim + iPanel*nDim + 0]<<" y: "<<surface_geo[(iLocSample-1)*nPanel*nDim + iPanel*nDim + 1]<<" z: "<<surface_geo[(iLocSample-1)*nPanel*nDim + iPanel*nDim + 2]<<endl;
	}



}//End of the subroutine F1A






void FWHSolver::F1A_SourceTimeDominant( CConfig *config, CGeometry *geometry){

 #ifdef HAVE_MPI
	int rank, nProcessor;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
 #endif
 su2double PP[3]= {0.0,0.0,0.0};
 su2double TT1[3]={0.0,0.0,0.0};
 su2double TT2[3]={0.0,0.0,0.0};
 su2double TT3[3]={0.0,0.0,0.0};
 su2double TT4[3]={0.0,0.0,0.0};
 su2double tau[3]= {0.0,0.0,0.0};
 unsigned long iObserver, iSample, iSample2, iDim, iLocSample, iPanel,iqSample, Smallest_iLocSample;
 su2double delta, dtt;
 su2double pp_n, TT1_n,TT2_n,TT3_n,TT4_n;
 su2double r1,r2,r3,R;
 su2double pp_mean, SPL_iObserver;
 char buffer[32];
 unsigned long res_factor=1;
 su2double dt_src,t;
 unsigned long iObserver_inner, iObserver_outer;
 unsigned long nObserver_inner, nObserver_outer;
 unsigned long Pre_Mot=0;
 dt_src=config->GetDelta_UnstTime()*SamplingFreq;
 ofstream pp_FWH_file ;
 ofstream pp_fft_file ;
 if (ceil(log2(nSample))!=floor(log2(nSample))){
        if (rank==MASTER_NODE){
                cout<<"Number of samples are not power of 2 ---> no FFT computation "<<endl;
        }
 }
 su2double dTT = config->GetDelta_UnstTime()*SamplingFreq; //Time Step to calculate sampling frequency
 su2double Fs = 1/dTT; //Sampling frequency to be used in fft calculation

 su2double **t_interp_inner= NULL,**pp_TimeDomain_inner = NULL,**pp_TimeDomainGlobal_inner = NULL;
 su2double **T1_I=NULL,**T2_I=NULL,**T3_I=NULL,**T4_I=NULL;
 su2double *t_interp_outer = NULL,*pp_TimeDomain_outer = NULL,*pp_TimeDomainGlobal_outer = NULL;
 su2double *T1_O=NULL,*T2_O=NULL,*T3_O=NULL,*T4_O=NULL;


 if(config->GetAcoustic_Inner_ObsLoop()){
	
	delete[] t_interp_outer;
        delete[] pp_TimeDomain_outer;
        delete[] pp_TimeDomainGlobal_outer;
 	t_interp_inner= new su2double *[nObserver];
 	pp_TimeDomain_inner=new su2double *[nObserver];
 	if (rank==MASTER_NODE) pp_TimeDomainGlobal_inner=new su2double *[nObserver];
	T1_I=new su2double *[nObserver];
	T2_I=new su2double *[nObserver];
	T3_I=new su2double *[nObserver];
	T4_I=new su2double *[nObserver];
 	for (iObserver=0; iObserver<nObserver; iObserver++){
 		t_interp_inner[iObserver]=new su2double [nSample];
        	pp_TimeDomain_inner[iObserver]=new su2double [nSample];
        	if (rank==MASTER_NODE) pp_TimeDomainGlobal_inner[iObserver]=new su2double [nSample];
		T1_I[iObserver]=new su2double [nSample];
		T2_I[iObserver]=new su2double [nSample];
		T3_I[iObserver]=new su2double [nSample];
		T4_I[iObserver]=new su2double [nSample];
        	for (iSample=0; iSample<nSample; iSample++){
        		t_interp_inner[iObserver][iSample]=0.0;
                	pp_TimeDomain_inner[iObserver][iSample]=0.0;
                	if (rank==MASTER_NODE) pp_TimeDomainGlobal_inner[iObserver][iSample]=0.0;
			T1_I[iObserver][iSample]=0.0;
			T2_I[iObserver][iSample]=0.0;
			T3_I[iObserver][iSample]=0.0;
			T4_I[iObserver][iSample]=0.0;
     		}
 	}
	
 }else{
	delete[] t_interp_inner;
	delete[] pp_TimeDomain_inner;
	delete[] pp_TimeDomainGlobal_inner;
	t_interp_outer= new su2double [nSample];
        pp_TimeDomain_outer=new su2double [nSample];
	if (rank==MASTER_NODE) pp_TimeDomainGlobal_outer=new su2double [nSample];
	T1_O=new su2double [nSample];
	T2_O=new su2double [nSample];
	T3_O=new su2double [nSample];
	T4_O=new su2double [nSample];
 }

 if(config->GetAcoustic_Inner_ObsLoop()){
 	for (iObserver=0; iObserver<nObserver; iObserver++){
		t_interp_inner[iObserver][0]=StartTime[iObserver];
 	}
 }
 
 if(config->GetAcoustic_Inner_ObsLoop()){
	nObserver_inner=nObserver;
	nObserver_outer=1;
 }else{
	nObserver_outer=nObserver;
	nObserver_inner=1;
 }

 for (iObserver_outer=0; iObserver_outer<nObserver_outer; iObserver_outer++){	
	if(!(config->GetAcoustic_Inner_ObsLoop())){
		for (iSample=0; iSample<nSample; iSample++){
                	t_interp_outer[iSample]=0.0;
                	pp_TimeDomain_outer[iSample]=0.0;
                	if (rank==MASTER_NODE) pp_TimeDomainGlobal_outer[iSample]=0.0;
			T1_O[iSample]=0.0;
			T2_O[iSample]=0.0;
			T3_O[iSample]=0.0;
			T4_O[iSample]=0.0;
        	}
		t_interp_outer[0]=StartTime[iObserver_outer];
	}

	//if ((config->GetOutput_FileFormat() == PARAVIEW)) {
        	pp_out = new su2double [nPanel];
        	for (iPanel=0; iPanel<nPanel; iPanel++){
                	pp_out[iPanel]=0.0;
        	}
 	//}


 	iSample=0;
 	for (iLocSample=0; iLocSample<nqSample; iLocSample++){
		for (iPanel=0; iPanel < nPanel; iPanel++){
			for (iDim=0; iDim < nDim; iDim++){
				n[iLocSample*nPanel*nDim + iPanel*nDim + iDim ]=0.0;
			}
		}
		if(config->GetAcoustic_Prescribed_Motion()) Pre_Mot=1;
		Read_surface_csv_files( config, geometry, iSample, iLocSample, Pre_Mot);
		ComputeNorms(config, geometry, iSample, iLocSample,iObserver_outer);
	}
	if(!(config->GetAcoustic_Prescribed_Motion())){
 		for (iLocSample=1; iLocSample<nqSample-1; iLocSample++){
			ComputeVelocities(config, iSample, iLocSample);
 		}
	}
 	iLocSample=nqSample/2;
 	iSample2=0;
 	for (iSample=0; iSample<nSample; iSample++){
		if (rank==MASTER_NODE) cout<<"iSample-"<<iSample<<endl;
		Smallest_iLocSample=999999;
		for (iPanel=0; iPanel<nPanel; iPanel++){
			for (iObserver_inner=0; iObserver_inner<nObserver_inner; iObserver_inner++){
				if(config->GetAcoustic_Inner_ObsLoop()){
					iObserver=iObserver_inner;
					t=t_interp_inner[iObserver][iSample];
				}else{
					iObserver=iObserver_outer;
					t=t_interp_outer[iSample];
				}
				for (int i=0; i<3 ; i++){
					//cout<<"Hey-15 iSample-"<<iSample<<" i: "<<i<<" iLocSample: "<<iLocSample<<endl;
                        		ComputeObserverTime(config, iObserver, iPanel, iSample2, iLocSample-1+i,i);
                         		tau[i] = RadVec[5*i+4];
                  		}
				//cout<<"Hey-17 iSample-"<<iSample<<" iPanel="<<iPanel<<" t_interp:"<<t<<" taus: "<<tau[0]<<" "<<tau[1]<<" "<<tau[2]<<" iLocSample-"<<iLocSample<<endl;		
				while ( (t<tau[1]) || (t>tau[2]) ){
					if ( t<tau[1] ){
						if ( (iLocSample<3)  ) goto endloop;
						/*update tau and shift down*/
						iLocSample=iLocSample-1;
					}	
					if ( t>tau[2] ){
						if ( (iLocSample>nqSample-4) ) goto endloop;
						iLocSample=iLocSample+1;
					}
				
					for (int i=0; i<3 ; i++){
						ComputeObserverTime(config, iObserver, iPanel, iSample2, iLocSample-1+i,i);
						tau[i] = RadVec[5*i+4];
					}
				}
				endloop:
				if (iLocSample<Smallest_iLocSample){
                        		Smallest_iLocSample=iLocSample;
                      		}
				if (iPanel==0){
                                	if (iSample<nSample-1){
						if(config->GetAcoustic_Inner_ObsLoop()){
                                        		t_interp_inner[iObserver][iSample+1]=t_interp_inner[iObserver][iSample]+dt[iObserver];
						}else{
							t_interp_outer[iSample+1]=t_interp_outer[iSample]+dt[iObserver];
						}
                                	}
                        	}
				for (int i=0; i<3 ; i++){
					F1A_Formulation ( config,iObserver, iPanel, iSample,  iLocSample-1+i, i);
					PP[i]=pp_t;
					TT1[i]=T1;TT2[i]=T2;TT3[i]=T3;TT4[i]=T4+T5;
				}
				delta=t-tau[1];
				dtt=tau[2]-tau[1];
				/*2nd order taylor expansion for polynamial interpolation*/
				pp_n=PP[1]+delta*(PP[2]-PP[0])*0.5/dtt+delta*delta*0.5*(PP[2]-2*PP[1]+PP[0])/dtt/dtt;
				TT1_n=TT1[1]+delta*(TT1[2]-TT1[0])*0.5/dtt+delta*delta*0.5*(TT1[2]-2*TT1[1]+TT1[0])/dtt/dtt;
				TT2_n=TT2[1]+delta*(TT2[2]-TT2[0])*0.5/dtt+delta*delta*0.5*(TT2[2]-2*TT2[1]+TT2[0])/dtt/dtt;
				TT3_n=TT3[1]+delta*(TT3[2]-TT3[0])*0.5/dtt+delta*delta*0.5*(TT3[2]-2*TT3[1]+TT3[0])/dtt/dtt;
				TT4_n=TT4[1]+delta*(TT4[2]-TT4[0])*0.5/dtt+delta*delta*0.5*(TT4[2]-2*TT4[1]+TT4[0])/dtt/dtt;
				if(config->GetAcoustic_Inner_ObsLoop()){			
					pp_TimeDomain_inner[iObserver][iSample] = pp_TimeDomain_inner[iObserver][iSample] + pp_n;
					T1_I[iObserver][iSample]=T1_I[iObserver][iSample]+TT1_n;
					T2_I[iObserver][iSample]=T2_I[iObserver][iSample]+TT2_n;
					T3_I[iObserver][iSample]=T3_I[iObserver][iSample]+TT3_n;
					T4_I[iObserver][iSample]=T4_I[iObserver][iSample]+TT4_n;
				}else{
					pp_TimeDomain_outer[iSample] = pp_TimeDomain_outer[iSample] + pp_n;
					T1_O[iSample]=T1_O[iSample]+TT1_n;
					T2_O[iSample]=T2_O[iSample]+TT2_n;
					T3_O[iSample]=T3_O[iSample]+TT3_n;
					T4_O[iSample]=T4_O[iSample]+TT4_n;	
				}

				//if ((config->GetOutput_FileFormat() == PARAVIEW)) {
					pp_out[iPanel]=pp_n;
					R=RadVec[5*1+3]+delta*(RadVec[5*2+3]-RadVec[5*0+3])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+3]-2*RadVec[5*1+3]+RadVec[5*0+3])/dtt/dtt;
					r1=RadVec[5*1+0]+delta*(RadVec[5*2+0]-RadVec[5*0+0])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+0]-2*RadVec[5*1+0]+RadVec[5*0+0])/dtt/dtt;
					r2=RadVec[5*1+1]+delta*(RadVec[5*2+1]-RadVec[5*0+1])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+1]-2*RadVec[5*1+1]+RadVec[5*0+1])/dtt/dtt;
					r3=RadVec[5*1+2]+delta*(RadVec[5*2+2]-RadVec[5*0+2])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+2]-2*RadVec[5*1+2]+RadVec[5*0+2])/dtt/dtt;
					RetSurf[iPanel*nDim+0]=-(r1*R)+Observer_Locations[iObserver][0]-U1*t;
					RetSurf[iPanel*nDim+1]=-(r2*R)+Observer_Locations[iObserver][1]-U2*t;
					RetSurf[iPanel*nDim+2]=-(r3*R)+Observer_Locations[iObserver][2]-U3*t;			
				//}

			}//iObserver_inner Loop
		}//iPanel Loop
		//cout<<"rank"<<rank<<"hey-25"<<endl;	
		//PARAVIEW Section for Output

		//if ((config->GetOutput_FileFormat() == PARAVIEW)) {
                	if (rank==MASTER_NODE)  cout<<"Extracting paraview file (.vtk) for post-processing"<<endl;
                	Paraview_Output(config,iObserver,iSample, 0, 0, 1);
                	if (rank==MASTER_NODE)  cout<<"End of paraview extraction"<<endl;
        	//}//if Paraview

		#ifdef HAVE_MPI
       			SU2_MPI::Allreduce(&Smallest_iLocSample,&Smallest_iLocSample,1,MPI_UNSIGNED_LONG,MPI_MIN,MPI_COMM_WORLD);
      		#endif
		if ( (Smallest_iLocSample>5) ){
        		for (iqSample=0; iqSample<nqSample-1; iqSample++){
				for (iPanel=0; iPanel<nPanel; iPanel++){
                			for (iDim=0; iDim<nDim; iDim++){
						if(config->GetAcoustic_Prescribed_Motion()){
							Q[iqSample*nPanel*nDim + iPanel*nDim + iDim]=Q[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];
						}else{
							if ( (iqSample>0) && (iqSample<nqSample-2) && (!(config->GetAcoustic_Prescribed_Motion())) ) {
                        					Q[iqSample*nPanel*nDim + iPanel*nDim + iDim]=Q[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];
							}
						}
						n[iqSample*nPanel*nDim + iPanel*nDim + iDim]=n[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];
						if (iqSample==nqSample-2) n[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim]=0;
                                		surface_geo[iqSample*nPanel*nDim + iPanel*nDim + iDim]=surface_geo[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];
               				}
					F[iqSample*nPanel+iPanel]=F[(iqSample+1)*nPanel+iPanel];	
				}
         		}
			iSample2++;
			if (iSample2< (nSample-nqSample)){
				if(config->GetAcoustic_Prescribed_Motion()) Pre_Mot=1;
				Read_surface_csv_files( config, geometry, iSample2, nqSample-1, Pre_Mot);
				ComputeNorms(config, geometry, iSample2, nqSample-1,iObserver_outer);
			}
			if ( (iSample2< (nSample-nqSample-1)) && (!(config->GetAcoustic_Prescribed_Motion())) ){
				ComputeVelocities(config, iSample2, nqSample-2);
			}
		}

 	}//iSample Loop		
				
 	#ifdef HAVE_MPI
		for (iObserver_inner=0; iObserver_inner<nObserver_inner; iObserver_inner++){
			for (iSample=0; iSample<nSample; iSample++){
				if(config->GetAcoustic_Inner_ObsLoop()){
        				SU2_MPI::Reduce( &pp_TimeDomain_inner[iObserver_inner][iSample] , &pp_TimeDomainGlobal_inner[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM ,MASTER_NODE, MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &t_interp_inner[iObserver_inner][iSample] , &t_interp_inner[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T1_I[iObserver_inner][iSample] , &T1_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T2_I[iObserver_inner][iSample] , &T2_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T3_I[iObserver_inner][iSample] , &T3_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T4_I[iObserver_inner][iSample] , &T4_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
				}else{
					SU2_MPI::Reduce( &pp_TimeDomain_outer[iSample] , &pp_TimeDomainGlobal_outer[iSample] , 1 , MPI_DOUBLE , MPI_SUM ,MASTER_NODE, MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &t_interp_outer[iSample] , &t_interp_outer[iSample] , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T1_O[iSample] , &T1_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T2_O[iSample] , &T2_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T3_O[iSample] , &T3_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
					SU2_MPI::Allreduce( &T4_O[iSample] , &T4_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
				}
			}
		}
 	#endif

 	if (rank==MASTER_NODE){
		for (iObserver_inner=0; iObserver_inner<nObserver_inner; iObserver_inner++){
			pp_mean=0;
   			for (iSample=0; iSample<nSample; iSample++){
				if(config->GetAcoustic_Inner_ObsLoop()){
					iObserver=iObserver_inner;
            				pp_mean = (pp_TimeDomainGlobal_inner[iObserver][iSample]+pp_mean*iSample)/(iSample+1);
				}else{
					iObserver=iObserver_outer;
					pp_mean = (pp_TimeDomainGlobal_outer[iSample]+pp_mean*iSample)/(iSample+1);
				}
                	}
                	if (ceil(log2(nSample))==floor(log2(nSample))){
                       		cout<<"FFT Transformation to obtain dB-Freq[Hz] "<<endl;
				if(config->GetAcoustic_Inner_ObsLoop()){
                       			FFT_AcousticPressureSignal(pp_mean,pp_TimeDomainGlobal_inner[iObserver],res_factor);
					snprintf(buffer, sizeof(char) * 32, "pp_FWH_fft_%03d", (iObserver_inner+1));
				}else{
					FFT_AcousticPressureSignal(pp_mean,pp_TimeDomainGlobal_outer,res_factor);
                       			snprintf(buffer, sizeof(char) * 32, "pp_FWH_fft_%03d", (iObserver_outer+1));
                       		}
				pp_fft_file.open(buffer);
                       		cout<<"Frequency[Hz] vs dB data can be found in "<<buffer<<" for Observer-"<<iObserver+1<<endl;
                	}

                	SPL_iObserver = 0.0;
                	snprintf(buffer, sizeof(char) * 32, "pp_FWH_%03d_Zone_%i", (iObserver+1), (iZone));
                	pp_FWH_file.open(buffer);

                	for (iSample=0; iSample<nSample; iSample++){
                        	if (iSample==0){
                                	pp_FWH_file <<  "Time"<<' '<<"P-Fluctuation"<<' '<<"Acoustic-Pressure"<<' '<<"Thickness"<<' '<<"Loading"<<' '<<"Term-1"<<' '<<"Term-2"<<' '<<"Term-3"<<' '<<"Term-4"<<endl;
                        	}
				if(config->GetAcoustic_Inner_ObsLoop()){
                        		pp_FWH_file << std::setprecision(15) <<t_interp_inner[iObserver][iSample] <<' '<<pp_TimeDomainGlobal_inner[iObserver][iSample]-pp_mean<<' '<<pp_TimeDomainGlobal_inner[iObserver][iSample]<< ' ' <<T1_I[iObserver][iSample]+T2_I[iObserver][iSample]<<' '<< T3_I[iObserver][iSample]+T4_I[iObserver][iSample]<<' '<<T1_I[iObserver][iSample]<<' '<<T2_I[iObserver][iSample]<<' '<<T3_I[iObserver][iSample]<<' '<<T4_I[iObserver][iSample]<<' '<<endl;
                        	}else{
					pp_FWH_file << std::setprecision(15) <<t_interp_outer[iSample] <<' '<<pp_TimeDomainGlobal_outer[iSample]-pp_mean<<' '<<pp_TimeDomainGlobal_outer[iSample]<< ' ' <<T1_O[iSample]+T2_O[iSample]<<' '<<T3_O[iSample]+T4_O[iSample]<<' '<<T1_O[iSample]<<' '<<T2_O[iSample]<<' '<<T3_O[iSample]<<' '<<T4_O[iSample]<<' '<<endl;
				}
				if (ceil(log2(nSample))==floor(log2(nSample))){
                                	if (iSample==0){
                                       		pp_fft_file << "Frequency[Hz]"<<" "<<"dB"<<endl;
                                	}
                                	if ((iSample<(nSample/2+1)) && (iSample>0)){
                                       		pp_fft_file << std::setprecision(15) << Fs*iSample/nSample<< ' '<< real(pp_fft[iSample])<< ' '<<endl;
                                	}
                        	}
				if(config->GetAcoustic_Inner_ObsLoop()){
                        		SPL_iObserver= SPL_iObserver + (pp_TimeDomainGlobal_inner[iObserver][iSample]-pp_mean)*(pp_TimeDomainGlobal_inner[iObserver][iSample]-pp_mean);
               			}else{
					SPL_iObserver= SPL_iObserver + (pp_TimeDomainGlobal_outer[iSample]-pp_mean)*(pp_TimeDomainGlobal_outer[iSample]-pp_mean);
				}
			}
                	SPL_iObserver = sqrt(SPL_iObserver/nSample);
                	SPL = SPL + SPL_iObserver;
                	pp_FWH_file.close();
                	pp_fft_file.close();
                	cout<<std::setprecision(15)<<"RMS(p')-FWH over observer#"<<iObserver+1<<" -------------------> "<<SPL_iObserver <<"  *******"<<endl;
                	cout<<"F1A Computation is finished for Observer#"<<iObserver+1<<". Tabulated data (Time-domain) can be found in "<< buffer<<" ."<<endl;
		}
 	}
 }//iObserver_Outer Loop
 if (rank==MASTER_NODE){
	SPL = SPL/nObserver;
	cout<<endl<<std::setprecision(15)<<"****** RMS(p') averaged over "<<nObserver<<" observer locations = "<<SPL <<"  **************"<<endl;
 }


}//End of the function F1A_SourceTimeDominant



void FWHSolver::FFT_AcousticPressureSignal(su2double pp_mean, su2double *pp_TimeDomain, unsigned long res_factor){
	unsigned long iPanel, iSample;
	FFT* FFT_container = new FFT() ;
	su2double pLa;

	/*perform FFT on Acoustic Pressure Signal*/
	for (iSample=0; iSample<nSample/res_factor; iSample++){
		pp_fft[iSample].real(pp_TimeDomain[iSample*res_factor]-pp_mean);
     	}	
	/* Writing the complex array data*/
	CArray dataPP(pp_fft ,nSample/res_factor);
	/* Call the FFT function for performing the forward fft */
	FFT_container->fft_r2(dataPP);
	for (iSample=0; iSample<nSample/res_factor; iSample++){
		pLa = abs(real(dataPP[iSample]));
		pp_fft[iSample] = 20*log10(pLa/0.00002);
	}
}






void FWHSolver::Write_Sensitivities(CSolver *solver, CConfig *config, CGeometry *geometry){
    

    unsigned long iVar, iPanel, iSample,iExtIter,Max_nPanel, Tot_nPanel,nVar,Global_Index ;
      ofstream CAA_AdjointFile;
    char buffer [50];

    //start_iter  =  config->GetUnst_RestartIter();
    unsigned long start_iter = config->GetRestart_Iter();

    int rank, iProcessor, nProcessor, iFace;

//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

    unsigned long Buffer_Send_nPanel[1], *Buffer_Recv_nPanel = NULL;
    if (rank == MASTER_NODE) Buffer_Recv_nPanel= new unsigned long [nProcessor];

     Buffer_Send_nPanel[0]=nPanel;
#ifdef HAVE_MPI
      SU2_MPI::Gather(&Buffer_Send_nPanel, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPanel, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD); //send the number of vertices at each process to the master
      SU2_MPI::Allreduce(&nPanel,&Max_nPanel,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
      SU2_MPI::Reduce(&nPanel,&Tot_nPanel,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (panels)
#endif

      nVar = nDim+3;


      /* Loop through all time steps */
  for (iSample=0; iSample<nSample; iSample++){
      iExtIter = start_iter + iSample*SamplingFreq;

      /* pack sensitivity values in each processor and send to root */
      su2double *Buffer_Send_dJdU = new su2double [Max_nPanel*nVar];
      su2double *Buffer_Send_dJdX = new su2double [Max_nPanel*nDim];
      unsigned long *Buffer_Send_GlobalIndex = new unsigned long [Max_nPanel];
      //zero send buffers
      for (int i=0; i <Max_nPanel*nVar; i++){
       Buffer_Send_dJdU[i]=0.0;
      }
      for (int i=0; i <Max_nPanel*nDim; i++){
       Buffer_Send_dJdX[i]=0.0;
      }
      for (int i=0; i <Max_nPanel; i++){
       Buffer_Send_GlobalIndex[i]=0;
      }
      su2double *Buffer_Recv_dJdU = NULL;
      su2double *Buffer_Recv_dJdX = NULL;
      unsigned long *Buffer_Recv_GlobalIndex = NULL;

      if (rank == MASTER_NODE) {
       Buffer_Recv_dJdU = new su2double [nProcessor*Max_nPanel*nVar];
       Buffer_Recv_dJdX = new su2double [nProcessor*Max_nPanel*nDim];
       Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*Max_nPanel];
      }

      for (iVar=0; iVar<nVar; iVar++){
          for(iPanel=0; iPanel<nPanel; iPanel++){
              Buffer_Send_dJdU[iVar*nPanel+iPanel] = dJdU[iVar][iPanel][iSample];
	      if (iVar<3){
	      	Buffer_Send_dJdX[iVar*nPanel+iPanel] = dJdX[iVar][iPanel][iSample];

	      }
            }
        }

      for (iPanel=0; iPanel<nPanel; iPanel++){
         Buffer_Send_GlobalIndex[iPanel] = PointID[iPanel];
        }


#ifdef HAVE_MPI
     SU2_MPI::Gather(Buffer_Send_dJdU, Max_nPanel*nVar, MPI_DOUBLE, Buffer_Recv_dJdU,  Max_nPanel*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
     SU2_MPI::Gather(Buffer_Send_dJdX, Max_nPanel*nDim, MPI_DOUBLE, Buffer_Recv_dJdX,  Max_nPanel*nDim , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
     SU2_MPI::Gather(Buffer_Send_GlobalIndex,Max_nPanel, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, Max_nPanel , MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

     	/* root opens a file at each time step and write out the merged dJdU values at that time step into the file */
      if (rank == MASTER_NODE){
      	char cstr [200];

      	SPRINTF (cstr, "Adj_CAA");
      	if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
      	if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(iExtIter));
      	if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(iExtIter));
      	if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(iExtIter));
      	if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
      	strcat (cstr, buffer);
      	cout<<cstr<<endl;
      	CAA_AdjointFile.precision(15);
      	CAA_AdjointFile.open(cstr, ios::out);

      	/*--- Loop through all of the collected data and write each node's values ---*/
      	unsigned long Total_Index, Total_Index_2;
      	for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iPanel = 0; iPanel < Buffer_Recv_nPanel[iProcessor]; iPanel++) {
            	Global_Index = Buffer_Recv_GlobalIndex[iProcessor*Max_nPanel+iPanel ];
            	CAA_AdjointFile  << scientific << Global_Index << "\t";

           	for (iVar = 0; iVar < nVar+nDim; iVar++){

          		/*--- Current index position and global index ---*/
     			Total_Index  = iProcessor*Max_nPanel*nVar + iVar*Buffer_Recv_nPanel[iProcessor]  + iPanel;
 
          		/*--- Write to file---*/
			if(iVar<nVar){ 
          			CAA_AdjointFile << scientific <<  Buffer_Recv_dJdU[Total_Index]   << "\t";
			}else{
				Total_Index_2= iProcessor*Max_nPanel*nDim + (iVar-nVar)*Buffer_Recv_nPanel[iProcessor]  + iPanel;
				CAA_AdjointFile << scientific <<  Buffer_Recv_dJdX[Total_Index_2]   << "\t";
			}
           	}
          	CAA_AdjointFile  << endl;

    	  }
   	}

      	CAA_AdjointFile.close();

	delete [] Buffer_Recv_dJdU;
      	delete [] Buffer_Recv_dJdX;
      	delete [] Buffer_Recv_GlobalIndex;
      }//if for masternode
      delete [] Buffer_Send_dJdU;
      delete [] Buffer_Send_dJdX;
      delete [] Buffer_Send_GlobalIndex;
  }


	if (rank == MASTER_NODE) cout<<"Calling Paraview"<<endl;
	Paraview_Output(config,0,iSample, 1, 1, 1);	



}



void FWHSolver::Paraview_Output(CConfig *config, unsigned long iObserver, unsigned long iSample, unsigned long SURFACE_TYPE, unsigned long DERIVATIVES, unsigned long FIA_TERMS){

/* User Guideline

Paraview_Output(SURFACE_TYPE, DERIVATIVES, FIA_TERMS)

SURFACE_TYPE:  1 (if sigma surface) ,     0 (if retarded surface)
DERIVATIVES :  1 (if derivatives asked),  0 (not asked)
FIA_TERMS   :  1 (if F1A Terms required), 0 (not required)

e.g. Paraview_Output(1, 0, 1) -->  F1A terms on sigma surfaces (no derivative terms)

*/

int iProcessor, nProcessor, iFace;

//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

unsigned long iVar, iPanel, start_iter,iExtIter,nVar,Global_Index ;
char buffer [50];
unsigned long Total_Index, Total_Index_2;
unsigned long iConn, iDim;


if (iSample==0){
 unsigned long iConnSize, nConnSize, nConnSizeMaxTot;
 unsigned long nFaceMaxTot ;

 if (rank == MASTER_NODE) nFaceSizeMaster   = new unsigned long [nProcessor];
 if (rank == MASTER_NODE) nConnSizeMaster   = new unsigned long [nProcessor];
 nPanelMaster      = new unsigned long [nProcessor];


 iConnSize=0;
 for (iFace = 0; iFace < nFace+nHaloFaces ; iFace++){
	if (iFace<nFace){
                iConnSize=iConnSize+Conn[iFace][0]+1;
        }else{
                iConnSize=iConnSize+HaloLayer2[iFace-nFace][0]+1;
        }
 }

 nConnSize=iConnSize;
 unsigned long nAllFaces = nFace+nHaloFaces;

 #ifdef HAVE_MPI
        SU2_MPI::Allreduce(&nConnSize,&nConnSizeMax,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD);//Find max value of connectivity
        SU2_MPI::Gather(&nConnSize, 1, MPI_UNSIGNED_LONG, nConnSizeMaster, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);//Gather each conn size into master node
        SU2_MPI::Reduce(&nFace,&nFaceTot,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);//Find total number of faces
        SU2_MPI::Allreduce(&nAllFaces,&nFaceMax,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD);//Find max value of the face sizes
        SU2_MPI::Gather(&nAllFaces, 1, MPI_UNSIGNED_LONG, nFaceSizeMaster, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);//bring face size into master
        SU2_MPI::Allreduce(&nConnSizeMax,&nConnSizeMaxTot,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
	SU2_MPI::Allreduce(&nConnSize,&nConnSizeTot,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&nFaceMax,&nFaceMaxTot,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
        SU2_MPI::Allgather(&nPanel, 1, MPI_UNSIGNED_LONG, nPanelMaster, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&nPanel,&Max_nPanel,1,MPI_UNSIGNED_LONG,MPI_MAX,MPI_COMM_WORLD); //find the max num of vertices over all processes
        SU2_MPI::Reduce(&nPanel,&Tot_nPanel,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of vertices (panels)
        SU2_MPI::Reduce(&nHaloFaces,&nHaloFacesTot,1,MPI_UNSIGNED_LONG,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD); //find the total num of halofaces
 #endif

 unsigned long *Conn_single= new unsigned long [nConnSizeMax];
 unsigned long *CellType= new unsigned long [nFaceMax];
 unsigned long ext_panel;
 unsigned long *ext_panel2 = new unsigned long [nProcessor];

 for (int iProcessor=0; iProcessor<nProcessor; iProcessor++){
        ext_panel=0;
        for (int i=0; i<iProcessor ; i++){
                ext_panel=ext_panel+nPanelMaster[i];
        }
        ext_panel2[iProcessor]=ext_panel;
 }

 iConnSize=0;
 for (iFace = 0; iFace < nFace+nHaloFaces ; iFace++){
        if (iFace<nFace){
                Conn_single[iConnSize]=Conn[iFace][0];
                iConnSize++;
                for (int i=0; i<Conn[iFace][0]; i++){
                        Conn_single[iConnSize]=Conn[iFace][i+1]-1+ext_panel2[rank];
                        iConnSize++;
                }
                if (Conn[iFace][0]==3){;
                        CellType[iFace]=5;
                }
                if (Conn[iFace][0]==4){;
                        CellType[iFace]=9;
                }
        }else{
                Conn_single[iConnSize]=HaloLayer2[iFace-nFace][0];
                iConnSize++;
                for (int i=0; i<HaloLayer2[iFace-nFace][0] ; i++){
                        Conn_single[iConnSize]=HaloLayer2[iFace-nFace][2*(i+1)]-1+ext_panel2[HaloLayer2[iFace-nFace][2*i+1]];
                        iConnSize++;
                }
                if (HaloLayer2[iFace-nFace][0]==3){;
                        CellType[iFace]=5;
                }
                if (HaloLayer2[iFace-nFace][0]==4){;
                        CellType[iFace]=9;
                }
        }
 }

 Conn_single_master = new unsigned long [nProcessor*nConnSizeMax];
 for (int i=0; i< (nProcessor*nConnSizeMax); i++){
        Conn_single_master[i]=0;
        if (nProcessor==1){
                Conn_single_master[i]=Conn_single[i];
        }
 }
 CellTypeMaster = new unsigned long [nFaceMaxTot];
 for (int i=0; i<nFaceMaxTot; i++){
        CellTypeMaster[i]=0;
        if (nProcessor==1){
                CellTypeMaster[i]=CellType[i];
        }
 }



 #ifdef HAVE_MPI
        SU2_MPI::Gather(Conn_single, nConnSizeMax, MPI_UNSIGNED_LONG, Conn_single_master, nConnSizeMax, MPI_UNSIGNED_LONG,MASTER_NODE,MPI_COMM_WORLD);
        SU2_MPI::Gather(CellType, nFaceMax, MPI_UNSIGNED_LONG, CellTypeMaster, nFaceMax, MPI_UNSIGNED_LONG,MASTER_NODE,MPI_COMM_WORLD);
 #endif

 nVar = nDim+3;

 surfaceCoordsMaster= new su2double [nProcessor*Max_nPanel*nDim];
 SingleMaster_dJdU  = new su2double [nProcessor*Max_nPanel*nVar];
 SingleMaster_dJdX  = new su2double [nProcessor*Max_nPanel*nDim];

 Single_dJdU  = new su2double [Max_nPanel*nVar];
 Single_dJdX  = new su2double [Max_nPanel*nDim];
 surfaceCoords= new su2double [Max_nPanel*nDim];


 for (int i=0; i <Max_nPanel*nVar; i++){
        Single_dJdU[i]=0;
        if (nVar<nDim){
                Single_dJdX[i]=0;
        }
 }


 SingleMaster_pp = new su2double [nProcessor*Max_nPanel];
 Single_pp  = new su2double [Max_nPanel];

 for (int i=0; i <Max_nPanel; i++){
        Single_pp[i]=0;
 }
}
//for (iSample=0; iSample<nSample; iSample++){

        for (iPanel=0; iPanel<nPanel ; iPanel++){
                for (iDim=0; iDim<nDim ; iDim++){
                        if (SURFACE_TYPE==1){//Sigma Surfaces
                                //surfaceCoords[nDim*iPanel+iDim]=surface_geo[iPanel][iSample][iDim];
                                //surfaceCoords[nDim*iPanel+iDim]=surface_geo[iPanel][0][iDim];
				surfaceCoords[nDim*iPanel+iDim]=surface_geo[0*nPanel*nDim + iPanel*nDim + iDim];
                        }else{//retarded surfaces
                                surfaceCoords[nDim*iPanel+iDim]=RetSurf[iPanel*nDim+iDim];
                        }
                }
        }

        #ifdef HAVE_MPI
                SU2_MPI::Gather(surfaceCoords,Max_nPanel*nDim,MPI_DOUBLE, surfaceCoordsMaster, Max_nPanel*nDim, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        #endif


        if (DERIVATIVES==1){
                for (iVar=0; iVar<nVar; iVar++){
                        for(iPanel=0; iPanel<nPanel; iPanel++){
                                Single_dJdU[iVar*nPanel+iPanel] = dJdU[iVar][iPanel][iSample];
                                if (iVar<3){
                                        Single_dJdX[iVar*nPanel+iPanel] = dJdX[iVar][iPanel][iSample];
                                }
                        }
                }

                #ifdef HAVE_MPI
                        SU2_MPI::Gather(Single_dJdU, Max_nPanel*nVar, MPI_DOUBLE, SingleMaster_dJdU,  Max_nPanel*nVar , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                        SU2_MPI::Gather(Single_dJdX, Max_nPanel*nDim, MPI_DOUBLE, SingleMaster_dJdX,  Max_nPanel*nDim , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                #endif
        }

        if (FIA_TERMS==1){
                for(iPanel=0; iPanel<nPanel; iPanel++){
                        //Single_pp[iPanel]=pp_interp[iSample][iPanel];
			Single_pp[iPanel]=pp_out[iPanel];
                        if (nProcessor==1){
                                //SingleMaster_pp[iPanel]=pp_interp[iSample][iPanel];
                                SingleMaster_pp[iPanel]=pp_out[iPanel];
                        }
                }
                #ifdef HAVE_MPI
                        SU2_MPI::Gather(Single_pp , Max_nPanel, MPI_DOUBLE, SingleMaster_pp,  Max_nPanel , MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                #endif
        }


        if (rank == MASTER_NODE){

                ofstream paraview_file;
                snprintf(buffer, sizeof(char) * 32, "Adj_out_%06d.vtk", (iSample+1));
                if (FIA_TERMS==1){
                        snprintf(buffer, sizeof(char) * 32, "F1A_iZ_%i_obs_%02d_t_%04d.vtk", (iZone),(iObserver),(iSample));
                }
                paraview_file.open(buffer);
                cout<<"Paraview output - "<<buffer<<endl;
                paraview_file<<"# vtk DataFile Version 3.0"<<endl;
                paraview_file<<"vtk output"<<endl;
                paraview_file<<"ASCII"<<endl;
                paraview_file<<"DATASET UNSTRUCTURED_GRID"<<endl;
                paraview_file<<"POINTS "<<Tot_nPanel<<" double"<<endl;
                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                        for (iPanel=0; iPanel<nPanelMaster[iProcessor]*nDim; iPanel++){
                                Total_Index  = iProcessor*Max_nPanel*nDim +  iPanel;
                                paraview_file<<std::setprecision(15)<<surfaceCoordsMaster[Total_Index]<<'\t';
                                if (surfaceCoordsMaster[Total_Index] ==0){
                                        cout<<"here we have 0 coords--rank-"<<iProcessor<<" iPanel-"<<iPanel<<" Max_nPanel-"<<Max_nPanel<<" Total Index-"<<Total_Index<<endl;
                                }
                        }
                }
                paraview_file<<endl;
                paraview_file<<"CELLS "<<nFaceTot+nHaloFacesTot<<' '<<nConnSizeTot<<endl;
                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                        for (iConn = 0; iConn < nConnSizeMaster[iProcessor] ; iConn++){
                                Total_Index  = iProcessor*nConnSizeMax+iConn;
                                paraview_file<<Conn_single_master[Total_Index]<<'\t';
                        }
                }

                paraview_file<<endl;
                paraview_file<<"CELL_TYPES "<<nFaceTot+nHaloFacesTot<<endl;
                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                        for (iConn = 0; iConn < nFaceSizeMaster[iProcessor] ; iConn++){
                                paraview_file<<CellTypeMaster[iProcessor*nFaceMax+iConn]<<'\t';
                        }
                }
                paraview_file<<endl;

                if (DERIVATIVES==1){
                        paraview_file<<"POINT_DATA "<<Tot_nPanel<<endl;
                        for (iVar = 0; iVar < nVar+nDim; iVar++){
                                if (iVar==0) paraview_file<<"SCALARS dJ/d(rho) double"<<endl;
                                if (iVar==1) paraview_file<<"SCALARS dJ/d(rhoUx) double"<<endl;
                                if (iVar==2) paraview_file<<"SCALARS dJ/d(rhoUy) double"<<endl;
                                if (iVar==3) paraview_file<<"SCALARS dJ/d(rhoUz) double"<<endl;
                                if (iVar==4) paraview_file<<"SCALARS dJ/d(rhoE) double"<<endl;
                                if (iVar==5) paraview_file<<"SCALARS dJ/d(TKE) double"<<endl;
                                if (iVar==6) paraview_file<<"SCALARS dJ/dX double"<<endl;
                                if (iVar==7) paraview_file<<"SCALARS dJ/dY double"<<endl;
                                if (iVar==8) paraview_file<<"SCALARS dJ/dZ double"<<endl;
                                paraview_file<<"LOOKUP_TABLE default"<<endl;
                                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                                        for (iPanel = 0; iPanel < nPanelMaster[iProcessor]; iPanel++) {
                                                if(iVar<nVar){
                                                        Total_Index  = iProcessor*Max_nPanel*nVar + iVar*nPanelMaster[iProcessor]  + iPanel;
                                                        paraview_file<<std::setprecision(15)<<SingleMaster_dJdU[Total_Index]<< "\t";
                                                }else{
                                                        Total_Index_2= iProcessor*Max_nPanel*nDim + (iVar-nVar)*nPanelMaster[iProcessor]  + iPanel;
                                                        paraview_file<<std::setprecision(15)<< SingleMaster_dJdX[Total_Index_2]   << "\t";
                                                }
                                        }
                                }
                                paraview_file<<endl;
                        }//iVar loop
                }


                if (FIA_TERMS==1){
                        paraview_file<<"POINT_DATA "<<Tot_nPanel<<endl;
                        paraview_file<<"SCALARS pp double"<<endl;
                        paraview_file<<"LOOKUP_TABLE default"<<endl;
                        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                                for (iPanel = 0; iPanel < nPanelMaster[iProcessor]; iPanel++) {
                                        Total_Index  = iProcessor*Max_nPanel +  iPanel;
                                        paraview_file<<std::setprecision(15)<<SingleMaster_pp[Total_Index]<< "\t";
                                }
                        }
                        paraview_file<<endl;
                }

        }//if Master
 //}//iSample loop


}//end of function

void FWHSolver::CombineZones( CConfig *config){

 string data;
 string ttR,ppR,pAcR,TNR,LNR,T1R,T2R,T3R,T4R;
 su2double tt,pAc,TT1,TT2,TT3,TT4;
 su2double pp_mean;
 ifstream ppFWH;
 ofstream CombinedppFWH;
 unsigned long iSample,iObserver,iZ;

 char buffer[32]; // The filename buffer.

 su2double *PP    = new su2double [nSample];
 su2double *t     = new su2double [nSample];
 su2double *Term1 = new su2double [nSample];
 su2double *Term2 = new su2double [nSample];
 su2double *Term3 = new su2double [nSample];
 su2double *Term4 = new su2double [nSample];

 for (iObserver=0; iObserver<nObserver; iObserver++){
	for (iSample=0; iSample<nSample; iSample++){
		PP[iSample]=0.0;t[iSample]=0.0;
		Term1[iSample]=0.0;Term2[iSample]=0.0;Term3[iSample]=0.0;Term4[iSample]=0.0;
	}
	for (iZ=0; iZ<nZone; iZ++){

 		snprintf(buffer, sizeof(char) * 32, "pp_FWH_%03d_Zone_%i", (iObserver+1), (iZ));
 		ppFWH.open(buffer);
        	cout<<"Reading -> "<<buffer<<endl;
        	if (ppFWH.fail()) {
                	cout << "There is no file!!! " <<  buffer  << "."<< endl;
                	exit(EXIT_FAILURE);
        	}
 		getline(ppFWH,data); //to get rid of first line that is headerline
 		ppFWH.close();
 		ppFWH.open(buffer);
 		getline(ppFWH,data); //to get rid of first line that is headerline
		iSample=0;
 		while(getline(ppFWH,data)){
        		istringstream point_data(data);
        		getline(point_data,ttR,' ');
        		getline(point_data,ppR,' ');
        		getline(point_data,pAcR,' ');
        		getline(point_data,TNR,' ');
        		getline(point_data,LNR,' ');
			getline(point_data,T1R,' ');
			getline(point_data,T2R,' ');
			getline(point_data,T3R,' ');
			getline(point_data,T4R,' ');
        		stringstream geek(ttR); geek >> tt;
			if (iZ==0) t[iSample]=tt;
			stringstream geek2(pAcR); geek2 >> pAc; PP[iSample]=PP[iSample]+pAc;
			stringstream geek3(T1R); geek3 >> TT1; Term1[iSample]=Term1[iSample]+TT1;
			stringstream geek4(T2R); geek4 >> TT2; Term2[iSample]=Term2[iSample]+TT2;
			stringstream geek5(T3R); geek5 >> TT3; Term3[iSample]=Term3[iSample]+TT3;
			stringstream geek6(T4R); geek6 >> TT4; Term4[iSample]=Term4[iSample]+TT4;
			iSample++;
 		}
 		ppFWH.close();
	}
	pp_mean=0;
   	for (iSample=0; iSample<nSample; iSample++){
    		pp_mean = (PP[iSample]+pp_mean*iSample)/(iSample+1);
	}
	snprintf(buffer, sizeof(char) * 32, "pp_FWH_%03d_Combined", (iObserver+1));
	CombinedppFWH.open(buffer);
	for (iSample=0; iSample<nSample; iSample++){
            	if (iSample==0){
                	CombinedppFWH  <<  "Time"<<' '<<"P-Fluctuation"<<' '<<"Acoustic-Pressure"<<' '<<"Thickness"<<' '<<"Loading"<<' '<<"Term-1"<<' '<<"Term-2"<<' '<<"Term-3"<<' '<<"Term-4"<<endl;
        	}
         	CombinedppFWH << std::setprecision(15) <<t[iSample] <<' '<<PP[iSample]-pp_mean<<' '<<PP[iSample]<< ' ' <<Term1[iSample]+Term2[iSample]<<' '<< Term3[iSample]+Term4[iSample]<<' '<<Term1[iSample]<<' '<<Term2[iSample]<<' '<<Term3[iSample]<<' '<<Term4[iSample]<<endl;
	}
	CombinedppFWH.close();
 }
	



 delete [] PP;
 delete [] t;
 delete [] Term1;
 delete [] Term2;
 delete [] Term3;
 delete [] Term4;


}
//end of the subroutine: CombineZones



/*--- new implementation for reading ASCII surface tecplot files in v7: remains a work in progress ---*/

F1A::F1A(CConfig *config, CGeometry *geometry) {

    nDim                    = 0;
    nSurfaceNodes           = 0;
    nqSample                = 0;
    SPL                     = 0.0;

    Observer_Locations      = nullptr;
    globalIndexContainer    = nullptr;
    Normal                  = nullptr;
    UnitaryNormal           = nullptr;
    Area                    = nullptr;
    surface_geo             = nullptr;
    FWH_Surf_Vel            = nullptr;
    FWH_Surf_pp             = nullptr;
    Momentum                = nullptr;
    Q                       = nullptr;
    F                       = nullptr;
    RHO                     = nullptr;
    RadVec                  = nullptr;
    StartTime               = nullptr;
    EndTime                 = nullptr;
    dt                      = nullptr;
    pp_out                  = nullptr;
    pp_fft                  = nullptr;
      
    /*--- Store MPI rank and size ---*/

    rank = SU2_MPI::GetRank();
    size = SU2_MPI::GetSize();

}

F1A::~F1A(void) {

  delete [] Observer_Locations;
  delete [] Normal;
  delete [] UnitaryNormal;
  delete [] Area;
  delete [] globalIndexContainer;
  delete [] surface_geo;
  delete [] FWH_Surf_Vel;
  delete [] FWH_Surf_pp;
  delete [] Momentum;
  delete [] Q;
  delete [] F;
  delete [] RHO;
  delete [] RadVec;
  delete [] StartTime;
  delete [] EndTime;
  delete [] dt;
  delete [] pp_out;
  delete [] pp_fft;
}

void F1A::Initialize(CConfig *config, CGeometry *geometry) {

    unsigned long iVertex, iPoint, iMarker, Global_Index;
    unsigned long index = 0;

    unsigned long localnSurfaceNodes = 0;

    nDim                = geometry->GetnDim();
    nqSample            = config->GetAcoustic_nqSamples();
    SamplingFreq        = config->GetWrt_Sol_Freq_DualTime();
    nSample             = config->GetIter_Avg_Objective()/SamplingFreq;

    /*--- Get the number of surface nodes ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        /* --- Loop over boundary markers to select those on the FWH surface --- */
        if (config->GetMarker_All_KindBC(iMarker) == ACOUSTIC_BOUNDARY) {
            /*--- Get the total number of FWH surface nodes ---*/
            //nSurfaceNodes += geometry->nVertex[iMarker]; //does not work in parallel

            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
                if (geometry->nodes->GetDomain(iPoint)) {
                    localnSurfaceNodes++;
                }
            }

        }
    }

    SU2_MPI::Allreduce(&localnSurfaceNodes, &nSurfaceNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    cout << "rank " << rank << ", localnSurfaceNodes " << localnSurfaceNodes << ", nSurfaceNodes " << nSurfaceNodes << endl;

    globalIndexContainer = new unsigned long[nSurfaceNodes];

    /*--- Build an array containing the global indices for the surface nodes ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        /* --- Loop over boundary markers to select those on the FWH surface --- */
        if (config->GetMarker_All_KindBC(iMarker) == ACOUSTIC_BOUNDARY) {

            /*--- Loop over all the vertices on this boundary marker ---*/
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

                /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
                if (geometry->nodes->GetDomain(iPoint)) {

                    Global_Index = geometry->nodes->GetGlobalIndex(iPoint);
                    globalIndexContainer[index] = Global_Index;
                    index++;

                }
            }
        }
    }

    /*--- Sort the global index container in ascending order ---*/
    unsigned long tmp = 0;
    for (unsigned long i = 0; i < nSurfaceNodes; i++) {
        for (unsigned long j = i + 1; j < nSurfaceNodes; j++)
        {
            if (globalIndexContainer[j] < globalIndexContainer[i])
            {
                tmp = globalIndexContainer[i];
                globalIndexContainer[i] = globalIndexContainer[j];
                globalIndexContainer[j] = tmp;
            }
        }
    }

    nPanel      = nSurfaceNodes;
    surface_geo = new su2double [nSurfaceNodes*nqSample*nDim];
    FWH_Surf_Vel= new su2double [nSurfaceNodes*nqSample*nDim]; 
    Momentum    = new su2double [nSurfaceNodes*nqSample*nDim];
    Q           = new su2double [nSurfaceNodes*nqSample*nDim];
    F           = new su2double [nSurfaceNodes*nqSample*nDim];       
    RHO         = new su2double [nSurfaceNodes*nqSample];
    FWH_Surf_pp = new su2double [nSurfaceNodes*nqSample];
    Normal          = new su2double[nqSample*nSurfaceNodes*nDim];
    UnitaryNormal   = new su2double[nqSample*nSurfaceNodes*nDim];
    Area            = new su2double[nqSample*nSurfaceNodes];

    RadVec      = new su2double [15];
    for (int i=0; i<15; i++) {
        RadVec[i]=0.0;
    }

    //add this as a config option later
    string  text_line;
    ifstream Observer_LocationFile;
    string filename = "Observer_Locations.dat";
    Observer_LocationFile.open(filename.data() , ios::in);
    if (Observer_LocationFile.fail()) {
        cout << "There is no file!!! " <<  filename.data()  << "."<< endl;
        exit(EXIT_FAILURE);
    }
    getline (Observer_LocationFile, text_line);
    istringstream point_line(text_line);
    point_line >> nObserver ;
    Observer_Locations = new su2double [nObserver*nDim];

    unsigned long iObserver=0;
    while (getline(Observer_LocationFile, text_line)){
        istringstream point_line2(text_line);
        point_line2 >> Observer_Locations[iObserver*nDim + 0] >> Observer_Locations[iObserver*nDim + 1] >> Observer_Locations[iObserver*nDim + 2];
        iObserver++;
    }

    StartTime   = new su2double [nObserver];
    EndTime     = new su2double [nObserver];
    dt          = new su2double [nObserver];

    for(iObserver = 0; iObserver<nObserver ; iObserver++){
        StartTime[iObserver]= 0.0;
        EndTime[iObserver]= 9999999999.0;
        dt[iObserver]= 0.0;
    }

    su2double R = 287.058;
    FreeStreamPressure  =   config->GetPressure_FreeStream();
    FreeStreamDensity   =   FreeStreamPressure/R/config->GetTemperature_FreeStream();;

    su2double M     = config->GetMach();
    su2double AOA   = config->GetAoA()*PI_NUMBER/180.0;
    su2double AOS   = config->GetAoS()*PI_NUMBER/180.0;

    a_inf   = sqrt(config->GetGamma()*FreeStreamPressure / FreeStreamDensity);
    U1      = M*a_inf*cos(AOA)*cos(AOS) ;
    U2      = M*a_inf*sin(AOS) ;
    U3      = M*a_inf*sin(AOA)*cos(AOS);

    pp_fft  = new complex <su2double> [nSample];
    for (unsigned long iSample=0; iSample < nSample; iSample++){
        pp_fft[iSample]= 0.0;
    }

    pp_out = new su2double [nPanel];
    for (unsigned long iPanel=0; iPanel<nPanel; iPanel++){
        pp_out[iPanel]=0.0;
    }

}

void F1A::Read_TECPLOT_ASCII( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample){

#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == MASTER_NODE) {

        ifstream surface_file;
        string text_line, useless, Tag;
        string filename = "surface_flow";
        string ext = ".dat";

        int iVar;
        int nFields;
        int x_index = 0, y_index = 1, z_index = 2;
        int density_index = 3;
        int momentum_x_index = 4, momentum_y_index = 5, momentum_z_index = 6;
        int grid_vel_x_index = 9, grid_vel_y_index = 10, grid_vel_z_index = 11;
        int pressure_index = 12;
        int normal_x_index = 13, normal_y_index = 14, normal_z_index = 15;

        bool moving_surface = false;
        bool normals_output = false;

        vector<string> fields;
        fields.clear();

        passivedouble *Aux_DataStruct = nullptr;

        unsigned long ExtIter = (iSample + iLocSample) * SamplingFreq + config->GetRestart_Iter();

        filename = config->GetFilename(filename, ext, ExtIter);

        /*--- Open the restart file ---*/

        surface_file.open(filename.data(), ios::in);

        /*--- In case there is no restart file ---*/

        if (surface_file.fail()) {
            cout << ("TECPLOT ASCII surface file ") + string(filename) + string(" not found.\n");
            exit(EXIT_FAILURE);
        } else {
            cout << "Reading Tecplot surface file " << filename << endl;
        }

        /*--- Identify the number of fields (and names) in the restart file ---*/

        getline(surface_file, useless);
        getline(surface_file, text_line);

        char delimiter = ',';
        fields = PrintingToolbox::split(text_line, delimiter);

        //at index 0 trim tecplot field: VARIABLES =
        fields[0].erase(0, 12);

        for (auto &field : fields) {
            PrintingToolbox::trim(field);
            field.erase(remove(field.begin(), field.end(), '"'), field.end());
        }

        for (int iField = 0; iField < fields.size(); iField++) {
            if (fields[iField] == "x") x_index = iField;
            if (fields[iField] == "y") y_index = iField;
            if (fields[iField] == "z") z_index = iField;
            if (fields[iField] == "Density") density_index = iField;
            if (fields[iField] == "Pressure") pressure_index = iField;
            if (fields[iField] == "Momentum_x") momentum_x_index = iField;
            if (fields[iField] == "Momentum_y") momentum_y_index = iField;
            if (fields[iField] == "Momentum_z") momentum_z_index = iField;

            if (fields[iField] == "Grid_Velocity_x") {
                grid_vel_x_index = iField;
                moving_surface = true;
            }
            if (fields[iField] == "Grid_Velocity_y") grid_vel_y_index = iField;
            if (fields[iField] == "Grid_Velocity_z") grid_vel_z_index = iField;

            if (fields[iField] == "Normal_x") {
                normal_x_index = iField;
                normals_output = true;
            }
            if (fields[iField] == "Normal_y") normal_y_index = iField;
            if (fields[iField] == "Normal_z") normal_z_index = iField;
        }

        /*--- Set the number of variables, one per field in the restart file  ---*/

        nFields = (int) fields.size();

        /*--- Allocate memory for the restart data. ---*/

        Aux_DataStruct = new passivedouble[nFields * nSurfaceNodes];

        //skip line
        getline(surface_file, useless);

        /*--- Read all lines in the restart file and extract data. ---*/

        for (int iNode = 0; iNode < nSurfaceNodes; iNode++) {

            getline(surface_file, text_line);

            delimiter = '\t';
            vector<string> point_line = PrintingToolbox::split(text_line, delimiter);

            for (iVar = 0; iVar < nFields; iVar++) {
                Aux_DataStruct[iNode * nFields + iVar] = SU2_TYPE::GetValue(PrintingToolbox::stod(point_line[iVar]));
            }
        }

        //assign the auxillary data to the F1A global variables
        for (int iNode = 0; iNode < nSurfaceNodes; iNode++) {

            surface_geo[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 0] = Aux_DataStruct[iNode * nFields + x_index];
            surface_geo[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 1] = Aux_DataStruct[iNode * nFields + y_index];
            surface_geo[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 2] = Aux_DataStruct[iNode * nFields + z_index];

            RHO[iLocSample * nSurfaceNodes + iNode] = Aux_DataStruct[iNode * nFields + density_index];
            FWH_Surf_pp[iLocSample * nSurfaceNodes + iNode] =
                    Aux_DataStruct[iNode * nFields + pressure_index] - FreeStreamPressure;

            if (moving_surface) {
                FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 0] = Aux_DataStruct[iNode * nFields + grid_vel_x_index];
                FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 1] = Aux_DataStruct[iNode * nFields + grid_vel_y_index];
                FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 2] = Aux_DataStruct[iNode * nFields + grid_vel_z_index];
            } else {
                FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 0] = 0.0;
                FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 1] = 0.0;
                FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 2] = 0.0;
            }

            Momentum[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 0] = Aux_DataStruct[iNode * nFields + momentum_x_index];
            Momentum[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 1] = Aux_DataStruct[iNode * nFields + momentum_y_index];
            Momentum[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 2] = Aux_DataStruct[iNode * nFields + momentum_z_index];

            if (normals_output) {
                Normal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 0] = Aux_DataStruct[iNode * nFields + normal_x_index];
                Normal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 1] = Aux_DataStruct[iNode * nFields + normal_y_index];
                Normal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 2] = Aux_DataStruct[iNode * nFields + normal_z_index];
            } else {
                Normal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 0] = 0.0;
                Normal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 1] = 0.0;
                Normal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 2] = 0.0;
            }
        }

        surface_file.close();

        delete[] Aux_DataStruct;

    }

    SU2_MPI::Bcast(&surface_geo[iLocSample * nSurfaceNodes * nDim], nSurfaceNodes*nDim, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&Normal[iLocSample * nSurfaceNodes * nDim], nSurfaceNodes*nDim, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&RHO[iLocSample * nSurfaceNodes], nSurfaceNodes, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&Momentum[iLocSample * nSurfaceNodes * nDim], nSurfaceNodes*nDim, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&FWH_Surf_pp[iLocSample * nSurfaceNodes], nSurfaceNodes, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim], nSurfaceNodes*nDim, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

}

void F1A::LoadNormal(CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long iObserver) {

    ifstream surface_file;
    string text_line, useless, Tag;
    string filename = "surface_normals";
    string ext = ".dat";

    int iVar;
    int nFields;

    passivedouble *Aux_DataStruct = nullptr;

    unsigned long ExtIter = config->GetRestart_Iter();

    filename = config->GetFilename(filename, ext, ExtIter);

    /*--- Open the restart file ---*/

    surface_file.open(filename.data(), ios::in);

    /*--- In case there is no restart file ---*/

    if (surface_file.fail()) {
        cout << ("TECPLOT ASCII surface file ") + string(filename) + string(" not found.\n") ;
        exit(EXIT_FAILURE);
    }
    else{
        cout << "Reading Tecplot surface file " << filename << endl;
    }

    /*--- Skip the unnecessary data at the beginning ---*/
    for (int i = 0; i < 13; i++){
        getline (surface_file, useless);
    }

    nFields = 8;
    Aux_DataStruct = new passivedouble[nFields * nSurfaceNodes];

    char delimiter = ' ';

    /*--- Load the data ---*/
    for (int iNode = 0; iNode < nSurfaceNodes; iNode++ ) {

        getline(surface_file, text_line);

        vector<string> point_line = PrintingToolbox::split(text_line, delimiter);

        for (iVar = 0; iVar < nFields; iVar++) {
            Aux_DataStruct[iNode * nFields + iVar] = SU2_TYPE::GetValue(PrintingToolbox::stod(point_line[iVar]));
        }
    }

    for (int i = 0; i < nqSample; i++) {
        for (int iNode = 0; iNode < nSurfaceNodes; iNode++) {

            UnitaryNormal[i * nSurfaceNodes * nDim + iNode * nDim + 0] = Aux_DataStruct[iNode * nFields + 4];
            UnitaryNormal[i * nSurfaceNodes * nDim + iNode * nDim + 1] = Aux_DataStruct[iNode * nFields + 5];
            UnitaryNormal[i * nSurfaceNodes * nDim + iNode * nDim + 2] = Aux_DataStruct[iNode * nFields + 6];

            Area[i * nSurfaceNodes + iNode] = Aux_DataStruct[iNode * nFields + 7];
        }
    }

    surface_file.close();

    delete [] Aux_DataStruct;

}


void F1A::ComputeNormal( CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long iObserver){
/*
    unsigned long iVertex, iPoint, iMarker;
    unsigned short iDim;

    su2double *normal = nullptr;
    normal = new su2double[nDim];

    for (iPoint = 0; iPoint < nSurfaceNodes; iPoint++) {

        //allocate memory to the the heap to access across routines!
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
            // --- Loop over boundary markers to select those on the FWH surface --- //
            if (config->GetMarker_All_KindBC(iMarker) == ACOUSTIC_BOUNDARY) {

                iVertex = geometry->nodes->GetVertex(globalIndexContainer[iPoint], iMarker);

                //--- Check the vertex is contained by the marker ---//
                if (iVertex != -1) {

                    //--- Normal vector for this vertex (negate for outward convention) ---//
                    geometry->vertex[iMarker][iVertex]->GetNormal(normal);

                    for (iDim = 0; iDim < nDim; iDim++) {
                        //Normal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim] = -normal[iDim];
                        Normal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim] = normal[iDim];
                    }

                    Area[iLocSample * nSurfaceNodes + iPoint] = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Area[iLocSample * nSurfaceNodes + iPoint] += pow(Normal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim],2);
                    }

                    Area[iLocSample * nSurfaceNodes + iPoint] = sqrt(Area[iLocSample * nSurfaceNodes + iPoint]);

                    for (iDim = 0; iDim < nDim; iDim++) {
                        UnitaryNormal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim] =
                                Normal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim] / Area[iLocSample * nSurfaceNodes + iPoint];
                    }

                }
            }
        }
    }

    delete [] normal;
*/

    unsigned long iPoint;
    unsigned short iDim;

    for (iPoint = 0; iPoint < nSurfaceNodes; iPoint++) {
        Area[iLocSample * nSurfaceNodes + iPoint] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
            Area[iLocSample * nSurfaceNodes + iPoint] += pow( Normal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim], 2);
        }

        Area[iLocSample * nSurfaceNodes + iPoint] = sqrt(Area[iLocSample * nSurfaceNodes + iPoint]);

        for (iDim = 0; iDim < nDim; iDim++) {
            UnitaryNormal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim] =
                    Normal[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim] /
                    Area[iLocSample * nSurfaceNodes + iPoint];
        }
    }

}

void F1A::SetSurfaceGeom(CConfig *config, CGeometry *geometry, unsigned long iSample, unsigned long iLocSample, unsigned long iObserver) {

    unsigned long  iPoint;
    unsigned short iDim;

    for (iPoint = 0; iPoint < nSurfaceNodes; iPoint++) {
        su2double Coord[3];

        for (iDim = 0; iDim < nDim; iDim++) {
            Coord[iDim] = surface_geo[iLocSample * nSurfaceNodes * nDim + iPoint * nDim + iDim];
        }

        for (iDim = 0; iDim < nDim; iDim++) {
            geometry->nodes->SetCoord(globalIndexContainer[iPoint], iDim, Coord[iDim]);
        }
    }

    UpdateDualGrid(geometry, config);
}

void F1A::UpdateDualGrid(CGeometry *geometry, CConfig *config){

    /*--- After moving all nodes, update the dual mesh. Recompute the edges and
     dual mesh control volumes in the domain and on the boundaries. ---*/

    geometry->SetCoord_CG();
    geometry->SetControlVolume(config, UPDATE);
    geometry->SetBoundControlVolume(config, UPDATE);
    geometry->SetMaxLength(config);

}

void F1A::ComputeModifiedVelocity(CConfig *config, unsigned long iSample, unsigned long iLocSample) {

    su2double *u  = nullptr;
    su2double *v  = nullptr;
    su2double *Li = nullptr;
    su2double *U  = nullptr;

    u  = new su2double [nDim];
    v  = new su2double [nDim];
    Li = new su2double [nDim];
    U  = new su2double [nDim];

    U[0] = U1;
    U[1] = U2;
    U[2] = U3;

    // Compute the modified grid velocity and modified stress tensor
    for (int iNode = 0; iNode < nSurfaceNodes; iNode++ ) {
        for (int iDim=0; iDim<nDim; iDim++){
            v[iDim] = FWH_Surf_Vel[iLocSample * nSurfaceNodes * nDim + iNode * nDim + iDim] - U[iDim];
            u[iDim] = Momentum[iLocSample * nSurfaceNodes * nDim + iNode * nDim + iDim] / RHO[iLocSample * nSurfaceNodes + iNode];

            //modified grid velocity
            Q[iLocSample * nSurfaceNodes * nDim + iNode * nDim + iDim] = u[iDim] + ( (RHO[iLocSample * nSurfaceNodes + iNode] / FreeStreamDensity) - 1.0 ) * (u[iDim] - v[iDim]);
        }

        // Lij modified stress tensor
        for (int i=0; i<nDim; i++){
            for (int j=0; j<nDim; j++){
                Li[j] = Momentum[iLocSample * nSurfaceNodes * nDim + iNode * nDim + i]*(u[j]-v[j]);
                if (i==j) Li[j] += FWH_Surf_pp[iLocSample * nSurfaceNodes + iNode];
            }
            // F = Lij * nj
            F[iLocSample * nSurfaceNodes * nDim + iNode * nDim + i] =
                      Li[0] * UnitaryNormal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 0]
                    + Li[1] * UnitaryNormal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 1]
                    + Li[2] * UnitaryNormal[iLocSample * nSurfaceNodes * nDim + iNode * nDim + 2];
        }
    }

    delete [] u;
    delete [] v;
    delete [] Li;
    delete [] U;

}

void F1A::ComputeMinMaxInc_Time(CConfig *config, CGeometry *geometry) {

    unsigned long iLocSample, iObserver, iPanel;

#ifdef HAVE_MPI
    int rank, nProcessor;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

    su2double tau;
    iLocSample = 0;

    Read_TECPLOT_ASCII(config, geometry, 0, iLocSample);
    Read_TECPLOT_ASCII(config, geometry, nSample - 1, iLocSample + 1);

    //SetSurfaceGeom(config, geometry, 0, iLocSample, 0);
    ComputeNormal(config, geometry, 0, iLocSample, 0);

    for (iObserver=0; iObserver<nObserver; iObserver++){
        for (iLocSample=0; iLocSample<2; iLocSample++){
            for (iPanel=0; iPanel<nPanel; iPanel++){
                ComputeObserverTime( config, iObserver, iPanel, iLocSample*(nSample-1), iLocSample,0);
                tau= RadVec[4];
                if (iLocSample==0){
                    if (tau>StartTime[iObserver]){
                        StartTime[iObserver]=tau;
                    }
                }else{
                    if (tau<EndTime[iObserver]){
                        EndTime[iObserver]=tau;
                    }
                }
            }
        }
#ifdef HAVE_MPI
        if (EndTime[iObserver]==0.0){
                        	EndTime[iObserver]=999999999;
                	}
			SU2_MPI::Allreduce(&StartTime[iObserver],&StartTime[iObserver],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			SU2_MPI::Allreduce(&EndTime[iObserver],&EndTime[iObserver],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
        dt[iObserver]=(EndTime[iObserver]-StartTime[iObserver])/nSample;
        StartTime[iObserver]=StartTime[iObserver]+5*dt[iObserver];
        EndTime[iObserver]=EndTime[iObserver]-5*dt[iObserver];
        dt[iObserver]=(EndTime[iObserver]-StartTime[iObserver])/nSample;
        if (rank==MASTER_NODE) cout<<"For Observer-"<<iObserver<<", Endtime="<<EndTime[iObserver]<<" and Startime="<<StartTime[iObserver]<<endl;
        if (rank==MASTER_NODE) cout<<"dt="<<dt[iObserver]<<endl;
    }
}

void F1A::ComputeObserverTime( CConfig *config, unsigned long iObserver, unsigned long iPanel,unsigned long iSample, unsigned long iLocSample, unsigned long i){

/*This subroutine computes observer time and radiation vector for each panel */
    unsigned long iDim;
    su2double x,y,z;
    su2double xp[3]= {0.0,0.0,0.0};
    su2double xo[3]= {0.0,0.0,0.0};
    su2double U[3] = {U1, U2, U3};
    su2double dtt=0.1;
    su2double eps=1.0;
    su2double dtnew,diff,R,r1,r2,r3,r_mag,Time;

    unsigned long start_iter = config->GetRestart_Iter();

    Time= config->GetDelta_UnstTime()*(start_iter+(iSample+iLocSample)*SamplingFreq);

    x=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 0];
    y=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 1];
    z=surface_geo[iLocSample*nPanel*nDim + iPanel*nDim + 2];

    x=x-U1*Time;
    y=y-U2*Time;
    z=z-U3*Time;

    for (iDim=0; iDim<nDim; iDim++){
        xo[iDim]=Observer_Locations[iObserver*nDim + iDim]-U[iDim]*Time;
        xp[iDim]=xo[iDim]-U[iDim]*dtt;
    }
    while (eps>0.0000000001){
        R=sqrt((xp[0]-x)*(xp[0]-x)+(xp[1]-y)*(xp[1]-y)+(xp[2]-z)*(xp[2]-z));
        dtnew=R/a_inf;
        diff=dtnew-dtt;
        dtt=dtnew;
        eps=diff/dtnew*10000;
        eps=sqrt(eps*eps);
        for (iDim=0; iDim<nDim; iDim++) xp[iDim]=xo[iDim]-U[iDim]*dtt;
    }
    r1 = xp[0]-x;
    r2 = xp[1]-y;
    r3 = xp[2]-z;
    r_mag = sqrt(r1*r1+r2*r2+r3*r3);
    r1 = r1/r_mag; r2 = r2/r_mag;r3 = r3/r_mag;

    RadVec[5*i+0]=r1;
    RadVec[5*i+1]=r2;
    RadVec[5*i+2]=r3;
    RadVec[5*i+3]=r_mag;
    RadVec[5*i+4]=Time+dtt;//When the observer hears the signal

}


void F1A::F1A_SourceTimeDominant( CConfig *config, CGeometry *geometry){

#ifdef HAVE_MPI
    int rank, nProcessor;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
    su2double PP[3]= {0.0,0.0,0.0};
    su2double TT1[3]={0.0,0.0,0.0};
    su2double TT2[3]={0.0,0.0,0.0};
    su2double TT3[3]={0.0,0.0,0.0};
    su2double TT4[3]={0.0,0.0,0.0};
    su2double tau[3]= {0.0,0.0,0.0};
    unsigned long iObserver, iSample, iSample2, iDim, iLocSample, iPanel,iqSample, Smallest_iLocSample;
    su2double delta, dtt;
    su2double pp_n, TT1_n,TT2_n,TT3_n,TT4_n;
    su2double r1,r2,r3,R;
    su2double pp_mean, SPL_iObserver;
    char buffer[32];
    unsigned long res_factor=1;
    su2double dt_src,t;
    unsigned long iObserver_inner, iObserver_outer;
    unsigned long nObserver_inner, nObserver_outer;
    unsigned long Pre_Mot=0;
    dt_src=config->GetDelta_UnstTime()*SamplingFreq;
    ofstream pp_FWH_file ;
    ofstream pp_fft_file ;
    if (ceil(log2(nSample))!=floor(log2(nSample))){
        if (rank==MASTER_NODE){
            cout<<"Number of samples are not power of 2 ---> no FFT computation "<<endl;
        }
    }
    su2double dTT = config->GetDelta_UnstTime()*SamplingFreq; //Time Step to calculate sampling frequency
    su2double Fs = 1/dTT; //Sampling frequency to be used in fft calculation

    su2double **t_interp_inner= NULL,**pp_TimeDomain_inner = NULL,**pp_TimeDomainGlobal_inner = NULL;
    su2double **T1_I=NULL,**T2_I=NULL,**T3_I=NULL,**T4_I=NULL;
    su2double *t_interp_outer = NULL,*pp_TimeDomain_outer = NULL,*pp_TimeDomainGlobal_outer = NULL;
    su2double *T1_O=NULL,*T2_O=NULL,*T3_O=NULL,*T4_O=NULL;


    if(config->GetAcoustic_Inner_ObsLoop()){

        delete[] t_interp_outer;
        delete[] pp_TimeDomain_outer;
        delete[] pp_TimeDomainGlobal_outer;
        t_interp_inner= new su2double *[nObserver];
        pp_TimeDomain_inner=new su2double *[nObserver];
        if (rank==MASTER_NODE) pp_TimeDomainGlobal_inner=new su2double *[nObserver];
        T1_I=new su2double *[nObserver];
        T2_I=new su2double *[nObserver];
        T3_I=new su2double *[nObserver];
        T4_I=new su2double *[nObserver];
        for (iObserver=0; iObserver<nObserver; iObserver++){
            t_interp_inner[iObserver]=new su2double [nSample];
            pp_TimeDomain_inner[iObserver]=new su2double [nSample];
            if (rank==MASTER_NODE) pp_TimeDomainGlobal_inner[iObserver]=new su2double [nSample];
            T1_I[iObserver]=new su2double [nSample];
            T2_I[iObserver]=new su2double [nSample];
            T3_I[iObserver]=new su2double [nSample];
            T4_I[iObserver]=new su2double [nSample];
            for (iSample=0; iSample<nSample; iSample++){
                t_interp_inner[iObserver][iSample]=0.0;
                pp_TimeDomain_inner[iObserver][iSample]=0.0;
                if (rank==MASTER_NODE) pp_TimeDomainGlobal_inner[iObserver][iSample]=0.0;
                T1_I[iObserver][iSample]=0.0;
                T2_I[iObserver][iSample]=0.0;
                T3_I[iObserver][iSample]=0.0;
                T4_I[iObserver][iSample]=0.0;
            }
        }

    }else{
        delete[] t_interp_inner;
        delete[] pp_TimeDomain_inner;
        delete[] pp_TimeDomainGlobal_inner;
        t_interp_outer= new su2double [nSample];
        pp_TimeDomain_outer=new su2double [nSample];
        if (rank==MASTER_NODE) pp_TimeDomainGlobal_outer=new su2double [nSample];
        T1_O=new su2double [nSample];
        T2_O=new su2double [nSample];
        T3_O=new su2double [nSample];
        T4_O=new su2double [nSample];
    }

    if(config->GetAcoustic_Inner_ObsLoop()){
        for (iObserver=0; iObserver<nObserver; iObserver++){
            t_interp_inner[iObserver][0]=StartTime[iObserver];
        }
    }

    if(config->GetAcoustic_Inner_ObsLoop()){
        nObserver_inner=nObserver;
        nObserver_outer=1;
    }else{
        nObserver_outer=nObserver;
        nObserver_inner=1;
    }

    for (iObserver_outer=0; iObserver_outer<nObserver_outer; iObserver_outer++){
        if(!(config->GetAcoustic_Inner_ObsLoop())){
            for (iSample=0; iSample<nSample; iSample++){
                t_interp_outer[iSample]=0.0;
                pp_TimeDomain_outer[iSample]=0.0;
                if (rank==MASTER_NODE) pp_TimeDomainGlobal_outer[iSample]=0.0;
                T1_O[iSample]=0.0;
                T2_O[iSample]=0.0;
                T3_O[iSample]=0.0;
                T4_O[iSample]=0.0;
            }
            t_interp_outer[0]=StartTime[iObserver_outer];
        }

        iSample=0;
        for (iLocSample=0; iLocSample<nqSample; iLocSample++){
            for (iPanel=0; iPanel < nPanel; iPanel++){
                for (iDim=0; iDim < nDim; iDim++){
                    UnitaryNormal[iLocSample*nPanel*nDim + iPanel*nDim + iDim ]=0.0;
                }
            }

            Read_TECPLOT_ASCII( config, geometry, iSample, iLocSample);
            //SetSurfaceGeom(config, geometry, iSample, iLocSample, iObserver_outer);
            ComputeNormal(config, geometry, iSample, iLocSample, iObserver_outer);
            ComputeModifiedVelocity(config, iSample, iLocSample);

        }
        iLocSample=nqSample/2;
        iSample2=0;
        for (iSample=0; iSample<nSample; iSample++){
            if (rank==MASTER_NODE) cout<<"Computing the Acoustics at iSample-"<<iSample<<endl;
            Smallest_iLocSample=999999;
            for (iPanel=0; iPanel<nPanel; iPanel++){
                for (iObserver_inner=0; iObserver_inner<nObserver_inner; iObserver_inner++){
                    if(config->GetAcoustic_Inner_ObsLoop()){
                        iObserver=iObserver_inner;
                        t=t_interp_inner[iObserver][iSample];
                    }else{
                        iObserver=iObserver_outer;
                        t=t_interp_outer[iSample];
                    }
                    for (int i=0; i<3 ; i++){
                        //cout<<"Hey-15 iSample-"<<iSample<<" i: "<<i<<" iLocSample: "<<iLocSample<<endl;
                        ComputeObserverTime(config, iObserver, iPanel, iSample2, iLocSample-1+i,i);
                        tau[i] = RadVec[5*i+4];
                    }
                    //cout<<"Hey-17 iSample-"<<iSample<<" iPanel="<<iPanel<<" t_interp:"<<t<<" taus: "<<tau[0]<<" "<<tau[1]<<" "<<tau[2]<<" iLocSample-"<<iLocSample<<endl;
                    while ( (t<tau[1]) || (t>tau[2]) ){
                        if ( t<tau[1] ){
                            if ( (iLocSample<3)  ) goto endloop;
                            /*update tau and shift down*/
                            iLocSample=iLocSample-1;
                        }
                        if ( t>tau[2] ){
                            if ( (iLocSample>nqSample-4) ) goto endloop;
                            iLocSample=iLocSample+1;
                        }

                        for (int i=0; i<3 ; i++){
                            ComputeObserverTime(config, iObserver, iPanel, iSample2, iLocSample-1+i,i);
                            tau[i] = RadVec[5*i+4];
                        }
                    }
                    endloop:
                    if (iLocSample<Smallest_iLocSample){
                        Smallest_iLocSample=iLocSample;
                    }
                    if (iPanel==0){
                        if (iSample<nSample-1){
                            if(config->GetAcoustic_Inner_ObsLoop()){
                                t_interp_inner[iObserver][iSample+1]=t_interp_inner[iObserver][iSample]+dt[iObserver];
                            }else{
                                t_interp_outer[iSample+1]=t_interp_outer[iSample]+dt[iObserver];
                            }
                        }
                    }
                    for (int i=0; i<3 ; i++){
                        F1A_Formulation ( config,iObserver, iPanel, iSample,  iLocSample-1+i, i);
                        PP[i]=pp_t;
                        TT1[i]=T1;TT2[i]=T2;TT3[i]=T3;TT4[i]=T4+T5;
                    }
                    delta=t-tau[1];
                    dtt=tau[2]-tau[1];
                    /*2nd order taylor expansion for polynamial interpolation*/
                    pp_n=PP[1]+delta*(PP[2]-PP[0])*0.5/dtt+delta*delta*0.5*(PP[2]-2*PP[1]+PP[0])/dtt/dtt;
                    TT1_n=TT1[1]+delta*(TT1[2]-TT1[0])*0.5/dtt+delta*delta*0.5*(TT1[2]-2*TT1[1]+TT1[0])/dtt/dtt;
                    TT2_n=TT2[1]+delta*(TT2[2]-TT2[0])*0.5/dtt+delta*delta*0.5*(TT2[2]-2*TT2[1]+TT2[0])/dtt/dtt;
                    TT3_n=TT3[1]+delta*(TT3[2]-TT3[0])*0.5/dtt+delta*delta*0.5*(TT3[2]-2*TT3[1]+TT3[0])/dtt/dtt;
                    TT4_n=TT4[1]+delta*(TT4[2]-TT4[0])*0.5/dtt+delta*delta*0.5*(TT4[2]-2*TT4[1]+TT4[0])/dtt/dtt;
                    if(config->GetAcoustic_Inner_ObsLoop()){
                        pp_TimeDomain_inner[iObserver][iSample] = pp_TimeDomain_inner[iObserver][iSample] + pp_n;
                        T1_I[iObserver][iSample]=T1_I[iObserver][iSample]+TT1_n;
                        T2_I[iObserver][iSample]=T2_I[iObserver][iSample]+TT2_n;
                        T3_I[iObserver][iSample]=T3_I[iObserver][iSample]+TT3_n;
                        T4_I[iObserver][iSample]=T4_I[iObserver][iSample]+TT4_n;
                    }else{
                        pp_TimeDomain_outer[iSample] = pp_TimeDomain_outer[iSample] + pp_n;
                        T1_O[iSample]=T1_O[iSample]+TT1_n;
                        T2_O[iSample]=T2_O[iSample]+TT2_n;
                        T3_O[iSample]=T3_O[iSample]+TT3_n;
                        T4_O[iSample]=T4_O[iSample]+TT4_n;
                    }

                    //if ((config->GetOutput_FileFormat() == PARAVIEW)) {
                    pp_out[iPanel]=pp_n;
                    R=RadVec[5*1+3]+delta*(RadVec[5*2+3]-RadVec[5*0+3])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+3]-2*RadVec[5*1+3]+RadVec[5*0+3])/dtt/dtt;
                    r1=RadVec[5*1+0]+delta*(RadVec[5*2+0]-RadVec[5*0+0])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+0]-2*RadVec[5*1+0]+RadVec[5*0+0])/dtt/dtt;
                    r2=RadVec[5*1+1]+delta*(RadVec[5*2+1]-RadVec[5*0+1])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+1]-2*RadVec[5*1+1]+RadVec[5*0+1])/dtt/dtt;
                    r3=RadVec[5*1+2]+delta*(RadVec[5*2+2]-RadVec[5*0+2])*0.5/dtt+delta*delta*0.5*(RadVec[5*2+2]-2*RadVec[5*1+2]+RadVec[5*0+2])/dtt/dtt;
                    //RetSurf[iPanel*nDim+0]=-(r1*R)+Observer_Locations[iObserver][0]-U1*t;
                    //RetSurf[iPanel*nDim+1]=-(r2*R)+Observer_Locations[iObserver][1]-U2*t;
                    //RetSurf[iPanel*nDim+2]=-(r3*R)+Observer_Locations[iObserver][2]-U3*t;
                    //}

                }//iObserver_inner Loop
            }//iPanel Loop
            //cout<<"rank"<<rank<<"hey-25"<<endl;
            //PARAVIEW Section for Output

            //if ((config->GetOutput_FileFormat() == PARAVIEW)) {
            //if (rank==MASTER_NODE)  cout<<"Extracting paraview file (.vtk) for post-processing"<<endl;
            //Paraview_Output(config,iObserver,iSample, 0, 0, 1);
            //if (rank==MASTER_NODE)  cout<<"End of paraview extraction"<<endl;
            //}//if Paraview

#ifdef HAVE_MPI
            SU2_MPI::Allreduce(&Smallest_iLocSample,&Smallest_iLocSample,1,MPI_UNSIGNED_LONG,MPI_MIN,MPI_COMM_WORLD);
#endif
            if ( (Smallest_iLocSample>5) ){
                for (iqSample=0; iqSample<nqSample-1; iqSample++){
                    for (iPanel=0; iPanel<nPanel; iPanel++){
                        for (iDim=0; iDim<nDim; iDim++){

                            Q[iqSample*nPanel*nDim + iPanel*nDim + iDim] = Q[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];
                            UnitaryNormal[iqSample*nPanel*nDim + iPanel*nDim + iDim]=UnitaryNormal[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];
                            surface_geo[iqSample*nPanel*nDim + iPanel*nDim + iDim]=surface_geo[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];
                            F[iqSample*nPanel*nDim + iPanel*nDim + iDim]=F[(iqSample+1)*nPanel*nDim + iPanel*nDim + iDim];

                        }
                    }
                }
                iSample2++;
                if (iSample2< (nSample-nqSample)){

                    Read_TECPLOT_ASCII( config, geometry, iSample2, nqSample-1);
                    //SetSurfaceGeom(config, geometry, iSample2, nqSample - 1, iObserver_outer);
                    ComputeNormal(config, geometry, iSample2, nqSample - 1, iObserver_outer);
                    ComputeModifiedVelocity(config, iSample2, nqSample - 1);

                }
            }

        }//iSample Loop

        for (iObserver_inner=0; iObserver_inner<nObserver_inner; iObserver_inner++){
            for (iSample=0; iSample<nSample; iSample++){
                if(config->GetAcoustic_Inner_ObsLoop()){

                    su2double *source = nullptr;
                    if (rank == MASTER_NODE){
                        source = &pp_TimeDomainGlobal_inner[iObserver_inner][iSample];
                    }

                    SU2_MPI::Reduce( &pp_TimeDomain_inner[iObserver_inner][iSample] , source , 1 , MPI_DOUBLE , MPI_SUM ,MASTER_NODE, MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &t_interp_inner[iObserver_inner][iSample] , &t_interp_inner[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T1_I[iObserver_inner][iSample] , &T1_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T2_I[iObserver_inner][iSample] , &T2_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T3_I[iObserver_inner][iSample] , &T3_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T4_I[iObserver_inner][iSample] , &T4_I[iObserver_inner][iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                }else{

                    su2double *source = nullptr;
                    if (rank == MASTER_NODE){
                        source = &pp_TimeDomainGlobal_outer[iSample];
                    }

                    SU2_MPI::Reduce( &pp_TimeDomain_outer[iSample] , source , 1 , MPI_DOUBLE , MPI_SUM ,MASTER_NODE, MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &t_interp_outer[iSample] , &t_interp_outer[iSample] , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T1_O[iSample] , &T1_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T2_O[iSample] , &T2_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T3_O[iSample] , &T3_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                    SU2_MPI::Allreduce( &T4_O[iSample] , &T4_O[iSample] , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
                }
            }
        }

        if (rank==MASTER_NODE){
            for (iObserver_inner=0; iObserver_inner<nObserver_inner; iObserver_inner++){
                pp_mean=0;
                for (iSample=0; iSample<nSample; iSample++){
                    if(config->GetAcoustic_Inner_ObsLoop()){
                        iObserver=iObserver_inner;
                        pp_mean = (pp_TimeDomainGlobal_inner[iObserver][iSample]+pp_mean*iSample)/(iSample+1);
                    }else{
                        iObserver=iObserver_outer;
                        pp_mean = (pp_TimeDomainGlobal_outer[iSample]+pp_mean*iSample)/(iSample+1);
                    }
                }
                if (ceil(log2(nSample))==floor(log2(nSample))){
                    cout<<"FFT Transformation to obtain dB-Freq[Hz] "<<endl;
                    if(config->GetAcoustic_Inner_ObsLoop()){
                        FFT_AcousticPressureSignal(pp_mean,pp_TimeDomainGlobal_inner[iObserver],res_factor);
                        snprintf(buffer, sizeof(char) * 32, "pp_FWH_fft_%03d_Zone_%i", (iObserver_inner+1), (iZone));
                    }else{
                        FFT_AcousticPressureSignal(pp_mean,pp_TimeDomainGlobal_outer,res_factor);
                        snprintf(buffer, sizeof(char) * 32, "pp_FWH_fft_%03d_Zone_%i", (iObserver_outer+1), (iZone));
                    }
                    pp_fft_file.open(buffer);
                    cout<<"Frequency[Hz] vs dB data can be found in "<<buffer<<" for Observer-"<<iObserver+1<<endl;
                }

                SPL_iObserver = 0.0;
                snprintf(buffer, sizeof(char) * 32, "pp_FWH_%03d_Zone_%i", (iObserver+1), (iZone));
                pp_FWH_file.open(buffer);

                for (iSample=0; iSample<nSample; iSample++){
                    if (iSample==0){
                        pp_FWH_file <<  "Time"<<' '<<"P-Fluctuation"<<' '<<"Acoustic-Pressure"<<' '<<"Thickness"<<' '<<"Loading"<<' '<<"Term-1"<<' '<<"Term-2"<<' '<<"Term-3"<<' '<<"Term-4"<<endl;
                    }
                    if(config->GetAcoustic_Inner_ObsLoop()){
                        pp_FWH_file << std::setprecision(15) <<t_interp_inner[iObserver][iSample] <<' '<<pp_TimeDomainGlobal_inner[iObserver][iSample]-pp_mean<<' '<<pp_TimeDomainGlobal_inner[iObserver][iSample]<< ' ' <<T1_I[iObserver][iSample]+T2_I[iObserver][iSample]<<' '<< T3_I[iObserver][iSample]+T4_I[iObserver][iSample]<<' '<<T1_I[iObserver][iSample]<<' '<<T2_I[iObserver][iSample]<<' '<<T3_I[iObserver][iSample]<<' '<<T4_I[iObserver][iSample]<<' '<<endl;
                    }else{
                        pp_FWH_file << std::setprecision(15) <<t_interp_outer[iSample] <<' '<<pp_TimeDomainGlobal_outer[iSample]-pp_mean<<' '<<pp_TimeDomainGlobal_outer[iSample]<< ' ' <<T1_O[iSample]+T2_O[iSample]<<' '<<T3_O[iSample]+T4_O[iSample]<<' '<<T1_O[iSample]<<' '<<T2_O[iSample]<<' '<<T3_O[iSample]<<' '<<T4_O[iSample]<<' '<<endl;
                    }
                    if (ceil(log2(nSample))==floor(log2(nSample))){
                        if (iSample==0){
                            pp_fft_file << "Frequency[Hz]"<<" "<<"dB"<<endl;
                        }
                        if ((iSample<(nSample/2+1)) && (iSample>0)){
                            pp_fft_file << std::setprecision(15) << Fs*iSample/nSample<< ' '<< real(pp_fft[iSample])<< ' '<<endl;
                        }
                    }
                    if(config->GetAcoustic_Inner_ObsLoop()){
                        SPL_iObserver= SPL_iObserver + (pp_TimeDomainGlobal_inner[iObserver][iSample]-pp_mean)*(pp_TimeDomainGlobal_inner[iObserver][iSample]-pp_mean);
                    }else{
                        SPL_iObserver= SPL_iObserver + (pp_TimeDomainGlobal_outer[iSample]-pp_mean)*(pp_TimeDomainGlobal_outer[iSample]-pp_mean);
                    }
                }
                SPL_iObserver = sqrt(SPL_iObserver/nSample);
                SPL = SPL + SPL_iObserver;
                pp_FWH_file.close();
                pp_fft_file.close();
                cout<<std::setprecision(15)<<"RMS(p')-FWH over observer#"<<iObserver+1<<" -------------------> "<<SPL_iObserver <<"  *******"<<endl;
                cout<<"F1A Computation is finished for Observer#"<<iObserver+1<<". Tabulated data (Time-domain) can be found in "<< buffer<<" ."<<endl;
            }
        }
    }//iObserver_Outer Loop
    if (rank==MASTER_NODE){
        SPL = SPL/nObserver;
        cout<<endl<<std::setprecision(15)<<"****** RMS(p') averaged over "<<nObserver<<" observer locations = "<<SPL <<"  **************"<<endl;
    }

}



void F1A::F1A_Formulation ( CConfig *config,unsigned long iObserver, unsigned long iPanel,  unsigned long iSample, unsigned long iLocSample, unsigned long i){

#ifdef HAVE_MPI
    int rank, nProcessor;
        	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        	MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
    unsigned long iDim;
    su2double Lnr,LnM, L_dotr, Qn, n_dot_Q, Qdotn, M, Mr,Mr1, M_dotr,n_dot,M_dot,L_dot;
    su2double Qdot;
    su2double K,R,r,dS;
    su2double dt_src;
    dt_src=config->GetDelta_UnstTime()*SamplingFreq;
    Lnr=0;LnM=0; L_dotr=0; Qn=0; n_dot_Q=0;  Qdotn=0; M=0; Mr=0; M_dotr=0;
    Qdot = 0;
    Mr1 = 0;
    dS = Area[iPanel];

    for (iDim=0; iDim<nDim; iDim++){
        r = RadVec[5*i+iDim];
        M_dot=(FWH_Surf_Vel[(iLocSample+1)*nPanel*nDim+iPanel*nDim+iDim]-FWH_Surf_Vel[(iLocSample-1)*nPanel*nDim+iPanel*nDim+iDim])/2.0/dt_src/a_inf;
        M_dotr=M_dotr+M_dot*r;
        Mr=Mr+(FWH_Surf_Vel[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*r)/a_inf;
        Mr1=1-Mr;
        M=M + FWH_Surf_Vel[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*FWH_Surf_Vel[iLocSample*nPanel*nDim+iPanel*nDim+iDim]/a_inf/a_inf;
        Qdot=Qdot+ (Q[(iLocSample+1)*nPanel*nDim+iPanel*nDim+iDim]-Q[(iLocSample-1)*nPanel*nDim+iPanel*nDim+iDim])/2.0/dt_src;
        Qdotn=Qdotn+ Qdot*UnitaryNormal[iLocSample*nPanel*nDim+iPanel*nDim+iDim];
        n_dot=(UnitaryNormal[(iLocSample+1)*nPanel*nDim+iPanel*nDim+iDim]-UnitaryNormal[(iLocSample-1)*nPanel*nDim+iPanel*nDim+iDim])/2.0/dt_src;
        n_dot_Q=n_dot_Q+n_dot*Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim];
        Qn=Qn+Q[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*UnitaryNormal[iLocSample*nPanel*nDim+iPanel*nDim+iDim];
        L_dot=(F[(iLocSample+1)*nPanel*nDim+iPanel*nDim+iDim]-F[(iLocSample-1)*nPanel*nDim+iPanel*nDim+iDim])/2.0/dt_src;
        L_dotr= L_dotr+ L_dot*r;
        Lnr     =Lnr    + F[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*r;
        LnM     =LnM    + F[iLocSample*nPanel*nDim+iPanel*nDim+iDim]*FWH_Surf_Vel[iLocSample*nPanel*nDim+iPanel*nDim+iDim]/a_inf;
    }
    M  = sqrt( M );
    R  = RadVec[5*i+nDim];
    K=M_dotr*R+Mr*a_inf-M*M*a_inf;
    T1=(FreeStreamDensity*(Qdotn+n_dot_Q))/R/Mr1/Mr1*(dS/4.0/M_PI);
    T2= (FreeStreamDensity*Qn*K)/R/R/Mr1/Mr1/Mr1*(dS/4.0/M_PI);
    T3= L_dotr/R/Mr1/Mr1/a_inf*(dS/4.0/M_PI);
    T4= (Lnr-LnM)/R/R/Mr1/Mr1*(dS/4.0/M_PI);
    T5= (Lnr*K)/R/R/Mr1/Mr1/Mr1/a_inf*(dS/4.0/M_PI);
    pp_t= (T1+T2+T3+T4+T5);

}

void F1A::FFT_AcousticPressureSignal(su2double pp_mean, su2double *pp_TimeDomain, unsigned long res_factor){
    unsigned long iPanel, iSample;
    FFT* FFT_container = new FFT() ;
    su2double pLa;

    /*perform FFT on Acoustic Pressure Signal*/
    for (iSample=0; iSample<nSample/res_factor; iSample++){
        pp_fft[iSample].real(pp_TimeDomain[iSample*res_factor]-pp_mean);
    }
    /* Writing the complex array data*/
    CArray dataPP(pp_fft ,nSample/res_factor);
    /* Call the FFT function for performing the forward fft */
    FFT_container->fft_r2(dataPP);
    for (iSample=0; iSample<nSample/res_factor; iSample++){
        pLa = abs(real(dataPP[iSample]));
        pp_fft[iSample] = 20*log10(pLa/0.00002); // 0.00002 threshold of hearing
    }
}


void F1A::CombineZones( CConfig *config){

    string data;
    string ttR,ppR,pAcR,TNR,LNR,T1R,T2R,T3R,T4R;
    su2double tt,pAc,TT1,TT2,TT3,TT4;
    su2double pp_mean;
    ifstream ppFWH;
    ofstream CombinedppFWH;
    unsigned long iSample,iObserver,zone;

    char buffer[32]; // The filename buffer.

    su2double *PP    = new su2double [nSample];
    su2double *t     = new su2double [nSample];
    su2double *Term1 = new su2double [nSample];
    su2double *Term2 = new su2double [nSample];
    su2double *Term3 = new su2double [nSample];
    su2double *Term4 = new su2double [nSample];

    for (iObserver=0; iObserver<nObserver; iObserver++){
        for (iSample=0; iSample<nSample; iSample++){
            PP[iSample]=0.0;t[iSample]=0.0;
            Term1[iSample]=0.0;Term2[iSample]=0.0;Term3[iSample]=0.0;Term4[iSample]=0.0;
        }
        for (zone=0; zone < nZone; zone++){

            snprintf(buffer, sizeof(char) * 32, "pp_FWH_%03d_Zone_%i", (iObserver+1), (zone));
            ppFWH.open(buffer);
            cout<<"Reading -> "<<buffer<<endl;
            if (ppFWH.fail()) {
                cout << "There is no file!!! " <<  buffer  << "."<< endl;
                exit(EXIT_FAILURE);
            }
            getline(ppFWH,data); //to get rid of first line that is headerline
            ppFWH.close();
            ppFWH.open(buffer);
            getline(ppFWH,data); //to get rid of first line that is headerline
            iSample=0;
            while(getline(ppFWH,data)){
                istringstream point_data(data);
                getline(point_data,ttR,' ');
                getline(point_data,ppR,' ');
                getline(point_data,pAcR,' ');
                getline(point_data,TNR,' ');
                getline(point_data,LNR,' ');
                getline(point_data,T1R,' ');
                getline(point_data,T2R,' ');
                getline(point_data,T3R,' ');
                getline(point_data,T4R,' ');
                stringstream geek(ttR); geek >> tt;
                if (zone == 0) t[iSample]=tt;
                stringstream geek2(pAcR); geek2 >> pAc; PP[iSample]=PP[iSample]+pAc;
                stringstream geek3(T1R); geek3 >> TT1; Term1[iSample]=Term1[iSample]+TT1;
                stringstream geek4(T2R); geek4 >> TT2; Term2[iSample]=Term2[iSample]+TT2;
                stringstream geek5(T3R); geek5 >> TT3; Term3[iSample]=Term3[iSample]+TT3;
                stringstream geek6(T4R); geek6 >> TT4; Term4[iSample]=Term4[iSample]+TT4;
                iSample++;
            }
            ppFWH.close();
        }
        pp_mean=0;
        for (iSample=0; iSample<nSample; iSample++){
            pp_mean = (PP[iSample]+pp_mean*iSample)/(iSample+1);
        }
        snprintf(buffer, sizeof(char) * 32, "pp_FWH_%03d_Combined", (iObserver+1));
        CombinedppFWH.open(buffer);
        for (iSample=0; iSample<nSample; iSample++){
            if (iSample==0){
                CombinedppFWH  <<  "Time"<<' '<<"P-Fluctuation"<<' '<<"Acoustic-Pressure"<<' '<<"Thickness"<<' '<<"Loading"<<' '<<"Term-1"<<' '<<"Term-2"<<' '<<"Term-3"<<' '<<"Term-4"<<endl;
            }
            CombinedppFWH << std::setprecision(15) <<t[iSample] <<' '<<PP[iSample]-pp_mean<<' '<<PP[iSample]<< ' ' <<Term1[iSample]+Term2[iSample]<<' '<< Term3[iSample]+Term4[iSample]<<' '<<Term1[iSample]<<' '<<Term2[iSample]<<' '<<Term3[iSample]<<' '<<Term4[iSample]<<endl;
        }
        CombinedppFWH.close();
    }

    delete [] PP;
    delete [] t;
    delete [] Term1;
    delete [] Term2;
    delete [] Term3;
    delete [] Term4;

}

