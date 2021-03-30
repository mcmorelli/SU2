/*!
* \file precice.cpp
* \brief Adapter class for coupling SU2 with preCICE for FSI.
* \author Alexander Rusch
*/

#include "../include/precice.hpp"

Precice::Precice(
        const string &preciceConfigurationFileName,
        int solverProcessIndex,
        int solverProcessSize,
        CGeometry ***geometry_container,
        CSolver ****solver_container,
        CConfig **config_container,
        CVolumetricMovement **grid_movement)
        :
        solverProcessIndex(solverProcessIndex),
        solverProcessSize(solverProcessSize),
        solverInterface("SU2_CFD", preciceConfigurationFileName, solverProcessIndex, solverProcessSize),
        geometry_container(geometry_container),
        solver_container(solver_container),
        config_container(config_container),
        grid_movement(grid_movement),
        //For implicit coupling
        coric(precice::constants::actionReadIterationCheckpoint()),
        cowic(precice::constants::actionWriteIterationCheckpoint()) {

    nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();
    verbosityLevel_high = config_container[ZONE_0]->GetpreCICE_VerbosityLevel_High();
    nPoint = geometry_container[ZONE_0][MESH_0]->GetnPoint();
    nVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    globalNumberWetSurfaces = config_container[ZONE_0]->GetnMarker_PreCICE();

    localNumberWetSurfaces = 0;
    dt_savedState = 0;

    StopCalc_savedState = false;

    vertexIDs = nullptr;
    forceID = nullptr;
    displDeltaID = nullptr;
    meshID = nullptr;
    forces = nullptr;
    displacementDeltas = nullptr;
    processWorkingOnWetSurface = true;
    valueMarkerWet = nullptr;
    vertexSize = nullptr;
    indexMarkerWetMappingLocalToGlobal = nullptr;

    Coord_Saved = nullptr;
    Coord_n_Saved = nullptr;
    Coord_n1_Saved = nullptr;
    Coord_p1_Saved = nullptr;
    GridVel_Saved = nullptr;
    GridVel_Grad_Saved = nullptr;
    solution_Saved = nullptr;
    solution_time_n_Saved = nullptr;
    solution_time_n1_Saved = nullptr;

    Coord_Saved = new double *[nPoint];
    Coord_n_Saved = new double *[nPoint];
    Coord_n1_Saved = new double *[nPoint];
    Coord_p1_Saved = new double *[nPoint];
    GridVel_Saved = new double *[nPoint];
    GridVel_Grad_Saved = new double **[nPoint];
    solution_Saved = new double *[nPoint];
    solution_time_n_Saved = new double *[nPoint];
    solution_time_n1_Saved = new double *[nPoint];
    for (int iPoint = 0; iPoint < nPoint; iPoint++) {
        Coord_Saved[iPoint] = new double[nDim];
        Coord_n_Saved[iPoint] = new double[nDim];
        Coord_n1_Saved[iPoint] = new double[nDim];
        Coord_p1_Saved[iPoint] = new double[nDim];
        GridVel_Saved[iPoint] = new double[nDim];
        GridVel_Grad_Saved[iPoint] = new double *[nDim];
        for (int iDim = 0; iDim < nDim; iDim++) {
            GridVel_Grad_Saved[iPoint][iDim] = new double[nDim];
        }
        solution_Saved[iPoint] = new double[nVar];
        solution_time_n_Saved[iPoint] = new double[nVar];
        solution_time_n1_Saved[iPoint] = new double[nVar];
    }
}

Precice::~Precice(void) {

    for (int i = 0; i < localNumberWetSurfaces; i++) {
        delete[] vertexIDs[i];
    }
    delete[] vertexIDs;

    for (int iPoint = 0; iPoint < nPoint; iPoint++) {
        delete[] Coord_Saved[iPoint];
        delete[] Coord_n_Saved[iPoint];
        delete[] Coord_n1_Saved[iPoint];
        delete[] Coord_p1_Saved[iPoint];
        delete[] GridVel_Saved[iPoint];
        delete[] solution_Saved[iPoint];
        delete[] solution_time_n_Saved[iPoint];
        delete[] solution_time_n1_Saved[iPoint];
    }
    delete[] Coord_Saved;
    delete[] Coord_n_Saved;
    delete[] Coord_n1_Saved;
    delete[] Coord_p1_Saved;
    delete[] GridVel_Saved;
    delete[] solution_Saved;
    delete[] solution_time_n_Saved;
    delete[] solution_time_n1_Saved;

    for (int iPoint = 0; iPoint < nPoint; iPoint++) {
        for (int iDim = 0; iDim < nDim; iDim++) {
            delete[] GridVel_Grad_Saved[iPoint][iDim];
        }
        delete[] GridVel_Grad_Saved[iPoint];
    }
    delete[] GridVel_Grad_Saved;

    delete[] forceID;
    delete[] displDeltaID;
    delete[] valueMarkerWet;
    delete[] vertexSize;
    delete[] indexMarkerWetMappingLocalToGlobal;
    delete[] meshID;
}

double Precice::initialize(){

  /*--- preCICE timestep size ---*/
  double precice_dt;

  if (verbosityLevel_high) {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Initializing preCICE..." << endl;
  }

  CheckDimensionalConsistency();

  GetPreciceMeshID();

  SetnLocalPreciceMarkers();

  CheckWorkingProcess();

  SetMarkerMapping();

  SetMeshVertices();

  precice_dt = SetTimeStep();

  return precice_dt;
}

void Precice::CheckDimensionalConsistency() {

    //Checking for dimensional consistency of SU2 and preCICE - Exit if not consistent
    if(solverInterface.getDimensions() != geometry_container[ZONE_0][MESH_0]->GetnDim()){
        cout << "Dimensions of SU2 and preCICE are not equal! Now exiting..." << endl;
        exit(EXIT_FAILURE);
    }
}

void Precice::GetPreciceMeshID() {

    //Checking for number of wet surfaces - Exit if not cat least one wet surface defined
    if(globalNumberWetSurfaces < 1){
        cout << "There must be at least one wet surface! Now exiting..." << endl;
        exit(EXIT_FAILURE);
    } else {
        meshID = new int[globalNumberWetSurfaces];
        forceID = new int[globalNumberWetSurfaces];
        displDeltaID = new int[globalNumberWetSurfaces];
        for (int i = 0; i < globalNumberWetSurfaces; i++) {
            //Get preCICE meshIDs
            meshID[i] = solverInterface.getMeshID("SU2_Mesh" + to_string(i));
        }
    }
}

void Precice::SetnLocalPreciceMarkers() {

    unsigned short iMarker, jMarker;
    string Marker_Tag, PreCICE_Tag;

    for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
        if (config_container[ZONE_0]->GetMarker_All_Moving(iMarker) != YES) continue;

        Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);

        for (jMarker = 0; jMarker < config_container[ZONE_0]->GetnMarker_PreCICE(); jMarker++) {

            PreCICE_Tag = config_container[ZONE_0]->GetMarker_PreCICE_TagBound(jMarker);

            if (Marker_Tag != PreCICE_Tag) {
                continue;
            }

            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << " PreCICE_Tag: " << PreCICE_Tag << ", Marker_Tag: " <<  Marker_Tag << endl;
            localNumberWetSurfaces++;
        }
    }
}

void Precice::CheckWorkingProcess() {

    if (localNumberWetSurfaces < 1) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Does not work on the wet surface at all." << endl;
        processWorkingOnWetSurface = false;
    }
}

void Precice::SetMarkerMapping() {

    unsigned short iMarker, jMarker;
    string Marker_Tag, PreCICE_Tag;

    if (processWorkingOnWetSurface) {

        //Store the wet surface marker values in an array, which has the size equal to the number of wet surfaces actually being worked on by this process
        valueMarkerWet = new short[localNumberWetSurfaces];
        indexMarkerWetMappingLocalToGlobal = new short[localNumberWetSurfaces];
        int markerIndex = 0;

        for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {

            if (config_container[ZONE_0]->GetMarker_All_Moving(iMarker) != YES) continue;

            Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);

            for (jMarker = 0; jMarker < config_container[ZONE_0]->GetnMarker_PreCICE(); jMarker++) {

                PreCICE_Tag = config_container[ZONE_0]->GetMarker_PreCICE_TagBound(jMarker);

                if (Marker_Tag != PreCICE_Tag) {
                    continue;
                }

                valueMarkerWet[markerIndex] = config_container[ZONE_0]->GetMarker_All_TagBound(Marker_Tag);
                indexMarkerWetMappingLocalToGlobal[markerIndex] = config_container[ZONE_0]->GetMarker_PreCICE(PreCICE_Tag);
                markerIndex++;
            }
        }
        vertexIDs = new int*[localNumberWetSurfaces];
    }
}

void Precice::SetMeshVertices() {

    if (processWorkingOnWetSurface) {
        vertexSize = new unsigned long[localNumberWetSurfaces];
        for (int i = 0; i < localNumberWetSurfaces; i++) {
            vertexSize[i] = geometry_container[ZONE_0][MESH_0]->nVertex[valueMarkerWet[i]];

            unsigned long iNode;  /*--- variable for storing the node indices - one at the time ---*/
            su2double* Coord = nullptr;

            auto* coords = new double [vertexSize[i]*nDim];

            for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
                iNode = geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode();
                Coord = geometry_container[ZONE_0][MESH_0]->nodes->GetCoord(iNode);

                for (int iDim = 0; iDim < nDim; iDim++) {
                    coords[iVertex*nDim + iDim] = Coord[iDim];
                }
            }

            //preCICE internal
            vertexIDs[i] = new int[vertexSize[i]];
            solverInterface.setMeshVertices(meshID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], coords, vertexIDs[i]);
            forceID[indexMarkerWetMappingLocalToGlobal[i]] = solverInterface.getDataID("Forces" + to_string(indexMarkerWetMappingLocalToGlobal[i]), meshID[indexMarkerWetMappingLocalToGlobal[i]]);
            displDeltaID[indexMarkerWetMappingLocalToGlobal[i]] = solverInterface.getDataID("DisplacementDeltas" + to_string(indexMarkerWetMappingLocalToGlobal[i]), meshID[indexMarkerWetMappingLocalToGlobal[i]]);

            delete [] coords;

        }
        for (int i = 0; i < globalNumberWetSurfaces; i++) {
            bool flag = false;
            for (int j = 0; j < localNumberWetSurfaces; j++) {
                if (indexMarkerWetMappingLocalToGlobal[j] == i) {
                    flag = true;
                }
            }
            if (!flag) {
                solverInterface.setMeshVertices(meshID[i], 0, nullptr, nullptr);
                forceID[i] = solverInterface.getDataID("Forces" + to_string(i), meshID[i]);
                displDeltaID[i] = solverInterface.getDataID("DisplacementDeltas" + to_string(i), meshID[i]);
            }
        }
    } else {
        for (int i = 0; i < globalNumberWetSurfaces; i++) {
            solverInterface.setMeshVertices(meshID[i], 0, nullptr, nullptr);
            forceID[i] = solverInterface.getDataID("Forces" + to_string(i), meshID[i]);
            displDeltaID[i] = solverInterface.getDataID("DisplacementDeltas" + to_string(i), meshID[i]);
        }
    }
}

double Precice::SetTimeStep() {

    double time_step;

    time_step = solverInterface.initialize();
    if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": ...done initializing preCICE!" << endl;
    }

    return time_step;
}

double Precice::advance(double computedTimestepLength) {

    double max_precice_dt;

    if (processWorkingOnWetSurface) {

        ComputeForces();

        //3. Advance solverInterface
        max_precice_dt = ComputeMaxTimeStep(computedTimestepLength);

        SetDisplacements();

        return max_precice_dt;
    } else {
        //3. Advance solverInterface
        max_precice_dt = ComputeMaxTimeStep(computedTimestepLength);

        return max_precice_dt;
    }
}

void Precice::ComputeForces() {

    //Get physical simulation information
    bool incompressible = (config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
    bool viscous_flow = ((config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES) || (config_container[ZONE_0]->GetKind_Solver() == RANS));

    //Compute factorForces for redimensionalizing forces ("ND" = Non-Dimensional)
    double* Velocity_Real = config_container[ZONE_0]->GetVelocity_FreeStream();
    double Density_Real = config_container[ZONE_0]->GetDensity_FreeStream();
    double* Velocity_ND = config_container[ZONE_0]->GetVelocity_FreeStreamND();
    double Density_ND = config_container[ZONE_0]->GetDensity_FreeStreamND();
    double Velocity2_Real = 0.0;  /*--- denotes squared real velocity ---*/
    double Velocity2_ND = 0.0;  /*--- denotes squared non-dimensional velocity ---*/
    //Compute squared values
    for (int iDim = 0; iDim < nDim; iDim++){
        Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
        Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
    }
    //Compute factor for redimensionalizing forces
    double factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);
    if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Factor for (non-/re-)dimensionalization of forces: " << factorForces << endl;  /*--- for debugging purposes ---*/
    }

    for (int i = 0; i < localNumberWetSurfaces; i++) {
        if (verbosityLevel_high) {
            //1. Compute forces
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Computing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
        }

        auto* nodeVertex = new unsigned long [vertexSize[i]];

        auto** normalsVertex = new double* [vertexSize[i]];
        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++){
            normalsVertex[iVertex] = new double [nDim];
        }

        auto** normalsVertex_Unit = new double* [vertexSize[i]];
        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++){
            normalsVertex_Unit[iVertex] = new double [nDim];
        }

        double Area;
        double Pn = 0.0;  /*--- denotes pressure at a node ---*/
        double Pinf = 0.0;  /*--- denotes environmental (farfield) pressure ---*/
        double** Grad_PrimVar = nullptr; /*--- denotes (u.A. velocity) gradients needed for computation of viscous forces ---*/
        double Viscosity = 0.0;
        double Tau[3][3];
        double TauElem[3];
        //double forces_su2[vertexSize[i]][nDim];  /*--- forces will be stored such, before converting to simple array ---*/
        auto** forces_su2 = new double* [vertexSize[i]];
        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++){
            forces_su2[iVertex] = new double [nDim];
        }

        /*--- Loop over vertices of coupled boundary ---*/
        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {

            //Get node number (= index) to vertex (= node)
            nodeVertex[iVertex] = geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode(); /*--- Store all nodes (indices) in a vector ---*/
            // Get normal vector
            for (int iDim = 0; iDim < nDim; iDim++){
                normalsVertex[iVertex][iDim] = (geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNormal())[iDim];
            }
            // Unit normals
            Area = 0.0;
            for (int iDim = 0; iDim < nDim; iDim++) {
                Area += normalsVertex[iVertex][iDim]*normalsVertex[iVertex][iDim];
            }
            Area = sqrt(Area);
            for (int iDim = 0; iDim < nDim; iDim++) {
                normalsVertex_Unit[iVertex][iDim] = normalsVertex[iVertex][iDim]/Area;
            }
            // Get the values of pressure and viscosity
            CVariable* flow_nodes = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetNodes();
            Pinf = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetPressure_Inf();
            Pn = flow_nodes->GetPressure(nodeVertex[iVertex]);

            if (viscous_flow){
                Viscosity = flow_nodes->GetLaminarViscosity(nodeVertex[iVertex]);
                Grad_PrimVar = flow_nodes->GetGradient_Primitive(nodeVertex[iVertex]);

            }

            // Calculate the forces_su2 in the nodes for the inviscid term --> Units of force (non-dimensional).
            for (int iDim = 0; iDim < nDim; iDim++) {
                forces_su2[iVertex][iDim] = -(Pn-Pinf)*normalsVertex[iVertex][iDim];
            }
            // Calculate the forces_su2 in the nodes for the viscous term
            if (viscous_flow){
                // Divergence of the velocity
                double div_vel = 0.0;
                for (int iDim = 0; iDim < nDim; iDim++){
                    div_vel += Grad_PrimVar[iDim+1][iDim];
                }
                if (incompressible){
                    div_vel = 0.0;  /*--- incompressible flow is divergence-free ---*/
                }
                for (int iDim = 0; iDim < nDim; iDim++) {
                    for (int jDim = 0 ; jDim < nDim; jDim++) {
                        // Dirac delta
                        double Delta = 0.0;
                        if (iDim == jDim){
                            Delta = 1.0;
                        }
                        // Viscous stress
                        Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
                                          2.0/3.0*Viscosity*div_vel*Delta;
                        // Add Viscous component in the forces_su2 vector --> Units of force (non-dimensional).
                        forces_su2[iVertex][iDim] += Tau[iDim][jDim]*normalsVertex[iVertex][jDim];
                    }
                }
            }
            // Rescale forces_su2 to SI units
            for (int iDim = 0; iDim < nDim; iDim++) {
                forces_su2[iVertex][iDim] = forces_su2[iVertex][iDim]*factorForces;
            }
        }
        //convert forces_su2 into forces
        forces = new double[vertexSize[i]*nDim];


        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
            for (int iDim = 0; iDim < nDim; iDim++) {
                //Do not write forces for duplicate nodes! -> Check wether the color of the node matches the MPI-rank of this process. Only write forces, if node originally belongs to this process.
                if (geometry_container[ZONE_0][MESH_0]->nodes->GetColor(iVertex) == solverProcessIndex) {
                    forces[iVertex*nDim + iDim] = forces_su2[iVertex][iDim];
                }
                else{
                    forces[iVertex*nDim + iDim] = 0;
                }
            }
        }
        if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done computing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << endl;
        }

        //2. Write forces
        if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Writing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
        }
        //Load Ramping functionality: Reduce force vector before transferring by a ramping factor, which increases with the number of elapsed time steps; Achtung: ExtIter beginnt bei 0 (ohne Restart) und bei einem Restart (Startlösung) nicht bei 0, sondern bei der Startiterationsnummer
        if (config_container[ZONE_0]->GetpreCICE_LoadRamping() && ((config_container[ZONE_0]->GetnTime_Iter() - config_container[ZONE_0]->GetRestart_Iter()) < config_container[ZONE_0]->GetpreCICE_LoadRampingDuration())) {
            if (verbosityLevel_high) {
                cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Load ramping factor in preCICE: " << config_container[ZONE_0]->GetnTime_Iter() - config_container[ZONE_0]->GetRestart_Iter() + 1 << "/" << config_container[ZONE_0]->GetpreCICE_LoadRampingDuration() << endl;
            }
            *forces = *forces * ((config_container[ZONE_0]->GetnTime_Iter() - config_container[ZONE_0]->GetRestart_Iter()) + 1) / config_container[ZONE_0]->GetpreCICE_LoadRampingDuration();
        }
        solverInterface.writeBlockVectorData(forceID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], vertexIDs[i], forces);
        if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done writing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "." << endl;
        }

        delete [] forces;

        delete [] nodeVertex;

        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
            delete [] normalsVertex[iVertex];
        }
        delete [] normalsVertex;

        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
            delete [] normalsVertex_Unit[iVertex];
        }
        delete [] normalsVertex_Unit;

        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
            delete [] forces_su2[iVertex];
        }
        delete [] forces_su2;

    }
}

void Precice::SetDisplacements() {

    for (int i = 0; i < localNumberWetSurfaces; i++) {
        //4. Read displacements/displacementDeltas
        if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Reading displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
        }

        //double displacementDeltas_su2[vertexSize[i]][nDim]; /*--- displacementDeltas will be stored such, before converting to simple array ---*/
        auto** displacementDeltas_su2 = new double* [vertexSize[i]];
        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++){
            displacementDeltas_su2[iVertex] = new double [nDim];
        }

        displacementDeltas = new double[vertexSize[i]*nDim];

        solverInterface.readBlockVectorData(displDeltaID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], vertexIDs[i], displacementDeltas);
        if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done reading displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "." << endl;
        }

        //cout << "displacementDeltas" << endl;
        //5. Set displacements/displacementDeltas
        if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Setting displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
        }
        //convert displacementDeltas into displacementDeltas_su2
        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
            for (int iDim = 0; iDim < nDim; iDim++) {
                displacementDeltas_su2[iVertex][iDim] = displacementDeltas[iVertex*nDim + iDim];
            }
        }

        unsigned long iNode;
        //CVariable* nodes = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetNodes();
        CVariable* nodes = solver_container[ZONE_0][MESH_0][MESH_SOL]->GetNodes();
        su2double VarCoordAbs[3] = {0.0};

        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {

            iNode = geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode();

            for (int iDim = 0; iDim < nDim; iDim++) {
                VarCoordAbs[iDim] = nodes->GetBound_Disp(iNode, iDim) + displacementDeltas_su2[iVertex][iDim];
            }

            nodes->SetBound_Disp(iNode, VarCoordAbs);
        }

        if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done setting displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "." << endl;
        }

        delete [] displacementDeltas;

        for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
            delete[] displacementDeltas_su2[iVertex];
        }
        delete [] displacementDeltas_su2;

    }
}

double Precice::ComputeMaxTimeStep(double computedTimestepLength) {

    double max_timestep;

    max_timestep = solverInterface.advance( computedTimestepLength );

    return max_timestep;
}

bool Precice::isCouplingOngoing(){
  return solverInterface.isCouplingOngoing();
}

bool Precice::isActionRequired( const string& action ){
  return solverInterface.isActionRequired(action);
}

const string& Precice::getCowic(){
  return cowic;
}

const string& Precice::getCoric(){
  return coric;
}

void Precice::saveOldState( bool *StopCalc, double *dt ){

  CVariable* flow_nodes = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetNodes();
  auto geometry_nodes = geometry_container[ZONE_0][MESH_0]->nodes;
  auto& GridVelGrad = geometry_nodes->GetGridVel_Grad();

  for (int iPoint = 0; iPoint < nPoint; iPoint++) {

    //Save solutions at last and current time step
    solution_Saved[iPoint] = flow_nodes->GetSolution(iPoint);
    solution_time_n_Saved[iPoint]= flow_nodes->GetSolution_time_n(iPoint);
    solution_time_n1_Saved[iPoint] = flow_nodes->GetSolution_time_n1(iPoint);

    //Save coordinates at last, current and next time step
    Coord_Saved[iPoint] = geometry_nodes->GetCoord(iPoint);
    Coord_n_Saved[iPoint] = geometry_nodes->GetCoord_n(iPoint);
    Coord_n1_Saved[iPoint] = geometry_nodes->GetCoord_n1(iPoint);
    Coord_p1_Saved[iPoint] = geometry_nodes->GetCoord_p1(iPoint);

    for (int iDim = 0; iDim < nDim; iDim++) {
      //Save grid velocity
      GridVel_Saved[iPoint] = geometry_nodes->GetGridVel(iPoint);
      for (int jDim = 0; jDim < nDim; jDim++) {
        //Save grid velocity gradient
        GridVel_Grad_Saved[iPoint][iDim][jDim] = GridVelGrad[iPoint][iDim][jDim];
      }
    }
  }

  //Save wether simulation should be stopped after the current iteration
  StopCalc_savedState = *StopCalc;
  //Save the time step size
  dt_savedState = *dt;
  //Writing task has been fulfilled successfully
  solverInterface.markActionFulfilled(cowic);
}

void Precice::reloadOldState(bool *StopCalc, double *dt){

    CVariable* flow_nodes = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetNodes();
    auto geometry_nodes = geometry_container[ZONE_0][MESH_0]->nodes;
    auto& gridVelGrad = geometry_nodes->GetGridVel_Grad();

    for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    //Reload solutions at last and current time step
    flow_nodes->SetSolution(iPoint, solution_Saved[iPoint]);
    flow_nodes->Set_Solution_time_n(iPoint, solution_time_n_Saved[iPoint]);
    flow_nodes->Set_Solution_time_n1(iPoint, solution_time_n1_Saved[iPoint]);

    //Reload coordinates at last, current and next time step
    geometry_nodes->SetCoord(iPoint, Coord_Saved[iPoint]);
    geometry_nodes->SetCoord_n(iPoint, Coord_n_Saved[iPoint]);
    geometry_nodes->SetCoord_n1(iPoint, Coord_n1_Saved[iPoint]);
    geometry_nodes->SetCoord_p1(iPoint, Coord_p1_Saved[iPoint]);

    //Reload grid velocity
    geometry_nodes->SetGridVel(iPoint, GridVel_Saved[iPoint]);

    //Reload grid velocity gradient
    for (int iDim = 0; iDim < nDim; iDim++) {
      for (int jDim = 0; jDim < nDim; jDim++) {
        gridVelGrad[iPoint][iDim][jDim] =  GridVel_Grad_Saved[iPoint][iDim][jDim];
      }
    }

  }

  //Reload wether simulation should be stopped after current iteration
  *StopCalc = StopCalc_savedState;
  //Reload the time step size
  *dt = dt_savedState;
  //Reading task has been fulfilled successfully
  solverInterface.markActionFulfilled(coric);
}

void Precice::finalize(){
  cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Finalizing preCICE..." << endl;
  solverInterface.finalize();
  cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Done finalizing preCICE!" << endl;
}
