/*!
* \file precice.cpp
* \brief Adapter class for coupling SU2 with preCICE for FSI.
*/

#include "../include/precice.hpp"
#include "../../Common/include/toolboxes/geometry_toolbox.hpp"

// #include "../include/numerics/CNumerics.hpp"

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
    nGlobalMarkers = config_container[ZONE_0]->GetnMarker_PreCICE();

    nLocalMarkers = 0;
    dt_savedState = 0;

    StopCalc_savedState = false;

    VertexID = nullptr;
    ForceID = nullptr;
    DisplacementDeltaID = nullptr;
    MeshID = nullptr;
    Forces = nullptr;
    DisplacementDeltas = nullptr;
    workingProcess = true;
    Marker = nullptr;
    nVerticesOfMarker = nullptr;
    localToGlobalMapping = nullptr;

    Coord_Saved = nullptr;
    Coord_n_Saved = nullptr;
    Coord_n1_Saved = nullptr;
    Coord_p1_Saved = nullptr;
    GridVel_Saved = nullptr;
    GridVel_Grad_Saved = nullptr;
    solution_Saved = nullptr;
    solution_time_n_Saved = nullptr;
    solution_time_n1_Saved = nullptr;
    Displacement_Saved = nullptr;

    Coord_Saved = new double *[nPoint];
    Coord_n_Saved = new double *[nPoint];
    Coord_n1_Saved = new double *[nPoint];
    Coord_p1_Saved = new double *[nPoint];
    GridVel_Saved = new double *[nPoint];
    GridVel_Grad_Saved = new double **[nPoint];
    solution_Saved = new double *[nPoint];
    solution_time_n_Saved = new double *[nPoint];
    solution_time_n1_Saved = new double *[nPoint];
    Displacement_Saved = new double *[nPoint];
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
        Displacement_Saved[iPoint] = new double[nDim];
    }
}

Precice::~Precice(void) {

    for (int i = 0; i < nLocalMarkers; i++) {
        delete[] VertexID[i];
    }
    delete[] VertexID;

    for (int iPoint = 0; iPoint < nPoint; iPoint++) {
        delete[] Coord_Saved[iPoint];
        delete[] Coord_n_Saved[iPoint];
        delete[] Coord_n1_Saved[iPoint];
        delete[] Coord_p1_Saved[iPoint];
        delete[] GridVel_Saved[iPoint];
        delete[] solution_Saved[iPoint];
        delete[] solution_time_n_Saved[iPoint];
        delete[] solution_time_n1_Saved[iPoint];
        delete[] Displacement_Saved[iPoint];
    }
    delete[] Coord_Saved;
    delete[] Coord_n_Saved;
    delete[] Coord_n1_Saved;
    delete[] Coord_p1_Saved;
    delete[] GridVel_Saved;
    delete[] solution_Saved;
    delete[] solution_time_n_Saved;
    delete[] solution_time_n1_Saved;
    delete[] Displacement_Saved;

    for (int iPoint = 0; iPoint < nPoint; iPoint++) {
        for (int iDim = 0; iDim < nDim; iDim++) {
            delete[] GridVel_Grad_Saved[iPoint][iDim];
        }
        delete[] GridVel_Grad_Saved[iPoint];
    }
    delete[] GridVel_Grad_Saved;

    delete[] ForceID;
    delete[] DisplacementDeltaID;
    delete[] Marker;
    delete[] nVerticesOfMarker;
    delete[] localToGlobalMapping;
    delete[] MeshID;
}

double Precice::initialize() {

    // precice timestep size
    double precice_dt;

    //Check that both precice and SU2 are of the same dimensions.
    CheckDimensionalConsistency();

    //Get the ID belonging to the mesh with given name.
    GetPreciceMeshID();

    //Set the number of local precice markers for this process.
    SetnLocalPreciceMarkers();

    //Check if this process is the working process.
    CheckWorkingProcess();

    //Set the local and global marker mapping.
    SetMarkerMapping();

    //Set the mesh used for the coupling and data for exchange.
    SetMeshVertices();

    //Parallel communication to the coupling partner/s is setup.
    precice_dt = SetTimeStep();

    return precice_dt;
}

void Precice::CheckDimensionalConsistency() {

    //Checking for dimensional consistency of SU2 and preCICE - Exit if not consistent
    if (solverInterface.getDimensions() != geometry_container[ZONE_0][MESH_0]->GetnDim()) {
        cout << "Dimensions of SU2 and preCICE are not equal! Now exiting..." << endl;
        exit(EXIT_FAILURE);
    }
}

void Precice::GetPreciceMeshID() {

    //Checking for the total number of markers - Exit if not cat least one wet surface defined
    if (nGlobalMarkers < 1) {
        cout << "There must be at least one wet surface! Now exiting..." << endl;
        exit(EXIT_FAILURE);
    } else {
        MeshID = new int[nGlobalMarkers];
        ForceID = new int[nGlobalMarkers];
        DisplacementDeltaID = new int[nGlobalMarkers];
        for (int i = 0; i < nGlobalMarkers; i++) {
            //Get preCICE meshIDs
            MeshID[i] = solverInterface.getMeshID("SU2_Mesh" + to_string(i));
        }
    }
}

void Precice::SetnLocalPreciceMarkers() {

    //Checking for the number of local markers
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

            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize - 1 << " PreCICE_Tag: " << PreCICE_Tag
                 << ", Marker_Tag: " << Marker_Tag << endl;
            nLocalMarkers++;
        }
    }
}

void Precice::CheckWorkingProcess() {

    if (nLocalMarkers < 1) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize - 1
             << ": Does not work on the wet surface at all." << endl;
        workingProcess = false;
    }
}

void Precice::SetMarkerMapping() {

    unsigned short iMarker, jMarker;
    string Marker_Tag, PreCICE_Tag;

    if (workingProcess) {

        //Store the wet surface marker values in an array, which has the size equal to the number of wet surfaces actually being worked on by this process
        Marker = new short[nLocalMarkers];
        localToGlobalMapping = new short[nLocalMarkers];
        int markerIndex = 0;

        for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {

            if (config_container[ZONE_0]->GetMarker_All_Moving(iMarker) != YES) continue;

            Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);

            for (jMarker = 0; jMarker < config_container[ZONE_0]->GetnMarker_PreCICE(); jMarker++) {

                PreCICE_Tag = config_container[ZONE_0]->GetMarker_PreCICE_TagBound(jMarker);

                if (Marker_Tag != PreCICE_Tag) {
                    continue;
                }

                Marker[markerIndex] = config_container[ZONE_0]->GetMarker_All_TagBound(Marker_Tag);
                localToGlobalMapping[markerIndex] = config_container[ZONE_0]->GetMarker_PreCICE(PreCICE_Tag);
                markerIndex++;
            }
        }
        VertexID = new int *[nLocalMarkers];
    }
}

void Precice::SetMeshVertices() {

    if (workingProcess) {
        nVerticesOfMarker = new unsigned long[nLocalMarkers];
        for (int i = 0; i < nLocalMarkers; i++) {
            nVerticesOfMarker[i] = geometry_container[ZONE_0][MESH_0]->nVertex[Marker[i]];

            unsigned long iNode;  /*--- variable for storing the node indices - one at the time ---*/
            su2double *iCoord = nullptr;

            auto *Coord = new double[nVerticesOfMarker[i] * nDim];

            for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
                iNode = geometry_container[ZONE_0][MESH_0]->vertex[Marker[i]][iVertex]->GetNode();
                iCoord = geometry_container[ZONE_0][MESH_0]->nodes->GetCoord(iNode);

                for (int iDim = 0; iDim < nDim; iDim++) {
                    Coord[iVertex * nDim + iDim] = iCoord[iDim];
                }
            }

            //preCICE internal
            VertexID[i] = new int[nVerticesOfMarker[i]];
            solverInterface.setMeshVertices(MeshID[localToGlobalMapping[i]], nVerticesOfMarker[i], Coord, VertexID[i]);
            ForceID[localToGlobalMapping[i]] = solverInterface.getDataID("Forces" + to_string(localToGlobalMapping[i]),
                                                                         MeshID[localToGlobalMapping[i]]);
            DisplacementDeltaID[localToGlobalMapping[i]] = solverInterface.getDataID(
                    "DisplacementDeltas" + to_string(localToGlobalMapping[i]), MeshID[localToGlobalMapping[i]]);

            delete[] Coord;

        }
        for (int i = 0; i < nGlobalMarkers; i++) {
            bool flag = false;
            for (int j = 0; j < nLocalMarkers; j++) {
                if (localToGlobalMapping[j] == i) {
                    flag = true;
                }
            }
            if (!flag) {
                solverInterface.setMeshVertices(MeshID[i], 0, nullptr, nullptr);
                ForceID[i] = solverInterface.getDataID("Forces" + to_string(i), MeshID[i]);
                DisplacementDeltaID[i] = solverInterface.getDataID("DisplacementDeltas" + to_string(i), MeshID[i]);
            }
        }
    } else {
        for (int i = 0; i < nGlobalMarkers; i++) {
            solverInterface.setMeshVertices(MeshID[i], 0, nullptr, nullptr);
            ForceID[i] = solverInterface.getDataID("Forces" + to_string(i), MeshID[i]);
            DisplacementDeltaID[i] = solverInterface.getDataID("DisplacementDeltas" + to_string(i), MeshID[i]);
        }
    }
}

double Precice::SetTimeStep() {

    double time_step;

    time_step = solverInterface.initialize();
    if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize - 1 << ": ...done initializing preCICE!"
             << endl;
    }

    return time_step;
}

double Precice::advance(double computedTimestepLength) {

    double max_precice_dt;

    if (workingProcess) {

        // Compute the forces on the coupled surface using the SU2 solution.
        ComputeForces();

        //3. Advance solverInterface
        max_precice_dt = ComputeMaxTimeStep(computedTimestepLength);

        // Read and write the displacements given by the structural solver.
        SetDisplacements();

        return max_precice_dt;
    } else {
        //3. Advance solverInterface
        max_precice_dt = ComputeMaxTimeStep(computedTimestepLength);

        return max_precice_dt;
    }
}

void Precice::ComputeForces() {

    /*--- Get physical simulation information ---*/
    bool viscous_flow = ((config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES) ||
                         (config_container[ZONE_0]->GetKind_Solver() == RANS));

    auto config = config_container[ZONE_0];

    for (int i = 0; i < nLocalMarkers; i++) {

        unsigned long iNode;

        su2double *Normal = nullptr;
        su2double UnitNormal[3] = {0.0};

        su2double Area;
        su2double Pressure, Pressure_Freestream, Viscosity;
        su2double Tau[3][3], TauElem[3], Grad_Vel[3][3] = {{0.0}};

        auto **Total_Forces = new double *[nVerticesOfMarker[i]];
        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
            Total_Forces[iVertex] = new double[nDim];
        }

        auto **Pressure_Forces = new double *[nVerticesOfMarker[i]];
        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
            Pressure_Forces[iVertex] = new double[nDim];
        }

        auto **Friction_Forces = new double *[nVerticesOfMarker[i]];
        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
            Friction_Forces[iVertex] = new double[nDim];
        }

        CVariable *flow_nodes = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetNodes();
        Pressure_Freestream = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetPressure_Inf();

        /*--- Loop over vertices of coupled boundary ---*/
        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {

            iNode = geometry_container[ZONE_0][MESH_0]->vertex[Marker[i]][iVertex]->GetNode();

            Normal = geometry_container[ZONE_0][MESH_0]->vertex[Marker[i]][iVertex]->GetNormal();
            Area = GeometryToolbox::Norm(nDim, Normal);
            for (int iDim = 0; iDim < nDim; iDim++) {
                UnitNormal[iDim] = Normal[iDim] / Area;
            }

            /*--- Get the value of pressure on the node ---*/
            Pressure = flow_nodes->GetPressure(iNode);

            /*--- Pressure force computation ---*/
            for (int iDim = 0; iDim < nDim; iDim++) {
                Pressure_Forces[iVertex][iDim] = -(Pressure - Pressure_Freestream) * Normal[iDim];
            }

            if (viscous_flow) {

                /*--- Get the value of viscosity on the node ---*/
                Viscosity = flow_nodes->GetLaminarViscosity(iNode);

                for (int iDim = 0; iDim < nDim; iDim++) {
                    for (int jDim = 0; jDim < nDim; jDim++) {
                       Grad_Vel[iDim][jDim] = flow_nodes->GetGradient_Primitive(iNode, iDim + 1, jDim);
                    }
                }

                /*--- Evaluate Tau ---*/
                CNumerics::ComputeStressTensor(nDim, Tau, Grad_Vel, Viscosity);

                /*--- If necessary evaluate the QCR contribution to Tau ---*/
                bool QCR = config->GetQCR();
                if (QCR) CNumerics::AddQCR(nDim, Grad_Vel, Tau);

                /*--- Project Tau in each surface element ---*/
                for (int iDim = 0; iDim < nDim; iDim++) {
                    TauElem[iDim] = 0.0;
                    for (int jDim = 0; jDim < nDim; jDim++) {
                        TauElem[iDim] += Tau[iDim][jDim] * UnitNormal[jDim];
                    }
                }

                /*--- Viscous force computation ---*/
                for (int iDim = 0; iDim < nDim; iDim++) {
                    Friction_Forces[iVertex][iDim] = TauElem[iDim] * Area;
                }
            }
            else {
                for (int iDim = 0; iDim < nDim; iDim++) {
                    Friction_Forces[iVertex][iDim] = 0.0;
                }
            }

            /*--- Sum the total contribution of the forces ---*/
            for (int iDim = 0; iDim < nDim; iDim++) {
                Total_Forces[iVertex][iDim] = Pressure_Forces[iVertex][iDim] + Friction_Forces[iVertex][iDim];
            }

        }

        /*--- Write the forces to the precice solver interface ---*/
        Forces = new double[nVerticesOfMarker[i] * nDim];

        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
            iNode = geometry_container[ZONE_0][MESH_0]->vertex[Marker[i]][iVertex]->GetNode();

            for (int iDim = 0; iDim < nDim; iDim++) {
                if (geometry_container[ZONE_0][MESH_0]->nodes->GetColor(iNode) == solverProcessIndex) {
                    Forces[iVertex * nDim + iDim] = Total_Forces[iVertex][iDim];
                } else {
                    Forces[iVertex * nDim + iDim] = 0;
                }
            }
        }

        solverInterface.writeBlockVectorData(ForceID[localToGlobalMapping[i]], nVerticesOfMarker[i], VertexID[i],
                                             Forces);

        delete[] Forces;

        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
            delete[] Total_Forces[iVertex];
        }
        delete[] Total_Forces;

        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
            delete[] Pressure_Forces[iVertex];
        }
        delete[] Pressure_Forces;

        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {
            delete[] Friction_Forces[iVertex];
        }
        delete[] Friction_Forces;
    }
}

void Precice::SetDisplacements() {

    for (int i = 0; i < nLocalMarkers; i++) {

        DisplacementDeltas = new double[nVerticesOfMarker[i] * nDim];

        solverInterface.readBlockVectorData(DisplacementDeltaID[localToGlobalMapping[i]], nVerticesOfMarker[i],
                                            VertexID[i], DisplacementDeltas);

        unsigned long iNode;
        CVariable *nodes = solver_container[ZONE_0][MESH_0][MESH_SOL]->GetNodes();

        for (int iVertex = 0; iVertex < nVerticesOfMarker[i]; iVertex++) {

            su2double displacement[3] = {0.0, 0.0, 0.0};
            iNode = geometry_container[ZONE_0][MESH_0]->vertex[Marker[i]][iVertex]->GetNode();

            if (geometry_container[ZONE_0][MESH_0]->nodes->GetColor(iNode) == solverProcessIndex) {

                //Get the absolute boundary displacement relative to the initial undeformed mesh
                for (int iDim = 0; iDim < nDim; iDim++) {
                    displacement[iDim] = nodes->GetBound_Disp(iNode, iDim) + DisplacementDeltas[iVertex * nDim + iDim];
                }

                nodes->SetBound_Disp(iNode, displacement);
            }
        }

        delete[] DisplacementDeltas;
    }
}

double Precice::ComputeMaxTimeStep(double computedTimestepLength) {

    double max_timestep;

    max_timestep = solverInterface.advance(computedTimestepLength);

    return max_timestep;
}

bool Precice::isCouplingOngoing() {
    return solverInterface.isCouplingOngoing();
}

bool Precice::isTimeWindowComplete() {
    return solverInterface.isTimeWindowComplete();
}

bool Precice::isActionRequired(const string &action) {
    return solverInterface.isActionRequired(action);
}

void Precice::markActionFulfilled(const std::string &action) {
    solverInterface.markActionFulfilled(action);
}

const string &Precice::getCowic() {
    return cowic;
}

const string &Precice::getCoric() {
    return coric;
}

void Precice::saveOldState(bool *StopCalc, double *dt) {

    CVariable *flow_nodes = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetNodes();
    CVariable *mesh_nodes = solver_container[ZONE_0][MESH_0][MESH_SOL]->GetNodes();
    auto geometry_nodes = geometry_container[ZONE_0][MESH_0]->nodes;
    auto &GridVelGrad = geometry_nodes->GetGridVel_Grad();

    for (int iPoint = 0; iPoint < nPoint; iPoint++) {

        //Save displacement at the last time step
        for (int iDim = 0; iDim < nDim; iDim++) {
            Displacement_Saved[iPoint][iDim] = mesh_nodes->GetBound_Disp(iPoint, iDim);
        }

        //Save the flow solution
        solution_Saved[iPoint] = flow_nodes->GetSolution(iPoint);
        solution_time_n_Saved[iPoint] = flow_nodes->GetSolution_time_n(iPoint);
        solution_time_n1_Saved[iPoint] = flow_nodes->GetSolution_time_n1(iPoint);

        //Save grid velocity
        GridVel_Saved[iPoint] = geometry_nodes->GetGridVel(iPoint);

        for (int iDim = 0; iDim < nDim; iDim++) {
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

void Precice::reloadOldState(bool *StopCalc, double *dt, bool initialTimeStep) {

    CVariable *flow_nodes = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetNodes();
    CVariable *mesh_nodes = solver_container[ZONE_0][MESH_0][MESH_SOL]->GetNodes();
    auto geometry_nodes = geometry_container[ZONE_0][MESH_0]->nodes;
    auto &GridVelGrad = geometry_nodes->GetGridVel_Grad();

    for (int iPoint = 0; iPoint < nPoint; iPoint++) {

        //Reload the displacement at the last time step
        mesh_nodes->SetBound_Disp(iPoint, Displacement_Saved[iPoint]);

        if (!initialTimeStep) {
            //Reload the flow solution at the last time step
            flow_nodes->SetSolution(iPoint, solution_time_n1_Saved[iPoint]);
            flow_nodes->Set_Solution_time_n(iPoint, solution_time_n1_Saved[iPoint]);
            flow_nodes->Set_Solution_time_n1(iPoint, solution_time_n1_Saved[iPoint]);

            //Reload grid velocity
            for (int iDim = 0; iDim < nDim; iDim++) {
                geometry_nodes->SetGridVel(iPoint, iDim, GridVel_Saved[iPoint][iDim]);
            }

            //Reload grid velocity gradient
            for (int iDim = 0; iDim < nDim; iDim++) {
                for (int jDim = 0; jDim < nDim; jDim++) {
                    GridVelGrad[iPoint][iDim][jDim] = GridVel_Grad_Saved[iPoint][iDim][jDim];
                }
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

void Precice::finalize() {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize - 1 << ": Finalizing preCICE..." << endl;
    solverInterface.finalize();
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize - 1 << ": Done finalizing preCICE!" << endl;
}
