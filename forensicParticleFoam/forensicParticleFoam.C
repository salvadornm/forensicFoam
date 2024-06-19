/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoParticleFoam

Description
    Transient solver for the passive transport of a particle cloud.

    Uses a pre- calculated velocity field to evolve the cloud.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "parcelCloudList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "compressibleCourantNo.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        clouds.storeGlobalPositions();

        mesh.update();

        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Move the mesh
        mesh.move();

        if (mesh.changing())
        {
            U.correctBoundaryConditions();
        }
        
        scalar t = runTime.value();
        scalar dt = runTime.deltaTValue();
       
        Info<< "[SNM] time = " << t <<  " deltat=" << dt << endl;

        // check time

        if (t >= Time2)
        {
            indexTime1++;indexTime2++;
            Time1 = timevecDirs[indexTime1];
            Time2 = timevecDirs[indexTime2];
            word nuevotiempo =  timeDirs[indexTime2].name();

            // Info << " U2-> U1" << endl;
            forAll(mesh.C(),cellI)
            {
            U1[cellI] = U2[cellI]; 
            }
            Info << " Read New U2 time = " << nuevotiempo << endl;


            // Check Limit of files
            if (indexTime2 > nDirs)
            {
                Info << " [forensicPartcileFoam]: " << endl;
                Info << " No more files to read \n " << endl;
                return 0;
            }


            // Dynamically allocate a new vector
            volVectorField* U2p = new volVectorField
            (
                IOobject
                (
                "U",
                nuevotiempo,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
                ),
                mesh
            );            
            // U2p -> U2 
            forAll(mesh.C(),cellI)
            {
            U2[cellI] = (*U2p)[cellI]; 
            }
            // delete vector to free memory 
            delete U2p;
            
        }

        // interpolation weights
        scalar w2 = (t -Time1)/(Time2-Time1);
        scalar w1 = 1.0 - w2; 

        // debug lines
        // Info<< " ----  " << endl;
        // Info<< " w1= " << w1 << " and  w2=" << w2 << endl;
        // Info << " i1 and i2 = " << indexTime1 << " "  << indexTime2  << endl;
        // Info << " t1 and t2 = " << Time1 << " "  << Time2 << endl;
        
        // Create U by interpolation in time        
        forAll(mesh.C(),cellI)
        {
           U[cellI] = w1*U1[cellI] + w2*U2[cellI];
        }

        // debug
        Info << "U1 = " << U1[icellexample] << endl; 
        Info << "U  = " << U[icellexample] << endl; 
        Info << "U2 = " << U2[icellexample] << endl; 

        scalar Zref = 1800;
        // Info<< " Adjust top BC at z >  " << Zref << endl;
        vector Uw = vector(0.0, 0.0, 1.0);        
        forAll(mesh.C(),cellI)
        {
           const scalar z = mesh.C()[cellI].z();
           if (z > Zref) 
            { U[cellI] = Uw;}
        }

        // ---------------------------------------------------------//

        clouds.evolve();

        
        runTime.write();

        // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        //     << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        //     << nl << endl;
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
