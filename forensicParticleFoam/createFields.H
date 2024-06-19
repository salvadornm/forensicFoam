#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Creating U1 and U2 \n" << endl;

volVectorField U1
(
    IOobject
    (
        "U1",
        runTime.timeName(),
        mesh,        
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength/dimTime, Zero) // dimensions and Initialise 0
);

// volVectorField U2
// (
//     IOobject
//     (
//         "U2",
//         runTime.timeName(),
//         mesh,        
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedVector(dimLength/dimTime, Zero) // dimensions and Initialise 0
// );

Info<< "Calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho*U) & mesh.Sf()
);


// Info<< "Creating turbulence model\n" << endl;
// autoPtr<compressible::momentumTransportModel> turbulence
// (
//     compressible::momentumTransportModel::New
//     (
//         rho,
//         U,
//         phi,
//         thermo
//     )
// );

Info<< "Constructing clouds" << endl;
parcelCloudList clouds(rho, U, thermo.mu(), g);


//  U->U1

forAll(mesh.C(),cellI)
{
U1[cellI] = U[cellI]; 
}


Info << " ********* SNM ... \n";

// timeDirs is of type instantlist with the diredtories

instantList timeDirs = runTime.times();


//Info <<  timeDirs  << endl;

// Calculate ndirectories and allocate 
const int  nDirs = timeDirs.size() ;
scalar timevecDirs[nDirs+1];

Info << " nDirs = " << nDirs  << endl;

// convert list into an array
int timeindex = 1;
for (label timei=1; timei < timeDirs.size(); ++timei)
{
    timevecDirs[timeindex] = mag( timeDirs[timei].value() );
    timeindex++;
}

// for (int i=1;i < nDirs;i++)
// {
//    Info << " [debug] timevecDirs= " << timevecDirs[i] << endl; 
// }

// global index time
int indexTime1,indexTime2;
scalar Time1,Time2;

// find closest index lower than starting time
scalar t0 = runTime.value();

int ifound = -1;
for (int i=1;i < nDirs;i++)
{
   scalar diff = timevecDirs[i] -t0 + 1e-8;
   if (diff > 0)
   {
    ifound = i;
    break;
   }
}

Info << " Ifound = " << ifound << endl;

indexTime1 = ifound; Time1 = timevecDirs[indexTime1];
indexTime2 = indexTime1 + 1; Time2 =  timevecDirs[indexTime2];

Info << " i1 and i2 = " << indexTime1 << " "  << indexTime2  << endl;
Info << " t1 and t2 = " << Time1 << " "  << Time2 << endl;

Info << " Initial time===" << t0  << endl;


const int icellexample = 23;


Info << "U1 = " << U1[icellexample] << endl; 

word nuevotiempo =  timeDirs[indexTime2].name();

Info << " tiempo " << nuevotiempo << endl;

//  Create New Vector U2 and read it form directory at time 2
// volVectorField* U2p = new volVectorField
// (
//     IOobject
//     (
//         "U",        
//         nuevotiempo,
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );

volVectorField U2
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


//Info << "U2 = " << U2[icellexample] << endl; 
//Info << "U2p = " << (*U2p)[icellexample] << endl; 


// Info<< "Exit SNM\n" << endl;
// return 0;
