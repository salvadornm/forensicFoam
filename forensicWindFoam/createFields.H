#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidReactionThermo> pThermo(fluidReactionThermo::New(mesh));
fluidReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();


#include "compressibleCreatePhi.H"

mesh.schemes().setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidReactionThermophysicalTransportModel>
thermophysicalTransport
(
    fluidReactionThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> reaction(combustionModel::New(thermo, turbulence()));

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    thermo.incompressible()
);

mesh.schemes().setFluxRequired(p_rgh.name());

hydrostaticInitialisation
(
    p_rgh,
    p,
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);


// SNM: auxiliar vector
volVectorField Fg
(
    IOobject
    (
        "Fg",
        runTime.timeName(),
        mesh,        
        IOobject::MUST_READ
    ),
    mesh
);

// SNM: Global Wind restitutive force (overall pressure gradient)
Info<< "Creating field Fwind\n" << endl;
volVectorField Fwind
(
    IOobject
    (
        "Fwind",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength/dimTime/dimTime, Zero)
);

Info<< "Creating field buffer\n" << endl;
volScalarField buffer
(
    IOobject
    (
        "buffer",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength/dimLength, Zero) // no dimensions
);


// arrays to build the plane averages
const int nzprof=400;
vectorField Uavg(nzprof);    // z-profile         Uavg(z)
scalarField Vavg(nzprof);    // aux array storing sumV(z)
scalarField Tavg(nzprof);    // aux array storing Tavg(z)
const scalar Lz=max(mesh.C()).value()[2];  // size of domain in Z
const vector Lmax = max(mesh.C()).value();
const vector Lmin = min(mesh.C()).value();
const scalar XBUF = 500;const scalar YBUF = 400;const scalar ZBUF = 100;
Random auxrand(1);

Info<< " Reading Wind Data \n " <<  endl;
// arrays to store wind data
const int nmaxwindsteps=30; // max number of data stored
float wind_speed[nmaxwindsteps],wind_ang[nmaxwindsteps],wind_T[nmaxwindsteps];
float wind_time[nmaxwindsteps];
float wind_speedT,wind_angT,wind_TT;
int windsteps=0;   // number of wind steps in file (to read in file)
int currentwindstep = 0; // actual  wind step
int variable_wind = 0; 

// READ DATA

//manual adjust to test
windsteps = 6;

wind_time[0] = 0; wind_speed[0]=  3.124100; wind_ang[0]=  140.194430;wind_T[0]= 309.787196;
wind_time[1] = 1; wind_speed[1]=  3.929377; wind_ang[1]=  165.256440;wind_T[1]= 308.007197;
wind_time[2] = 2; wind_speed[2]=  4.909175; wind_ang[2]=  176.496480;wind_T[2]= 307.237196;
wind_time[3] = 3; wind_speed[3]=  4.925444; wind_ang[3]=  185.826340;wind_T[3]= 305.457198;
wind_time[4] = 4; wind_speed[4]=  4.780168; wind_ang[4]=  195.780750;wind_T[4]= 304.877198;
wind_time[5] = 5; wind_speed[5]=  3.883298;	wind_ang[5]=  214.508530;wind_T[5]=	304.497198;

//
wind_speedT = wind_speed[currentwindstep];
wind_angT   = wind_ang[currentwindstep];
wind_TT     = wind_T[currentwindstep];     

// std::ifstream file("wind.dat"); // Open the file
    
// if (!file.is_open()) {
//     std::cerr << "Error opening Wind file." << std::endl;
//     return 1;
// }

// float c1, c2, c3, c4;
    
// while (file >> c1 >> c2 >> c3 >> c4)
// { // Read until end of file
//    // Process the read values
//     std::cout << "C 1: " << c1 << ", C 2: " << c2 << ", C 3: " << c3 << " C 4: " << c4 << std::endl;
// }
// file.close(); // Close the file


// change from h to secs
for (int i = 0; i<windsteps; i++) { wind_time[i] = wind_time[i]*3600.0;}  // **  DEBUG 3600


// buffer array

forAll(mesh.C(),cellI)
{
    const scalar x = mesh.C()[cellI].x(); 
    const scalar y = mesh.C()[cellI].y(); 
    const scalar z = mesh.C()[cellI].z(); 

    buffer[cellI] = 0;
    if (((x - Lmin[0]) < XBUF) || ((Lmax[0]- x) < XBUF)  ||
        ((y - Lmin[1]) < YBUF) || ((Lmax[1]- y) < YBUF) ) 
    {
       if (z < ZBUF) {buffer[cellI] = 1;}
    }    

}


//********************************

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
