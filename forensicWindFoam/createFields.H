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


//manually adjust number of steps
windsteps = 24;  //<------------------


// READ DATA 
std::ifstream fin("wind.dat"); // Open the file
string line;
    
if (!fin.is_open()) {
    std::cerr << "Error opening Wind file." << std::endl;
    return 1;
}

const int ncolumns =4;
float c[ncolumns];
int i=0;

// discard first line (header)
getline(fin,line);
// read data
while (fin >> c[0] >> c[1] >> c[2] >> c[3])
{
wind_time[i]  = c[0];
wind_speed[i] = c[1];
wind_ang[i]   = c[2];
wind_T[i]     = c[3];
i++;
}
fin.close();

// check if data is read correctly
//for (int i = 0; i<windsteps; i++) 
//{
//Info << i << " speed= " <<  wind_speed[i] << " ang= " << wind_ang[i] << " T= " << wind_T[i] << endl;
//}


// change from h to secs
for (int i = 0; i<windsteps; i++) { wind_time[i] = wind_time[i]*3600.0;}  // **  DEBUG 3600


// assign current time
wind_speedT = wind_speed[currentwindstep];
wind_angT   = wind_ang[currentwindstep];
wind_TT     = wind_T[currentwindstep];     


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
