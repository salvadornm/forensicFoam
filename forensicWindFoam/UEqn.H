//  ...... Solve the Momentum equation    // 

// -----      DATA FROM FILE                    ----//
const scalar Vref = 3.124100; scalar windangle = 140.194430;
// -

// constants for ABL
const scalar zmax  = 500;
const scalar zg = 0; const scalar z0 = 0.5; const scalar kappa =0.41;
const scalar Zref = 100.0;  
const scalar Ustar =  kappa*Vref/std::log((Zref+z0)/z0);
const scalar deg2rad = 0.5*3.141592/90;
windangle = windangle*deg2rad;
// nudging model comstant
const scalar otauf =1.0/200.0; // tau=200 secs
scalar udir = std::sin(windangle);
scalar vdir = std::cos(windangle);

// Nudging model
forAll(mesh.C(),cellI)
{
    const scalar z = mesh.C()[cellI].z(); 
    scalar zbl = min(z,zmax) - zg + z0;
    scalar Uwind = Ustar*std::log(zbl/z0)/kappa;
    scalar up = Uwind*udir;
    scalar vp = Uwind*vdir;

    // correspond k-index to z coordinate
    int k = int(z*nzprof/Lz);
    vector uforce = Uavg[k]; // x-y average

    //vector uforce = U[cellI]; // DBG 

    Fwind[cellI][0]  = (up-uforce[0])*otauf;  
    Fwind[cellI][1]  = (vp-uforce[1])*otauf;  
    Fwind[cellI][2]  = 0.0;
}


MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
 ==
    fvModels.source(rho, U)
    + fvc::Sp(rho,Fwind)   // Nudging model 
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rho)  //SNM  coment this line 
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
