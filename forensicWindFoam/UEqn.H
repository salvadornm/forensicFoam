//  ...... Solve the Momentum equation    // 

// -----      DATA FROM FILE                    ----//
//const scalar Vref = 3.124100; scalar windangle = 140.194430;
const scalar Vref = wind_speedT; scalar windangle = wind_angT;

// -

// constants for ABL
const scalar zmax  = 500;
const scalar zg = 0; const scalar z0 = 0.5; const scalar kappa =0.41;
const scalar Zref = 100.0;  
const scalar Ustar =  kappa*Vref/std::log((Zref+z0)/z0);

windangle = windangle*deg2rad;
// nudging model comstant
const scalar otauf    = 1.0/10.0;  // tau=10 secs
const scalar otaufbuf = 1.0/5.0;   // tau=3 secs (buffer)

scalar udir = std::sin(windangle);
scalar vdir = std::cos(windangle);

if (nudging==1)
{
// Nudging model
forAll(mesh.C(),cellI)
{
    const scalar x = mesh.C()[cellI].x(); 
    const scalar y = mesh.C()[cellI].y(); 
    const scalar z = mesh.C()[cellI].z(); 

    const scalar damp = buffer[cellI];  // buffer area

    scalar zbl = min(z,zmax) - zg + z0;
    scalar Uwind = Ustar*std::log(zbl/z0)/kappa;
    vector Uw;
    Uw[0] = Uwind*udir;
    Uw[1] = Uwind*vdir;
    Uw[2] = 0.0;

    // correspond k-index to z coordinate
    int k = int(z*nzprof/Lz);

    // forcing term
    Fwind[cellI]  = (Uw -Uavg[k])*otauf*(1.0-damp);  

    // In BUFFER AREA NEAR BOUNDARIES in X-Y
    // force more
    Fwind[cellI] += (Uw -U[cellI])*otaufbuf*damp;                
    // add  small randomness (10% intensity) in buffer
    scalar auxr=auxrand.scalarNormal(); 
    Fwind[cellI] += auxr*Uw*0.1*damp*otauf;// dived by dt
}
//
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
