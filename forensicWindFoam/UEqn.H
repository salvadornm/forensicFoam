// Solve the Momentum equation




// WIND DATA
const scalar Vref = 3.12; scalar windangle = 90;// DATA FROM FILE...

// constants for ABL
const scalar zmax  = 500;
const scalar zg = 0; const scalar z0 = 0.5; const scalar kappa =0.41;
const scalar Zref = 100.0;  
const scalar Ustar =  kappa*Vref/log((Zref+z0)/z0);
const scalar deg2rad = 3.141592/90;
windangle = windangle*deg2rad;
// nudging miodel comstant
const scalar otauf =1.0/200.0; // tau=200 secs

// Nudging model
// forAll(mesh.C(),cellI)
// {
//     const scalar z = mesh.C()[cellI].z(); 
//     scalar zbl = min(z,zmax) - zg + z0;
//     scalar Uwind = Ustar*log(zbl/z0)/kappa;
//     scalar up = Uwind*sin(windangle);
//     scalar vp = Uwind*cos(windangle);
//     //scalar up = Vref*min(pow(z/Zref, 2),1); 
//     Fwind[cellI][0]  = (up-U[cellI][0])*otauf;  // use conditional Uxy
//     Fwind[cellI][1]  = (vp-U[cellI][1])*otauf;  
//     Fwind[cellI][2]  = 0.0;
// }


MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
 ==
    fvModels.source(rho, U)
 // + fvc::Sp(rho,Fwind)   //SNM
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
