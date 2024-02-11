// Solve the Momentum equation

forAll(mesh.C(),cellI)
{
    Fg[cellI][0]  = 0.0;
    Fg[cellI][1]  = 0.0;
    Fg[cellI][2]  = -9.81*rho[cellI];
}


MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
 ==
    fvModels.source(rho, U)
 // + Fg //SNM  
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
