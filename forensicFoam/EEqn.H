{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
//      + fvc::ddt(rho, K) + fvc::div(phi, K)
//      + (
//            he.name() == "e"
 //         ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
 //         : -dpdt
  //      )
      + thermophysicalTransport->divq(he)
     ==
    //    reaction->Qdot()
       fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info << "==================== Some Information " << endl;
    Info<< "min/max(T) = "   << min(T).value()   << ", " << max(T).value() << endl;
    Info<< "min/max(rho) = " << min(rho).value() << ", " << max(rho).value() << endl;

    // integral 
    int spec = 6;  //C2H2 index 6   H2 index 5 
    volScalarField& H2 = composition.Y(spec);

    Info<< " max(H2) = "  << max(H2).value() << endl; 
    const scalarField& V = mesh.V();
    float sumH2=0;
    forAll(H2,cellI)
    {
     sumH2 += H2[cellI]*V[cellI]*rho[cellI]; // H2 mass
    }
    //Info << " mass (core) " << sumH2 << endl;
    reduce(sumH2, sumOp<scalar>());
    Info << " mass= " << sumH2 << endl;
    Info << "===================================== " << endl;

}
