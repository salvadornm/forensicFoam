//  Compute planar averages over height Uavg(z)         // 
//  NOTE: The average is store in an arrays   Uavg(nzprof)   //
//  NOTE: Uavg(k) is a vector  Uavg[k][0]=Uavg[i].component(vector::X)

// reset arrays
for (int k = 0; k<nzprof; k++)
{
    Uavg[k][0] = 0.0;
    Uavg[k][1] = 0.0;
    Uavg[k][2] = 0.0;
    Vavg[k] = 0.0;
}

// loop over cells and store sum of velocities
const auto& volume = mesh.V();

forAll(mesh.C(),cellI)
{
    const scalar z = mesh.C()[cellI].z(); 
    
    // correspond k-index to z coordinate, NOTE: int roudns down  7.99 ~ 7
    int k = int(z*nzprof/Lz);

    Uavg[k][0]  += U[cellI][0]*volume[cellI];
    Uavg[k][1]  += U[cellI][1]*volume[cellI];
    Uavg[k][2]  += U[cellI][2]*volume[cellI];
    Vavg[k]     += volume[cellI];
    
}

vector umax,umin ;
umax[0] = -1e20;umax[1]= -1e20;umax[2]= -1e20;
umin[0] = +1e20;umin[1]= +1e20;umin[2]= +1e20;


// use OpenFOAM sum to sum over all cores  and compute averages
for (int k = 0; k<nzprof; k++)
{
    // global sums
    vector aux  = returnReduce(Uavg[k], sumOp<vector>());
    scalar auxV = returnReduce(Vavg[k], sumOp<scalar>()) + 0.00001;
    // average
    Uavg[k] = aux /auxV;
    // for info
    umax = max(Uavg[k],umax);
    umin = min(Uavg[k],umin);

    Uavg[k][2] = 0;  // z-component force to be 0

    // temp [DBG]
    //Info << k << " UAVG = " <<  Uavg[k]  << endl;
}


Info << "min/max  Uavg = " << umin << ", " <<  umax << endl; 
Info << "====================================== " << endl;

Info <<  " Lmax=" << Lmin  << " Lmax=" << Lmax  << endl;



//Info << "min/max  XYZ  =" << min(mesh.C()).value() << ", " <<  max(mesh.C()).value()  << endl;
//Info<< " min/max  U = "   << min(U).value()   << ", " << max(U).value()   << endl;

