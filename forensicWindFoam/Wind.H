//  Adjusting Wind  by interpolating dat       // 
Info << "====================================== " << endl;
Info << " Adjusting Wind " << endl;

scalar t=runTime.value();   //current time

// advance time
if (t > wind_time[currentwindstep + 1]) {currentwindstep++;}

scalar dt  = t -wind_time[currentwindstep];
scalar dt0 = wind_time[currentwindstep + 1]-  wind_time[currentwindstep];

scalar dwdt = (wind_speed[currentwindstep + 1] - wind_speed[currentwindstep])/dt0;
scalar dadt = (wind_ang[currentwindstep + 1]   -   wind_ang[currentwindstep])/dt0;
scalar dTdt = (wind_T[currentwindstep + 1]     -     wind_T[currentwindstep])/dt0;


wind_speedT = wind_speed[currentwindstep] + dwdt*dt;
wind_angT   = wind_ang[currentwindstep]   + dadt*dt;
wind_TT     = wind_T[currentwindstep]     + dTdt*dt;

Info << " currentwindstep  " << currentwindstep << " time=" << t << endl;
Info << " Wind  Speed= " << wind_speedT << " Angle=  " << wind_angT << " T= " << wind_TT << endl;
Info << "====================================== " << endl;
