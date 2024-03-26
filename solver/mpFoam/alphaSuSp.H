volScalarField Su1 = SuOut1/rho2; // check to see which rho to use for different solid/liquid densities
volScalarField Su2 = SuOut2/rho3;
//Su += fvc::div(phi)*min(max(alpha1, scalar(0)), scalar(1));
Info << "Su1 max/min: " << max(Su1).value() << ", " << min(Su1).value() << endl;
Info << "Su2 max/min: " << max(Su2).value() << ", " << min(Su2).value() << endl;
//zeroField Su;

zeroField Sp;

zeroField divU;
