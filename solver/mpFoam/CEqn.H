{
    Info<< "Compute concentration field ..." << endl;
    volScalarField& diffEff1 = diffCo1;
    const volScalarField& liquidPhase = alpha1;

    //- update diffusion coefficient basing on concentration and ion cloud
    
    if(ionEffect)
    {
        Info<< "Correct diffusion coefficients for ionic strength ..." << endl;
        forAll(mesh.C(),i)
        {
            if(C1[i]<=0.0)
            {
                diffEff1[i] = 0;
            }
            else
            {
                scalar CStar = Foam::sqrt(C1[i]/conSat);
                diffEff1[i] = diffInLiquid1.value()*(1.0-(ionEffectA*CStar/2.0/pow((1+ionEffectB*CStar),2)))*Foam::pow(liquidPhase[i],2);
            }
            if(diffEff1[i]<1e-10)
            {
                diffEff1[i]=0.0;
            }
        }
    }
    else
    {
        forAll(mesh.C(),i)
        {
            diffEff1[i] = diffInLiquid1.value()*pow(liquidPhase[i],2);
            if(diffEff1[i]<1e-10)
            {
                diffEff1[i]=0.0;
            }
        }
    }
 
    Info<< "min/max(diffEff1) = " << min(diffEff1).value() << ", " << max(diffEff1).value() << endl;

    fvScalarMatrix CEqn1
    (
        fvm::ddt(liquidPhase,C1)
      + fvm::div(phi, C1, "div(phi,C)")
      - fvm::laplacian(diffEff1, C1,  "laplacian(diffEff,C)")
      ==
        phaseChange1.preRect(C1, SuOut1)
      //- phaseChange.Kexp(C)/molW
      + fvOptions(C1)
    );

    CEqn1.relax();

    fvOptions.constrain(CEqn1);

    CEqn1.solve();

    // constrain C1 to prevent stability issue
    // concentration must be positive
    const dimensionedScalar CZero("CZero", dimMoles/dimVolume, 0.0);
    C1 = max(CZero, C1);
    // dissolution reaction must be lower than saturation concentration
    if(precDissSolid1 == "dissolution")
    {
        C1 = min(C1, Csat1);
    }
    else if(precDissSolid1 == "precipitation")
    {
        C1 = max(Csat1, min(C1, Cinit1));
    }

    Info<< "min/max(C1) = "
        << min(C1).value() << ", " << max(C1).value() << endl;

    volScalarField& diffEff2 = diffCo2;

    forAll(mesh.C(),i)
    {
        diffEff2[i] = diffInLiquid2.value()*pow(liquidPhase[i],2);
        if(diffEff2[i]<1e-10)
        {
            diffEff2[i]=0.0;
        }
    }

    Info<< "min/max(diffEff2) = " << min(diffEff2).value() << ", " << max(diffEff2).value() << endl;

    fvScalarMatrix CEqn2
    (
        fvm::ddt(liquidPhase,C2)
      + fvm::div(phi, C2, "div(phi,C)")
      - fvm::laplacian(diffEff2, C2,  "laplacian(diffEff,C)")
      ==
        phaseChange2.preRect(C2, SuOut2)
      //- phaseChange.Kexp(C)/molW
      + fvOptions(C2)
    );

    CEqn2.relax();

    fvOptions.constrain(CEqn2);

    CEqn2.solve();

    // constrain C2 to prevent stability issue
    // concentration must be positive
    C2 = max(CZero, C2);
    // dissolution reaction must be lower than saturation concentration
    if(precDissSolid2 == "dissolution")
    {
        C2 = min(C2, Csat2);
    }
    else if(precDissSolid2 == "precipitation")
    {
        C2 = max(Csat2, min(C2, Cinit2));
    }

    Info<< "min/max(C2) = "
        << min(C2).value() << ", " << max(C2).value() << endl;

    Info<< "Concentration calculation complete!!!" << endl;

}
