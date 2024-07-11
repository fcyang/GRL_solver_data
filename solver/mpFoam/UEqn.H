    MRF.correctBoundaryVelocity(U);

    const volScalarField& liquidPhase = alpha1;

    fvVectorMatrix UEqn
    (
        1./(liquidPhase+SMALL)*fvm::ddt(rho, U) + fvm::div((rhoPhi/(fvc::interpolate(liquidPhase)+SMALL)), U)
      //+ MRF.DDt(rho, U)
      //+ turbulence->divDevRhoReff(rho, U)
      - fvm::laplacian((rho*mixture.nuModel1().nu())/(liquidPhase+SMALL), U)
     ==
        fvOptions(rho, U)
    );

    phaseChange1.addInterfacePorosity(UEqn);

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
