{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh,
                mesh.ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if (nAlphaSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Add the optional shear compression contribution
    if (scAlpha > 0)
    {
        phic +=
            scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }


    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
    }

    if (true)//MULESCorr)
    {
        #include "alphaSuSp.H"

        Info<< "Time step: " << mesh.time().deltaTValue() << endl;
        forAll(mesh.C(), cellI)
        {
            alpha2[cellI] = alpha2[cellI]+Su1[cellI]*mesh.time().deltaTValue();
            alpha3[cellI] = alpha3[cellI]+Su2[cellI]*mesh.time().deltaTValue();
        }

        //fvScalarMatrix alpha1Eqn
        //(
        //    (
        //        LTS
        //      ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
        //      : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
        //    )
        //  + fv::gaussConvectionScheme<scalar>
        //    (
        //        mesh,
        //        phiCN,
        //        upwind<scalar>(mesh, phiCN)
        //    ).fvmDiv(phiCN, alpha1)
       //// - fvm::Sp(fvc::ddt(dimensionedScalar("1", dimless, 1), mesh)
       ////           + fvc::div(phiCN), alpha1)
        // ==
        //    Su + fvm::Sp(Sp + divU, alpha1)
        //);

        //alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        Info<< "Phase-2 volume fraction = "
            << alpha2.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
            << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
            << endl;
        
        Info<< "Phase-3 volume fraction = "
            << alpha3.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha3.name() << ") = " << min(alpha3).value()
            << "  Max(" << alpha3.name() << ") = " << max(alpha3).value()
            << endl;

        //tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        //alphaPhi10 = talphaPhi1UD();

        //if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        //{
        //    Info<< "Applying the previous iteration compression flux" << endl;
        //    MULES::correct
        //    (
        //        geometricOneField(),
        //        alpha1,
        //        alphaPhi10,
        //        talphaPhi1Corr0.ref(),
        //        oneField(),
        //        zeroField()
        //    );

        //    alphaPhi10 += talphaPhi1Corr0();
        //}

        //// Cache the upwind-flux
        //talphaPhi1Corr0 = talphaPhi1UD;

        // remove any solid phase error outside of interfacial region
        Info<< "Update Cmask ... " << endl;
        phaseChange1.updateCmask();

        // generating nucleation site for empty cell next to wall BC
        forAll(mesh.C(), cellI)
        {
            if((nuSite[cellI]==1.0) && (alpha1[cellI]>nuIntVofCut))
            {
                alpha2[cellI] += nuIntVof*mesh.time().deltaTValue();
                Cmask[cellI] = 1.0;
            }
        }

        // forAll(mesh.C(),cellI)
        // {
        //     if(Cmask[cellI]==0.0)
        //     {
        //         alpha1[cellI]=1.0;
        //     }
        // }   

        alpha1 = 1.0 - alpha2 - alpha3;

        mixture.correct();
    }


    //for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    //{
    //    #include "alphaSuSp.H"

    //    Info<< "Time step: " << mesh.time().deltaTValue() << endl;
    //    forAll(mesh.C(), cellI)
    //    {
    //        alpha1[cellI] = alpha1[cellI]+Su[cellI]*mesh.time().deltaTValue();
    //    }

    //    //surfaceScalarField phir(phic*mixture.nHatf());

    //    //tmp<surfaceScalarField> talphaPhi1Un
    //    //(
    //    //    fvc::flux
    //    //    (
    //    //        phiCN(),
    //    //        cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime(),
    //    //        alphaScheme
    //    //    )
    //    //  + fvc::flux
    //    //    (
    //    //       -fvc::flux(-phir, alpha2, alpharScheme),
    //    //        alpha1,
    //    //        alpharScheme
    //    //    )
    //    //);

    //    //if (MULESCorr)
    //    //{
    //    //    tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);
    //    //    volScalarField alpha10("alpha10", alpha1);

    //    //    MULES::correct
    //    //    (
    //    //        geometricOneField(),
    //    //        alpha1,
    //    //        talphaPhi1Un(),
    //    //        talphaPhi1Corr.ref(),
    //    //        Sp,
    //    //        (-Sp*alpha1)(),
    //    //        oneField(),
    //    //        zeroField()
    //    //    );

    //    //    // Under-relax the correction for all but the 1st corrector
    //    //    if (aCorr == 0)
    //    //    {
    //    //        alphaPhi10 += talphaPhi1Corr();
    //    //    }
    //    //    else
    //    //    {
    //    //        alpha1 = 0.5*alpha1 + 0.5*alpha10;
    //    //        alphaPhi10 += 0.5*talphaPhi1Corr();
    //    //    }
    //    //}
    //    //else
    //    //{
    //    //    alphaPhi10 = talphaPhi1Un;

    //    //    MULES::explicitSolve
    //    //    (
    //    //        geometricOneField(),
    //    //        alpha1,
    //    //        phiCN,
    //    //        alphaPhi10,
    //    //        Sp,
    //    //        (Su + divU*min(alpha1(), scalar(1)))(),
    //    //        oneField(),
    //    //        zeroField()
    //    //    );
    //    //}

    //    // remove any solid phase error outside of interfacial region
    //    Info<< "Update Cmask ... " << endl;
    //    phaseChange.updateCmask();

    //    forAll(mesh.C(),cellI)
    //    {
    //        if(Cmask[cellI]==0.0)
    //        {
    //            alpha1[cellI]=1.0;
    //        }
    //    }   

    //    alpha2 = 1.0 - alpha1;

    //    mixture.correct();
    //}

//    if (alphaApplyPrevCorr && MULESCorr)
//    {
//        talphaPhi1Corr0 = alphaPhi10 - talphaPhi1Corr0;
//        talphaPhi1Corr0.ref().rename("alphaPhi1Corr0");
//    }
//    else
//    {
//        talphaPhi1Corr0.clear();
//    }
//
//    #include "rhofs.H"
//
//    if
//    (
//        word(mesh.ddtScheme("ddt(rho,U)"))
//     == fv::EulerDdtScheme<vector>::typeName
//     || word(mesh.ddtScheme("ddt(rho,U)"))
//     == fv::localEulerDdtScheme<vector>::typeName
//    )
//    {
//        rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
//    }
//    else
//    {
//        if (ocCoeff > 0)
//        {
//            // Calculate the end-of-time-step alpha flux
//            alphaPhi10 =
//                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
//        }
//
//        // Calculate the end-of-time-step mass flux
//        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
//    }

    liquidMassOld_ = liquidMass_;
    liquidMass_ = gSum((mesh.V()*rho1*alpha1)());
    solidMass_ = gSum((mesh.V()*rho2*alpha2 + mesh.V()*rho3*alpha3)());
    dmdtTot_ = gSum((mesh.V()*SuOut1+mesh.V()*SuOut2)());
    deltaT_ = mesh.time().deltaTValue();

    Info<< "Total liquid phase weight = "
        << liquidMass_.value()
        << endl;

    Info<< "Total solid phase weight = "
        << solidMass_.value()
        << endl;

    Info<< "Total dmdt = "
        << dmdtTot_.value()
        << endl;

    if(Pstream::master() && runTime.outputTime())
    {
        string fileName = "outputResults.csv";
        std::ifstream file(fileName);
        std::fstream fout;
        fout.open(fileName, std::ios::out | std::ios::app);
        if (file.peek() == std::ifstream::traits_type::eof())
        {
            fout << "Time[s],liquidMass[kg],solidMass[kg],dmdt[kg/s],deltaT[s]" << "\n";
        }
        fout << mesh.time().timeName() << ", "
            << liquidMass_.value() << ", "
            << solidMass_.value() << ", "
            << dmdtTot_.value() << ", "
            << runTime.deltaT().value() << "\n";
    }

    // constrain liquid VOF to [0, 1]
    forAll(mesh.C(),i)
    {
        if(alpha1[i]<0.0)
        {
            alpha1[i]=0.0;
        }
        else if(alpha1[i]>1.0)
        {
            alpha1[i]=1.0;
            alpha2[i]=0.0;
            alpha3[i]=0.0;
        }
        else if(alpha2[i]>1.0)
        {
            alpha1[i]=0.0;
            alpha2[i]=1.0;
            alpha3[i]=0.0;
        }
        else if(alpha3[i]>1.0)
        {
            alpha1[i]=0.0;
            alpha2[i]=0.0;
            alpha3[i]=1.0;
        }

        //if(alpha1[i]>fluidOnly_)
        //{
        //    alpha1[i]=1.0;
        //    alpha2[i]=0.0;
        //    alpha3[i]=0.0;
        //}
    }
}


