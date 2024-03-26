/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseChangeReaction.H"
#include "surfaceDensity.H"
#include "surfaceReaction.H"
#include "solidNucleation.H"
#include "constants.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvcReconstruct.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"
#include <random>

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeReaction::phaseChangeReaction
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarField& C,
    const volScalarField& alpha,
    const volScalarField& solid,
    volScalarField& Cmask
)
:
    dict_(dict),
    C_(C),
    mesh_(mesh),
    solid_(solid),
    phaseName_(dict.lookup("phaseName")),
    dirGrowthFlag_(dict.lookupOrDefault<bool>("dirGrowth", false)),
    dirGrowthNum_(dict.lookupOrDefault<scalar>("dirGrowthNum", 1.0)),
    dirGrowthNorm_(dict.lookupOrDefault<List<vector>>("dirGrowthNorm", {vector::zero})),
    dirFactor_(dict.lookupOrDefault<List<scalar>>("dirFactor", {1.0})),
    Cu_(dict.lookupOrDefault<scalar>("Cu", 1e15)),
    pKsp_(dict.lookupOrDefault<scalar>("pKsp", 1.0)),
    lnA_(dict.lookupOrDefault<scalar>("lnA", 1.0)),
    B_(dict.lookupOrDefault<scalar>("B", 1.0)),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.0)),
    Vmol_(dict.lookupOrDefault<scalar>("Vmol", 1.0)),
    reacModel_(dict.lookupOrDefault<word>("reactionModel", "linear")),
    K_("K", dimMoles/dimArea/dimTime, dict),
    KList_(dict.lookupOrDefault<List<scalar>>("KList", {1.0})),
    Mv_("Mv", dimMass/dimMoles, dict),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 1.1)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0.01)),
    alphaSolidMin_(dict.lookupOrDefault<scalar>("alphaSolidMin", 0.0)),
    alphaRestMax_(dict.lookupOrDefault<scalar>("alphaRestMax", 0.01)),
    smoothSurface_(dict.lookupOrDefault<bool>("smoothSurface", false)),
    smoothAreaDensity_(dict.lookupOrDefault<scalar>("smoothAreaDensity", 1.0)),
    surfaceAreaDensity_(dict, mesh, alpha, solid),
    surfaceReacRate_(dict,C,solid,phaseName_),
    solidNucleation_(dict,C,alpha,solid),
    cornerCell_(dict.lookupOrDefault<bool>("cornerCell", false)),
    gradLim_(dict.lookupOrDefault<scalar>("gradientLimit", 1e20)),
    alpha_(alpha),
    Cmask_(Cmask),
    debug_(dict.lookupOrDefault<bool>("debug", false)),
    surfaceReaction_(dict_.lookupOrDefault("surfaceReaction",false)),
    nucleation_(dict.lookupOrDefault<bool>("nucleation", false)),
    nucleationBoxFlag_(dict.lookupOrDefault<bool>("nucleationBox", false)),
    nucMinPoint_(dict.lookupOrDefault<vector>("nucMinPoint", vector(0,0,0))),
    nucMaxPoint_(dict.lookupOrDefault<vector>("nucMaxPoint", vector(0,0,0)))
{
    // Convert from g/mol to Kg/mol
    Mv_.value() = Mv_.value() * 1e-3;
    info();
    Info<< "Constructor of phaseChangeReaction finished" << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeReaction::info()
{
    Info<< "Settings for phase " << phaseName_ << ":" << endl;
    Info<< "    " << dict_ << endl;
}

Foam::tmp<Foam::volScalarField> Foam::phaseChangeReaction::Kexp(const volScalarField& field)
{
    volScalarField from
    (
        "from",
        alpha_
    );
    from.correctBoundaryConditions();

    volScalarField to
    (
        "to",
        1.0-alpha_
    );
    to.correctBoundaryConditions();

    Info<< "alpha liquid max/min: "<< max(from).value() << ", " << min(from).value() << endl;
    Info<< "alpha solid max/min: "<< max(to).value() << ", " << min(to).value() << endl;
    const volVectorField gradFrom(fvc::grad(from));
    const volVectorField gradTo(fvc::grad(to));
    const volScalarField areaDensitySmooth
    (
        IOobject
        (
            "areaDensitySmooth",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("areaDensitySmoothInput", dimless/dimLength, smoothAreaDensity_)
    );

    //- Calculating surface area density for phase change reactions
    surfaceAreaDensity_.update();
    const volScalarField& areaDensity_ = surfaceAreaDensity_.areaDensity();
    
    const volScalarField gradAlphaf(gradFrom & gradTo);

    dimensionedScalar VSUnit
    (
        "VSUnit",
        dimensionSet(0,-1,0,0,0,0,0),
        VSMALL
    );

    // calculate reaction rate from class
    surfaceReacRate_.update();
    const volScalarField& reacRate = surfaceReacRate_.surfaceReacRate();
    scalar surfaceReacFlag = 0.0;
    if(surfaceReaction_)
    {
        surfaceReacFlag = 1.0;
    }

    volScalarField massFluxPrec
    (
        "massFluxPrec",
        reacRate
        * Mv_
        * pos(from-0.001)
        * surfaceReacFlag
    );

    
    if (mesh_.time().outputTime())
    {
        areaDensity_.write();
        Cmask_.write();
        //to.write();
    }

    dimensionedScalar totReactionRate_(0.0);

    // calculate total reaction rate from class
    // !!! Need to be implemented
    totReactionRate_ = gSum((mesh_.V()*massFluxPrec*areaDensity_)());

    Info<< "precipitate return: " << min(massFluxPrec*areaDensity_).value() << ", "<< max(massFluxPrec*areaDensity_).value() << endl;
    Info<< "Total reaction rate: " << totReactionRate_.value() << endl;

    if(smoothSurface_)
    {
        Info<< "Smooth surface areaDensity method used!" << endl;
        return massFluxPrec * areaDensitySmooth;
    }
    else
    {
        Info<< "Gradient surface areaDensity method used!" << endl;
        return massFluxPrec * areaDensity_ * Cmask_;
    }
}

Foam::tmp<Foam::fvScalarMatrix> 
Foam::phaseChangeReaction::preRect
(
    const volScalarField& C,
    volScalarField& SuOutPut
)
{
    Info<< "Start preRect calculation" << endl;
    
    tmp<fvScalarMatrix> cEqnPtr
    (
        new fvScalarMatrix(C, dimMoles/dimTime)
    );

    fvScalarMatrix& eqn = cEqnPtr.ref();
    // Net mass transfer from k to i phase
    tmp<volScalarField> cdmdtNetki
    (
        new volScalarField
        (
            IOobject
            (
                "cdmdtYki",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );
    volScalarField& dmdtNetki = cdmdtNetki.ref();

    // Explicit mass transfer rate
    Info<< "Calculate KExp" << endl;
    tmp<volScalarField> KExp = Kexp(C);
    Info<< "End of calculation of KExp" << endl;

    // Nucleation induced volumetric mass transfer
    // !!! Need to be implemented

    if (KExp.valid())
    {
        Info<< "Fluid to solid calculated!" << endl;
        dmdtNetki += KExp.ref();
        Info<< "Kexp output min/max: " << min(KExp.ref()) << ", " << max(KExp.ref()) << endl;
    }

    SuOutPut = dmdtNetki;
    eqn -= (dmdtNetki/Mv_);

    return cEqnPtr;
}

void Foam::phaseChangeReaction::addInterfacePorosity(fvVectorMatrix& UEqn)
{
    //- Current function will be moved to a new class in the future version

    const scalarField& Vc = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    //- Voller Prakash interfacial porosity model
    const volScalarField& liquidAlpha = alpha_;

    tmp<volScalarField> STerm 
    (
        new volScalarField
        (
            IOobject
            (
                "STerm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& STermRef = STerm.ref();

    STermRef = Cu_*sqr(1.0-liquidAlpha)/(pow3(liquidAlpha) + 1e-3)*Cmask_;

    Udiag += Vc*STermRef;
}

void Foam::phaseChangeReaction::updateCmask()
{
    //- Readjust in next version for better coding
    const volScalarField& from = alpha_;

    volScalarField to
    (
        "to",
        1.0-alpha_
    );
    to.correctBoundaryConditions();

    const volVectorField gradFrom(fvc::grad(from));
    const volVectorField gradTo(fvc::grad(to));
    const volScalarField gradAlphaf(gradFrom & gradTo);

    Info<< "calculate Cmask..." <<endl;

    // mark cells that is next to the solid phase cell
    // within same processor
    forAll(Cmask_, celli)
    {
        if (gradAlphaf[celli] < 0)
        {
            if (from[celli] > alphaMin_ && from[celli] < alphaMax_)
            {
                if(mag(gradFrom[celli])>gradLim_)
                {
                    Cmask_[celli] = 1.0;
                }
            }
        }
        //- check nearby cell within same processor
        bool flag = false;
        labelListList cellCellsHolder;
        cellCellsHolder.clear();

        forAll(mesh_.cellCells()[celli],cellj)
        {
            if (to[celli] > alphaSolidMin_)
            {
                Cmask_[mesh_.cellCells()[celli][cellj]] = 1.0;

                //- store cellCells lists of adjacent cells
                label dirNearbyCellMarker = mesh_.cellCells()[celli][cellj];
                cellCellsHolder.resize(cellCellsHolder.size()+1);
                labelList dirNearbyCellCells = mesh_.cellCells()[dirNearbyCellMarker];
                cellCellsHolder.append(dirNearbyCellCells);
                //Info<< dirNearbyCellCells << endl;

                flag = true;
            }
        }

        //- enable corner cell growth for cartisian grid
        if(flag == true && cornerCell_ == true)
        {
            labelList cornerCellList;
            cornerCellList.clear();
            forAll(cellCellsHolder, iterI)
            {
                forAll(cellCellsHolder, iterJ)
                {
                    if(iterI != iterJ)
                    {
                        forAll(cellCellsHolder[iterI], labelI)
                        {
                            if(cellCellsHolder[iterJ].found(cellCellsHolder[iterI][labelI]))
                            {
                                if(!cornerCellList.found(cellCellsHolder[iterI][labelI]))
                                {
                                    cornerCellList.append(cellCellsHolder[iterI][labelI]);
                                }
                            }                            
                        }
                    }
                }
            }
            //Info<< "corner cell list: " << cornerCellList << endl;


            forAll(cornerCellList, cellI)
            {
                Cmask_[cornerCellList[cellI]] = 1.0;
            }
        }

        //- nucleation site enable growth
        // if(nucleationFlag_)
        // {
        //     if(nuSite[celli]==1.0)
        //     {
        //         Cmask_[celli] = 1.0;
        //     }
        // }
        //if (flag)
        //{
        //    Cmask_[celli] = 1.0;
        //}
    }

    //- check mesh boundaries of the processor
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const fvPatchScalarField& pf = to.boundaryField()[patchI];
        const labelList& faceCells = pf.patch().faceCells();

        //- Coupled boundaries (processor, cylic, etc)
        if(pf.coupled())
        {
            scalarField neighbors = pf.patchNeighbourField();

            forAll(faceCells, faceI)
            {
                //if (neighbors[faceI]>alphaSolidMin_)
                //{
                //    Info<< faceI << endl;
                //    Info<< neighbors[faceI] << endl;
                //}
                
                if (neighbors[faceI] >= alphaSolidMin_)
                {
                    Cmask_[faceCells[faceI]] = 1.0;
                }
            }
        }
    } 

    Cmask_.correctBoundaryConditions();
}

Foam::phaseChangeReaction::~phaseChangeReaction(){};

// ************************************************************************* //
