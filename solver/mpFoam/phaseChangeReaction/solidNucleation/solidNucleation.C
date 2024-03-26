/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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

#include "solidNucleation.H"
#include "fvcGrad.H"
#include <random>
#include "constants.H"
#include <fstream>


using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidNucleation::
solidNucleation
(
    const dictionary& dict,
    const volScalarField& C,
    const volScalarField& liquid,
    const volScalarField& solid
)
:
    nucleationModel_(dict.lookupOrDefault<word>("nucleationModel", "wall")),
    C_(C),
    dict_(dict),
    liquid_(liquid),
    solid_(solid),
    speciesName_(dict.lookupOrDefault<word>("phaseName", "test")),
    nuSite_
    (
        IOobject
        (
            "nuSite_"+speciesName_,
            C_.time().timeName(),
            C_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C_.mesh(),
        dimensionedScalar("nuSite"+speciesName_, dimless, Zero)
    ),
    nuRateOut_
    (
        IOobject
        (
            "nuRateOut_"+speciesName_,
            C_.time().timeName(),
            C_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C_.mesh(),
        dimensionedScalar("nuRateOut"+speciesName_, dimless/dimArea/dimTime, Zero)
    ),
    nuTotal_(0.0),
    lnA_(dict.lookupOrDefault<scalar>("lnA", 1.0)),
    B_(dict.lookupOrDefault<scalar>("B", 1.0)),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.0)),
    Vmol_(dict.lookupOrDefault<scalar>("Vmol", 1.0)),
    nucleationBoxFlag_(dict.lookupOrDefault<bool>("nucleationBox", false)),
    nucMinPoint_(dict.lookupOrDefault<vector>("nucMinPoint", vector(0,0,0))),
    nucMaxPoint_(dict.lookupOrDefault<vector>("nucMaxPoint", vector(0,0,0))),
    debug_(dict.lookupOrDefault<bool>("debug", false))
{
    readParameters();
    info();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidNucleation::~solidNucleation()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::solidNucleation::readParameters()
{
    // Read parameters from dictionary dict_
}

void Foam::solidNucleation::info()
{
    Info<< "B parameters: " << B_ << endl;
    Info<< "lnA parameters: " << lnA_ << endl;
    Info<< "gamma parameters: " << gamma_ << endl;
    Info<< "Vmol parameters: " << Vmol_ << endl;
}

Foam::tmp<Foam::volScalarField> 
Foam::solidNucleation::calNucRate()
{
    tmp<volScalarField> nuRateTmp 
    (
        new volScalarField
        (
            IOobject
            (
                "nuRateTmp",
                C_.time().timeName(),
                C_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            C_.mesh(),
            dimensionedScalar(dimless/dimArea/dimTime, Zero)
        )
    );
    volScalarField& nuRateRef = nuRateTmp.ref();

    // obtain SI field from concentration
    tmp<volScalarField> SI = CtoSI();    
    Info<< "SI output min/max: " << min(SI.ref()) << ", " << max(SI.ref()) << endl;

    // define constants
    const scalar kb = 1.38064852e-23;
    const scalar NA = 6.02214076e23;

    forAll(C_,cellI)
    {
        // check if cell within the box defined for nucleation
        vector pos = C_.mesh().C()[cellI];

        // calculate 
        if(cellInsideTheBox(pos))
        {
            nuRateRef[cellI] = exp(lnA_)*exp((-B_*pow(gamma_,3.0)*pow(Vmol_/NA,2.0))/(2.303*2.303*pow(SI.ref()[cellI],2.0)*pow(kb*293,3.0)));
        }
        else
        {
            nuRateRef[cellI] = 0.0;
        }
        // Info<< SI.ref()[cellI] << endl;
        // Info<< nuRateRef[cellI] << endl;
    }

    return nuRateTmp;
}

void Foam::solidNucleation::nuSiteCal()
{
    // Create tmp nucleation field
    tmp<volScalarField> nuSitePerStep
    (
        new volScalarField
        (
            IOobject
            (
                "nuSitePerStep",
                C_.time().timeName(),
                C_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            C_.mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& nuSitePerStepRef = nuSitePerStep.ref();

    // Create tmp face area field
    tmp<volScalarField> faceAreaTmp
    (
        new volScalarField
        (
            IOobject
            (
                "faceAreaTmp",
                C_.time().timeName(),
                C_.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            C_.mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& faceAreaTmpRef = faceAreaTmp.ref();

    // Calculate the possibilities for nucleation
    tmp<volScalarField> nuRateTmp = calNucRate();
    nuRateOut_ = nuRateTmp.ref();
    Info<< "Max/Min nucleation rate(n/m2/s): " << max(nuRateOut_).value() << ", " << min(nuRateOut_).value()<< endl;

    // Loop through wall boundary patches to estimate nucleation
    const fvPatchList& patches = C_.mesh().boundary();
    labelList wallList;
    wallList.clear();

    scalarList faceAreaList;
    faceAreaList.clear();

    // Obtain face area from boundary patch
    forAll(patches, patchI)
    {
        const fvPatch& p = patches[patchI];
        if(p.type()==nucleationModel_)
        {
            forAll(p, pFaceI)
            {
                label faceCelli = p.faceCells()[pFaceI];
                wallList.resize(wallList.size()+1);
                wallList[wallList.size()-1]=faceCelli;
                scalar faceCellAreai = C_.mesh().magSf().boundaryField()[patchI][pFaceI];
                faceAreaList.resize(faceAreaList.size()+1);
                faceAreaList[faceAreaList.size()-1]=faceCellAreai;
            }
        }
    }
    //Pout << "Wall list size: "<< wallList.size() << endl;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::default_random_engine generator(rd());
    
    dimensionedScalar totSufArea(0.0);
    dimensionedScalar nuRateAccum(0.0);
    dimensionedScalar averNuRate(0.0);

    forAll(wallList, i)
    {
        if (nucleationBoxFlag_)
        {
            vector pos = C_.mesh().C()[wallList[i]];
            // Nucleation restriction using a box
            if(cellInsideTheBox(pos))
            {
                const cell& faces = C_.mesh().cells()[wallList[i]];
                faceAreaTmpRef[wallList[i]] = faceAreaList[i];

                // Generate random number using uniform_int_distribution
                std::random_device rd;  //Will be used to obtain a seed for the random number engine
                std::default_random_engine generator(rd());
                scalar upperLim = 1.0/(faceAreaTmpRef[wallList[i]]+VSMALL);
                std::uniform_int_distribution<long long unsigned> dis(1, static_cast<long long unsigned>(upperLim));
                scalar randNum = dis(generator);

                // Generate nucleation site basing on nucleation rate and random number
                if((nuRateOut_[wallList[i]])>randNum)
                {
                    nuSitePerStepRef[wallList[i]] = 1.0;
                    if(liquid_[wallList[i]] > 0.5)
                    {
                        //nuTotal_.value() += 1.0;
                        Info<< "New nucleation site found: " << nuRateOut_.ref()[wallList[i]] << " > " << randNum << endl;
                    }
                }
            }
        }
        else
        {
            const cell& faces = C_.mesh().cells()[wallList[i]];
            faceAreaTmpRef[wallList[i]] = faceAreaList[i];

            // Generate random number using uniform_int_distribution
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::default_random_engine generator(rd());
            scalar upperLim = 1.0/(faceAreaTmpRef[wallList[i]]+VSMALL);
            std::uniform_int_distribution<long long unsigned> dis(1, static_cast<long long unsigned>(upperLim));
            scalar randNum = dis(generator);

            // Generate nucleation site basing on nucleation rate and random number
            if((nuRateOut_[wallList[i]])>randNum)
            {
                nuSitePerStepRef[wallList[i]] = 1.0;
                if(liquid_[wallList[i]] > 0.5)
                {
                    //nuTotal_.value() += 1.0;
                    Info<< "New nucleation site found: " << nuRateOut_[wallList[i]] << " > " << randNum << endl;
                }
            }
        }


        if(debug_)
        {
            //Info<< "Max faceArea: " << faceArea << endl;
            //cout<< "Total cell counts per unit area: " << upper << endl;
            //Info<< "Random integer generated: " << randNum << endl;
            //Info<< "Nucleation rate: " << nuRate.ref()[wallList[i]] << endl;
            Info<< "Computational time: " << C_.mesh().time().value() << endl;
        }
    }

    // pass markers from tmp field to nuSite
    forAll(nuSite_,cellI)
    {
        nuSite_[cellI] = nuSitePerStepRef[cellI];
    }
    Info<< "Nucleation site within time step: " << max(nuSite_).value() << ", " << min(nuSite_).value()<< endl;
    
    Info<< "Max/Min cell surface Area(m2): " << max(faceAreaTmpRef).value() << ", " << min(faceAreaTmpRef).value()<< endl;

    nuTotal_.value() += gSum((nuSitePerStepRef)());

    totSufArea.value() = gSum((faceAreaTmpRef*liquid_*pos(liquid_-0.1))())+VSMALL;
    nuRateAccum.value() = gSum((faceAreaTmpRef*liquid_*nuRateOut_*pos(liquid_-0.1))());

    Info<< "Total nucleation sites: " << nuTotal_.value() << endl;
    averNuRate.value() = nuRateAccum.value()/totSufArea.value();

    // output available surface area
    if(Pstream::master() && C_.mesh().time().outputTime())
    {
        string fileName = "nucleationSurfaceArea.csv";
        std::ifstream file(fileName);
        std::fstream fout;
        fout.open(fileName, std::ios::out | std::ios::app);
        if (file.peek() == std::ifstream::traits_type::eof())
        {
            fout << "Time[s],totSufArea,averNuRate" << "\n";
        }
        fout << C_.mesh().time().timeName() << ", "
            << totSufArea.value() << ", "
            << averNuRate.value() << "\n";
    }

    if(debug_)
    {
        if (C_.mesh().time().outputTime())
        {
            faceAreaTmpRef.write();
            //wallMarkerRef.write();
        }
    }
}

Foam::tmp<Foam::volScalarField> 
Foam::solidNucleation::CtoSI
()
{
    // Read the constant reaction value k defined in dict
    const scalar pKsp_(readScalar(dict_.lookup("pKsp")));

    tmp<volScalarField> SItmp 
    (
        new volScalarField
        (
            IOobject
            (
                "SItmp",
                C_.time().timeName(),
                C_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            C_.mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& SItmpRef = SItmp.ref();

    scalar Ksp = pow(10,-pKsp_);
    Info<< "KSP value: " << Ksp << endl;

    // mol/m3 to SI
    forAll(C_, cellI)
    {
        SItmpRef[cellI] = max(0.001, log10((pow(C_[cellI]/1000,2.0)+VSMALL)/Ksp));
    };

    return SItmp;
}


bool Foam::solidNucleation::cellInsideTheBox(point& pos) const
{
    bool res = true;
    
    for(int i=0; i<3; i++)
    {
        res = res  &&  ( pos[i]>=nucMinPoint_[i] && pos[i]<=nucMaxPoint_[i] );
    }
    
    return res;
}

void Foam::solidNucleation::update()
{
    nuSiteCal();
}

// ************************************************************************* //
