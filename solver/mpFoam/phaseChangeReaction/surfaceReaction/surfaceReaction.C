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

#include "surfaceReaction.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceReaction::
surfaceReaction
(
    const dictionary& dict,
    const volScalarField& C,
    const volScalarField& solid,
    const word& speciesName
)
:
    reacModel_(dict.lookupOrDefault<word>("reacModel", "constant")),
    C_(C),
    dict_(dict),
    solid_(solid),
    speciesName_(speciesName),
    surfaceReacRate_
    (
        IOobject
        (
            "surfaceReacRate"+speciesName,
            C_.time().timeName(),
            C_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
            ),
        C_.mesh(),
        dimensionedScalar("surfaceReacRateInput"+speciesName, dimMoles/dimArea/dimTime, Zero)
    ),
    index_(dict.lookupOrDefault<scalar>("test",1.0))
{
    info();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceReaction::~surfaceReaction()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::surfaceReaction::info()
{
    Info<< "Surface reaction model for DBS >> " << reacModel_ << endl;
}

void Foam::surfaceReaction::calReacRate()
{
    if(reacModel_ == "constant")
    {
        // Read the constant reaction value k defined in dict
        const scalar K_(readScalar(dict_.lookup("K")));
        
        // Read the corresponding saturation concentration
        const scalar Csat_(readScalar(dict_.lookup("Csat")));

        // Calculate reaction rate basing on given parameters and concentration field
        forAll(surfaceReacRate_, cellI)
        {
            surfaceReacRate_[cellI] = K_;
        }

        Info<< "Surface reaction rate min/max : "
            << min(surfaceReacRate_).value() <<", "<< max(surfaceReacRate_).value()<<endl;
    }
    else if(reacModel_ == "linear")
    {
        // Read the constant reaction value k defined in dict
        const scalar K_(readScalar(dict_.lookup("K")));
        
        // Read the corresponding saturation concentration
        const scalar Csat_(readScalar(dict_.lookup("Csat")));

        // Calculate reaction rate basing on given parameters and concentration field
        forAll(surfaceReacRate_, cellI)
        {
            // Check if dissolution or precipitation
            if(C_[cellI]<Csat_)
            {
                // Dissolution only happens within cell with solid phase
                surfaceReacRate_[cellI] = K_*(C_[cellI]-Csat_)*pos(solid_[cellI]-0.01);
            }
            else
            {
                // Precipitation only happens within cell with over saturated concentration
                surfaceReacRate_[cellI] = K_*(C_[cellI]-Csat_)*pos(C_[cellI]-Csat_);
            }
        }

        Info<< "Surface reaction rate min/max : "
            << min(surfaceReacRate_).value() <<", "<< max(surfaceReacRate_).value()<<endl;
    }
    else if(reacModel_ == "second")
    {
        // Read the constant reaction value k defined in dict
        const scalar K_(readScalar(dict_.lookup("K")));
        
        // Read the corresponding saturation concentration
        const scalar Csat_(readScalar(dict_.lookup("Csat")));

        // Calculate reaction rate basing on given parameters and concentration field
        forAll(surfaceReacRate_, cellI)
        {
            surfaceReacRate_[cellI] = K_*sqr(C_[cellI]-Csat_)*sign(C_[cellI]-Csat_);
        }

        Info<< "Surface reaction rate min/max : "
            << min(surfaceReacRate_).value() <<", "<< max(surfaceReacRate_).value()<<endl;
    }
    else if(reacModel_ == "SI")
    {
        // Not yet implemented
    }
    else
    {
        SeriousErrorIn("surfaceReaction::calReacRate")
              << "There is no defined reacModel in dictionary; "
              << "Current usable model options:\n"
              << "constant\n"
              << "linear\n"
              << "second\n"
              << exit(FatalError);
    }
}

Foam::tmp<Foam::volScalarField> 
Foam::surfaceReaction::CtoSI
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

void Foam::surfaceReaction::update()
{
    calReacRate();
}

// ************************************************************************* //
