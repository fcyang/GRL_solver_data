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

#include "surfaceDensity.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceDensity::
surfaceDensity
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarField& alpha,
    const volScalarField& solid
)
:
    areaDensityModel_(dict.lookupOrDefault<word>("areaDensityModel", "constant")),
    alpha_(alpha),
    areaDensity_
    (
        IOobject
        (
            "areaDensity",
            alpha.time().timeName(),
            alpha.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha.mesh(),
        dimensionedScalar("areaDensityInput", dimless/dimLength, Zero)
    ),
    mesh_(mesh),
    dict_(dict)
{
    info();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceDensity::~surfaceDensity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::surfaceDensity::info()
{
    Info<< "Surface area density model for DBS >> " << areaDensityModel_ << endl;
}

void Foam::surfaceDensity::calAreaDensity()
{
    volScalarField liquidPhase 
    (
        "liquidPhase",
        alpha_
    );
    liquidPhase.correctBoundaryConditions();

    if(areaDensityModel_ == "constant")
    {
        // Read the constant smooth area density value (m^2/m^3) defined in dict
        const scalar areaDensityFactor_(readScalar(dict_.lookup("areaDensityFactor")));
        
        // Assign uniform value to area density field
        forAll(areaDensity_, cellI)
        {
            areaDensity_[cellI] = areaDensityFactor_;
        }

        Info<< "areaDensity min/max : "<<min(areaDensity_).value()<<", "<< max(areaDensity_).value()<<endl;
    }
    else if(areaDensityModel_ == "gradient")
    {
        const volVectorField gradLiquid(fvc::grad(liquidPhase));

        forAll(areaDensity_, cellI)
        {
            areaDensity_[cellI] = 2.0*mag(gradLiquid[cellI]);
        }

        Info<< "areaDensity min/max : "<<min(areaDensity_).value()<<", "<< max(areaDensity_).value()<<endl;
    }
    else if(areaDensityModel_ == "exponent")
    {
        // Read the constant smooth area density value (m^2/m^3) defined in dict
        const scalar areaDensityFactor_(readScalar(dict_.lookup("areaDensityFactor")));
        const scalar areaDensityExp_(readScalar(dict_.lookup("areaDensityExp")));

        forAll(areaDensity_, cellI)
        {            
            areaDensity_[cellI] = areaDensityFactor_*pow((1.0-liquidPhase[cellI]), areaDensityExp_);
        }

        Info<< "areaDensity min/max : "<<min(areaDensity_).value()<<", "<< max(areaDensity_).value()<<endl;
    }
    else if(areaDensityModel_ == "hemiSphere")
    {
        // Single cell hemisphere shape crystal surface area calculation
        // Read prefactors
        const scalar areaDensityFactor_(readScalar(dict_.lookup("areaDensityFactor")));
        const scalar tmpPi = 3.1415926;

        forAll(areaDensity_, cellI)
        {            
            areaDensity_[cellI] = areaDensityFactor_*2.0*tmpPi*pow((3.0/2.0/tmpPi*(1.0-liquidPhase[cellI])*mesh_.V()[cellI]), 2.0/3.0)/mesh_.V()[cellI];
        }

        Info<< "areaDensity min/max : "<<min(areaDensity_).value()<<", "<< max(areaDensity_).value()<<endl;
    }
    else if(areaDensityModel_ == "mixed")
    {
        // Not implemented yet
    }
    else
    {
        SeriousErrorIn("surfaceDensity::calAreaDensity")
              << "There is no defined areaDensityModel in dictionary; "
              << "Current usable model options:\n"
              << "constant\n"
              << "gradient\n"
              << exit(FatalError);
    }
}

void Foam::surfaceDensity::update()
{
    calAreaDensity();
}

// ************************************************************************* //
