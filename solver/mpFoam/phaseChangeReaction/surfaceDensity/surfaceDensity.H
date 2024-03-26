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

Class
    Foam::InterfaceCompositionModel

Description
    Base model for calculating surface area density (m^2/m^3)

SourceFiles
    surfaceDensity.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceDensity_H
#define surfaceDensity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "localMin.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
    typedef List<vector> vectorList;
/*---------------------------------------------------------------------------*\
                  Class surfaceDensity Declaration
\*---------------------------------------------------------------------------*/

class surfaceDensity
{
protected:

    // Protected data

        //- Reaction model switcher
        const word areaDensityModel_;

        //- Alpha field (dimless)
        const volScalarField& alpha_;

        //- Area density field (m^2/m^3)
        volScalarField areaDensity_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Dict reference
        const dictionary& dict_;

    // Protected member functions

        //- Inputs reading function
        //- Calculate volume-based area density from liquid phase
        virtual void calAreaDensity();

public:

    //- Runtime type information
    //TypeName("surfaceDensity");

    //- Constructors
    //- Construct from components
    surfaceDensity
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const volScalarField& alpha,
        const volScalarField& solid
    );

    //- Destructor
    virtual ~surfaceDensity();

    //- Member Functions

        //- Print inputs info
        virtual void info();
    
        //- Update the area density field
        virtual void update();

        //- Return the area density field
        virtual volScalarField& areaDensity()
        {
            return areaDensity_;
        }

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
