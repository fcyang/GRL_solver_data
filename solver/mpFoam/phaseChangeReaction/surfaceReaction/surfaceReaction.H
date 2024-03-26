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
    Base model for calculating surface reaction rate (mol/m^2/s)

SourceFiles
    surfaceReaction.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceReaction_H
#define surfaceReaction_H

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
                  Class surfaceReaction Declaration
\*---------------------------------------------------------------------------*/

class surfaceReaction
{
protected:

    // Protected data

        //- Reaction model switcher
        const word reacModel_;

        //- Concentration field (dimless)
        const volScalarField& C_;

        //- Dict reference
        const dictionary& dict_;

        //- Solid phase Vof
        const volScalarField& solid_;

        //- Species name
        const word& speciesName_;

        //- Surface reaction rate calculated basing on concentration field (mol/m^2/s)
        volScalarField surfaceReacRate_;

        //- Species index in input file
        const scalar index_;

    // Protected member functions

        //- Inputs reading function
        //- Calculate volume-based area density from liquid phase
        virtual void calReacRate();

        //- Convert concentration to SI
        virtual tmp<volScalarField> CtoSI();

public:

    //- Runtime type information
    //TypeName("surfaceReaction");

    //- Constructors
    //- Construct from components
    surfaceReaction
    (
        const dictionary& dict,
        const volScalarField& C,
        const volScalarField& solid,
        const word& speciesName
    );

    //- Destructor
    virtual ~surfaceReaction();

    //- Member Functions

        //- Print inputs info
        virtual void info();
    
        //- Update the area density field
        virtual void update();

        //- Return the area density field
        virtual volScalarField& surfaceReacRate()
        {
            return surfaceReacRate_;
        }

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
