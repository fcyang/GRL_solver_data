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
    solidNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef solidNucleation_H
#define solidNucleation_H

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
                  Class solidNucleation Declaration
\*---------------------------------------------------------------------------*/

class solidNucleation
{
protected:

    // Protected data

        //- Nucleation model switcher
        const word nucleationModel_;

        //- Concentration field (dimless)
        const volScalarField& C_;

        //- Dict reference
        const dictionary& dict_;

        //- liquid phase Vof
        const volScalarField& liquid_;

        //- Solid phase Vof
        const volScalarField& solid_;

        //- Species name
        const word speciesName_;

        //- Nucleation site marker
        volScalarField nuSite_;

        //- Local nucleation rate (n/m^2/s)
        volScalarField nuRateOut_;

        //- Total nucleation counter
        dimensionedScalar nuTotal_;

        //- lnA prefactor for nucleation rate
        const scalar lnA_;

        //- B geometric factor for nucleation rate
        const scalar B_;

        //- gamma interfacial energy for nucleation
        const scalar gamma_;

        //- molar volume for nucleation
        const scalar Vmol_;

        //- Nucleation box max/min point
        const bool nucleationBoxFlag_;
        vector nucMinPoint_;
        vector nucMaxPoint_;

        //- Debug flag
        const bool debug_;

    // Protected member functions

        //- Inputs reading function
        //- Calculate local surface nucleation rate
        virtual tmp<volScalarField> calNucRate();

        //- Convert molar concentration to SI
        virtual tmp<volScalarField> CtoSI();

        //- Randomize crystal nucleation process
        virtual void nuSiteCal();

        //- Check whether a cell is within defined box
        virtual bool cellInsideTheBox(point& pos) const;

        //- Reading parameters
        virtual void readParameters();

public:

    //- Runtime type information
    //TypeName("solidNucleation");

    //- Constructors
    //- Construct from components
    solidNucleation
    (
        const dictionary& dict,
        const volScalarField& C,
        const volScalarField& liquid,
        const volScalarField& solid
    );

    //- Destructor
    virtual ~solidNucleation();

    //- Member Functions

        //- Print inputs info
        virtual void info();
    
        //- Update the nucleation rate field basing on C field
        virtual void update();

        //- Return the nucleation rate field
        virtual volScalarField& nuRate()
        {
            return nuRateOut_;
        }

        //- Return the nucleation site field per step
        virtual volScalarField& nuSitePerStep()
        {
            return nuSite_;
        }

        //- Return the total nucleation number
        virtual dimensionedScalar& nuTotal()
        {
            return nuTotal_;
        }

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
