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
    Base class for interface composition models, templated on the two
    thermodynamic models either side of the interface.

SourceFiles
    phaseChangeReaction.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeReaction_H
#define phaseChangeReaction_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceDensity.H"
#include "surfaceReaction.H"
#include "solidNucleation.H"
#include "fvMatricesFwd.H"
#include "localMin.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
    typedef List<vector> vectorList;
/*---------------------------------------------------------------------------*\
                  Class phaseChangeReaction Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeReaction
{
protected:

    // Protected data

        //- Input dictionary
        const dictionary& dict_;

        //- Concentration field
        const volScalarField& C_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Solid phase Vof
        const volScalarField& solid_;

        //- Phase name
        const word phaseName_;
        
        //- Directional growth option
        bool dirGrowthFlag_;

        //- Directional growth number
        const scalar dirGrowthNum_;

        //- Directional growth norms
        const List<vector> dirGrowthNorm_;

        //- Directional growth similarity prefactor
        const List<scalar> dirFactor_;

        //- Constant interfacial porosity coefficient
        const scalar Cu_;

        //- pKsp for dissolution/precipitation species
        const scalar pKsp_;

        //- lnA prefactor for nucleation rate
        const scalar lnA_;

        //- B geometric factor for nucleation rate
        const scalar B_;

        //- gamma interfacial energy for nucleation
        const scalar gamma_;

        //- molar volume for nucleation
        const scalar Vmol_;

        //- Reaction model switcher
        const word reacModel_;

        //- Precipitation coefficient
        dimensionedScalar K_;

        //- Directional precipitation coefficient
        const List<scalar> KList_;

        //- Molar weight of the vapour in the continous phase
        dimensionedScalar Mv_;

        //- 'To' phase maximum value for the mass transfer
        scalar alphaMax_;

        //- 'To' phase minumum value for the mass transfer
        scalar alphaMin_;

        //- 'To' phase minumum value for the surface identification
        scalar alphaSolidMin_;

        //- Alpha maximum for the rest of phases
        scalar alphaRestMax_;

        //- Surface areaDensity option
        bool smoothSurface_;

        //- User input smooth surface areaDensity constant [m2/m3]
        scalar smoothAreaDensity_;

        //- Surface area density class
        surfaceDensity surfaceAreaDensity_;

        //- Surface reaction rate class
        surfaceReaction surfaceReacRate_;

        //- Solid Nucleation class
        solidNucleation solidNucleation_;

        //- Option for enable corner cell Cmask growth
        bool cornerCell_;

        //- User input vof gradient limit for enable growth [1/m]
        scalar gradLim_;

        //- Alpha field
        const volScalarField& alpha_;

        //- Cmask for solid phase region
        volScalarField& Cmask_;

        //- debug option
        bool debug_;

        //- Surface reaction option
        const bool surfaceReaction_;

        //- Nucleation site option
        const bool nucleation_;

        //- Nucleation box max/min point
        const bool nucleationBoxFlag_;
        vector nucMinPoint_;
        vector nucMaxPoint_;

    // Protected member functions

        //- Print inputs info
        virtual void info();

public:

    //- Runtime type information
    //TypeName("phaseChangeReaction");

    //- Constructors
        //- Construct from components
    phaseChangeReaction
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const volScalarField& C,
        const volScalarField& alpha,
        const volScalarField& solid,
        volScalarField& Cmask
    );

    //- Destructor
    virtual ~phaseChangeReaction();

    //- Member Functions

    //- Calculate volume-based total mass conversion rate
    virtual tmp<volScalarField> Kexp
    (
        const volScalarField& field
    );

    //- Add interfacial penalty term for UEqn
    void addInterfacePorosity(fvVectorMatrix& UEqn);

    //- Update heterogeneous reaction region for solid phase
    void updateCmask();

    //- Return mass reaction term for CEqn
    virtual tmp<fvScalarMatrix> preRect(const volScalarField& C, volScalarField& SuOutput);

    //- Return nucleation site basing on concentration field
    //virtual void nuSiteCal(volScalarField& nuSitei, volScalarField& nuRateOut, vectorList& nuSiteList, dimensionedScalar& nuTotal_);
    //virtual tmp<volScalarField> CtoSI();
    //virtual tmp<volScalarField> nuRateCal();

    //- Geometry constrain for crystal growth
    //virtual vectorList extractNuSite(const volScalarField& nuSite, vectorList& nuSiteList);

    //- Nucleation class return functions
    //- Return the nucleation rate field
    virtual volScalarField& nuRate()
    {
        return solidNucleation_.nuRate();
    }

    //- Return the nucleation site field per step
    virtual volScalarField& nuSitePerStep()
    {
        return solidNucleation_.nuSitePerStep();
    }

    //- Return the total nucleation number
    virtual dimensionedScalar& nuTotal()
    {
        return solidNucleation_.nuTotal();
    }

    //- Update the nucleation rate field
    void updateNucRate()
    {
        solidNucleation_.update();
    }
    
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
