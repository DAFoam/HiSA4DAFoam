/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2011-2012 OpenFOAM Foundation

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::characteristicInletOutletPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a characteristic inlet/outlet condition
    suitable for all Mach numbers.
    Inlet: Prescribed velocity and temperature (and pressure for M > 1)
    Outlet: Prescribed pressure (for M < 1)
    See characteristicBase for usage.

SourceFiles
    characteristicInletOutletPressureFvPatchScalarField.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef characteristicInletOutletPressureFvPatchScalarFields_H
#define characteristicInletOutletPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "characteristicBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class characteristicInletOutletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class characteristicInletOutletPressureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public characteristicBase
{
public:

    //- Runtime type information
    TypeName("characteristicInletOutletPressure");


    // Constructors

        //- Construct from patch and internal field
        characteristicInletOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        characteristicInletOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given characteristicInletOutletPressureFvPatchScalarField
        //  onto a new patch
        characteristicInletOutletPressureFvPatchScalarField
        (
            const characteristicInletOutletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        characteristicInletOutletPressureFvPatchScalarField
        (
            const characteristicInletOutletPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new characteristicInletOutletPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        characteristicInletOutletPressureFvPatchScalarField
        (
            const characteristicInletOutletPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new characteristicInletOutletPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
