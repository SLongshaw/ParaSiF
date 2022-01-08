/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::parabolicInletVelocityTwoDFvPatchVectorField

Description
    This boundary condition describes an inlet vector boundary condition in
    cylindrical co-ordinates given a central axis, central point, rpm, axial
    and radial velocity.

Usage
    \table
        Property     | Description             | Required    | Default value
        axis         | axis of rotation        | yes         |
        origin       | origin of rotation      | yes         |
        axialVelocity | axial velocity profile [m/s] | yes    |
        radialVelocity | radial velocity profile [m/s] | yes  |
        rpm          | rotational speed (revolutions per minute) | yes|
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            parabolicInletVelocity;
        axis            (0 0 1);
        origin          (0 0 0);
        axialVelocity   constant 30;
        radialVelocity  constant -10;
        rpm             constant 100;
    }
    \endverbatim

Note
    The \c axialVelocity, \c radialVelocity and \c rpm entries are Function1
    types, able to describe time varying functions.  The example above gives
    the usage for supplying constant values.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    parabolicInletVelocityTwoDFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

// No license
// ----------

// BAE-FSI
// parabolicInletVelocityTwoDFvPatchVectorField.H
//
// This is a part of the Partitioned Multi-physical Simulation Framework (parMupSiF)
//
// FEniCSV2019.1.0+ <-> MUI(mui4py) <-> MUI(C++) <-> OpenFOAMV6+ two way Coupling Module.
//
// Incompressible Navier-Stokes equations for fluid domain in OpenFOAM
// Structure dynamics equations for structure domain in FEniCS
//
// parabolicInletVelocityTwoDFvPatchVectorField.H is a C++ header file to create a 2D parabolic inlet velocity with a constant velocity magnitude
//  for the CFD solver of the OpenFOAM pert of the parMupSiF code
//
// located in the src/CFD/OpenFOAM/V6/parmupsifLibs/parmupsifBC/parabolicInletVelocityTwoD sub-folder of the parMupSiF folder
//
// Last changed: 25-September-2019
//
// Author(s): W.L
// Contact: wendi.liu@stfc.ac.uk
//
// Copyright 2019 UK Research and Innovation
// IBM Confidential
// OCO Source Materials
// 5747-SM3
// (c) Copyright IBM Corp. 2017, 2019
// The source code for this program is not published or otherwise 
// divested of its trade secrets, irrespective of what has 
// been deposited with the U.S. Copyright Office. 
//
// All rights reserved

#ifndef parabolicInletVelocityTwoDFvPatchVectorField_H
#define parabolicInletVelocityTwoDFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
         Class parabolicInletVelocityTwoDFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class parabolicInletVelocityTwoDFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Axial velocity
        const scalar maxVelocity_;

        //- Period to Maximum Velocity
 //       const scalar Tperiod_;

        //- Central point
        const vector centre_;

        //- Axis
        const vector axis_;
		
        //- Origin
        const vector origin_;

public:

   //- Runtime type information
   TypeName("parabolicInletVelocityTwoD");


   // Constructors

        //- Construct from patch and internal field
        parabolicInletVelocityTwoDFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicInletVelocityTwoDFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        parabolicInletVelocityTwoDFvPatchVectorField
        (
            const parabolicInletVelocityTwoDFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        parabolicInletVelocityTwoDFvPatchVectorField
        (
            const parabolicInletVelocityTwoDFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicInletVelocityTwoDFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicInletVelocityTwoDFvPatchVectorField
        (
            const parabolicInletVelocityTwoDFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicInletVelocityTwoDFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
