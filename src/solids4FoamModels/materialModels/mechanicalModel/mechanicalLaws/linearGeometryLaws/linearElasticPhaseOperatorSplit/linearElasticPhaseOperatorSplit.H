/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    linearElasticPhaseOperatorSplit

Description
    Hookean linear elastic mechanical law with brittle phase field fracture.
    
    - Based on approach developed by Miehe et Al. (Doi:  10.1016/j.cma.2010.04.011)

     For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
        or
        - Shear modulus (mu) and bulk modulus (K)

     For phase-field fracture model, the user can specify:
        - Critical fracture energy (Gc), characteristic length (lc) and a residual stiffness


SourceFiles
    linearElasticPhaseOperatorSplit.C

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD,
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef linearElasticPhaseOperatorSplit_H
#define linearElasticPhaseOperatorSplit_H

#include "mechanicalLaw.H"
#include "Eigen/Dense"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class linearElasticPhaseOperatorSplit Declaration
\*---------------------------------------------------------------------------*/

class linearElasticPhaseOperatorSplit
:
    public mechanicalLaw
{
    // Private data

        // Material parameters for phase field fracture evolution
        dimensionedScalar charLength_;
        dimensionedScalar viscousD_;
        dimensionedScalar Gc_;
        dimensionedScalar residualStiff_;

        //- Density
        dimensionedScalar rho_;

        //- Shear modulus (second Lame parameter)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- First Lame parameter
        dimensionedScalar lambda_;

        Switch monolithic_;

        //- Switch to enable solution of pressure Laplacian equation
        //  This can help quell oscillations in the hydrostatic stress
        Switch solvePressureEqn_;

        //- Optional: it is possible to scale the amount of smoothing in the
        //  pressure equation with this coefficient
        scalar pressureSmoothingCoeff_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigma)/3.0 == -p,
        //  where p is the hydrostatic pressure
        volScalarField sigmaHyd_;

        //Phase field variable d
        volScalarField Dpf_;

        // Elastic strain energy
        volScalarField strainEnergy_;

        // History strain energy
        volScalarField HStrainEnergy_;

        //- Total strain
        volSymmTensorField epsilon_;

        //- Positive strain tensor
        volSymmTensorField epsilonPositive_;

        // Negative strain tensor
        volSymmTensorField epsilonNegative_;

        //- Initial (residual) stress field
        //  Defaults to zero if not specified
        const volSymmTensorField sigma0_;

        //- Initial (residual) stress surface field
        mutable surfaceSymmTensorField* sigma0fPtr_;

    // Private Member Functions

        //- Make the initial stress surface field
        void makeSigma0f() const;

        //- Return a reference to the initial stress surface field
        const surfaceSymmTensorField& sigma0f() const;

       // Function for decomposing the strain tensor into it's positive and negative components
       void decompose
       (
         const symmTensor& epsilon,
         symmTensor& epsilonPositive,
         symmTensor& epsilonNegative
       );

       // Function for calculating the strain energy
       scalar calcStrainEnergy
       (
           const symmTensor& epsilon
       );

        //- Calculate hydrostatic component of the stress tensor
        void calculateHydrostaticStress
        (
            volScalarField& sigmaHyd,
            const volScalarField& trEpsilon
        );

        //- Calculate hydrostatic component of the stress tensor
        void calculateHydrostaticStress
        (
            surfaceScalarField& sigmaHyd,
            const surfaceScalarField& trEpsilon
        );

        //- Disallow default bitwise copy construct
        linearElasticPhaseOperatorSplit(const linearElasticPhaseOperatorSplit&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticPhaseOperatorSplit&);

public:

    //- Runtime type information
    TypeName("linearElasticPhaseOperatorSplit");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElasticPhaseOperatorSplit
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~linearElasticPhaseOperatorSplit();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Provide access to mu for the coupledUnsNonLinlinearElasticPhaseOperatorSplitSolid
        virtual const dimensionedScalar& mu() const;

        //- Bulk modulus
        virtual const dimensionedScalar& K() const;

        //- Young's modulus
        virtual const dimensionedScalar& E() const;

        //- Poisson's ratio
        virtual const dimensionedScalar& nu() const;

        //- Lambda
        virtual const dimensionedScalar& lambda() const;

        virtual void calcPhase();

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
