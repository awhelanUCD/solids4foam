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
    linearElasticPlasticPhaseFracturePhaseFracture

Description
    Ductile phase-field fracture model coupled with small strain elasto-plastic model. The phase field is described by the variable d with initial value (unbroken value) 1

     - Small strain elasticity with Mises/J2 plasticity  (Simo & Hughes (1998) in Box 3.2).

     - Ductile phase-field fracure model (Borden et Al. Doi: 10.1016/j.cma.2016.09.005)

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    The plastic hardening law is defined in a user file (see plasticStrainVsYieldStress in tutorial case) 

    For phase-field fracture model, the user can specify:
        - Critical fracture energy (Gc), characteristic length (lc), plastic work threshold threshold (w0)

SourceFiles
    linearElasticPlasticPhaseFracture.C

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticPlasticPhaseFracture_H
#define linearElasticPlasticPhaseFracture_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class linearElasticPlasticPhaseFracture Declaration
\*---------------------------------------------------------------------------*/

class linearElasticPlasticPhaseFracture
:
    public mechanicalLaw
{
    // Private data

        //- Initial density
        const dimensionedScalar rho_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        // parameters for phase-field fracture model
        dimensionedScalar charLength_;
        //dimensionedScalar viscousD_;
        dimensionedScalar Gc_;
        dimensionedScalar w0_;
        scalar residualStiff_;

        //- Table of post-yield stress versus plastic strain
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Switch to enable solution of pressure Laplacian equation
        //  This can help quell oscillations in the hydrostatic stress
        const Switch solvePressureEqn_;

        const Switch debug_;

        Switch monolithic_;

        //- Optional: it is possible to scale the amount of smoothing in the
        //  pressure equation with this coefficient
        scalar pressureSmoothingCoeff_;

        // phase field variable d
        volScalarField Dpf_;

        // elastic strain energy
        volScalarField strainEnergy_;

        // plastic strain energy
        volScalarField plasticStrainEnergy_;

        // history variable for the elastic strain energy
        volScalarField HStrainEnergy_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigma)/3.0 == -p,
        //  where p is the hydrostatic pressure
        volScalarField sigmaHyd_;

        // magnitude of the deviatoric stress
        volScalarField magS_;

        //- Yield stress
        volScalarField sigmaY_;

        //- Incremental change in sigmaY
        volScalarField DSigmaY_;

        //- Total strain
        volSymmTensorField epsilon_;

        //- Total plastic strain
        volSymmTensorField epsilonP_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        //- Equivalent plastic strain increment
        volScalarField DEpsilonPEq_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField DLambda_;

        volScalarField one_;

        volScalarField zero_;

        volScalarField plasticMult_;

        //- Equivalent plastic strain
        volScalarField epsilonPEq_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Tolerance for Newton loop
        static scalar LoopTol_;

        //- Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        //- finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticPlasticPhaseFracture(const linearElasticPlasticPhaseFracture&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticPlasticPhaseFracture&);

        //- Iteratively calculate plastic multiplier increment (DLambda)
        //  and current yield stress using Newton's method
        void newtonLoop
        (
          scalar& DLambda,             
          scalar& curSigmaY,
          const scalar& plasticMult,  
          const scalar& magSTrial,        
          const scalar& Dpf 
        ) const;

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


public:

    //- Runtime type information
    TypeName("linearElasticPlasticPhaseFracture");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElasticPlasticPhaseFracture
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~linearElasticPlasticPhaseFracture();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return the implicit stiffness as a diagTensor
        virtual tmp<volDiagTensorField> impKdiagTensor() const;

        // update the phase field variable 
        virtual void calcPhase();

        //- Update the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Update the stress surface field
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();

        //- Return the desired new time-step
        virtual scalar newDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
