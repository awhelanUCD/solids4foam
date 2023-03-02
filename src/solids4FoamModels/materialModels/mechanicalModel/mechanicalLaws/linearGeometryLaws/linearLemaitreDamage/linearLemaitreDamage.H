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
    linearLemaitreDamage

Description
    Elasto-plastic mechanical law as described by Simo & Hughes (1998)
    in Box 3.2 coupled with Lemaitre damage model. The user can choose to incorporate non-local damage aswell.

    Linear (Hookean) elasticity with Mises/J2 plasticity.

    Lemaitre danage model (Lemaitre ISBN 978-3-540- 21503-5)

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    The plastic hardening law is defined in a user file (see plasticStrainVsYieldStress in tutorial case)

    For damage parameters the user can specify 
    
        - The damage denominator (s0), the damage exponent (b) and the plastic strian threshold (epsilonD)

    For non-local gradient damage model, the user can specify:
        - the characteristic length (l) (for local damage model lc=0)


SourceFiles
    linearLemaitreDamage.C

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef linearLemaitreDamage_H
#define linearLemaitreDamage_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class linearLemaitreDamage Declaration
\*---------------------------------------------------------------------------*/

class linearLemaitreDamage
:
    public mechanicalLaw
{
    // Private data

        //characteristic length for Lemaitre damage law
        dimensionedScalar charLength_;

        //- Initial density
        const dimensionedScalar rho_;

        const scalar relaxationFactor_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Table of post-yield stress versus plastic strain
        interpolationTable<scalar> stressPlasticStrainSeries_;

        // Lemaitre damage exponent
        dimensionedScalar b_;
        
        // Lemaitre damage denominator
        dimensionedScalar s0_;

        // plastic strian threshold for damage evolution
        dimensionedScalar epsilonD_;

        //critical damage parameter
        scalar damageC_;

        //- Table of post-yield stress versus plastic strain
/*        interpolationTable<scalar> stressPlasticStrainSeries_;*/

        //- Switch to enable solution of pressure Laplacian equation
        //  This can help quell oscillations in the hydrostatic stress
        const Switch solvePressureEqn_;

        const Switch debug_;


        //- Optional: it is possible to scale the amount of smoothing in the
        //  pressure equation with this coefficient
        scalar pressureSmoothingCoeff_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigma)/3.0 == -p,
        //  where p is the hydrostatic pressure
        volScalarField sigmaHyd_;

        //- Yield stress
        volScalarField sigmaY_;

        //damage
        volScalarField damage_;
 
        // non-local damage
        volScalarField damageNonLocal_;

        //gradient of non-local damage
        volVectorField gradDamageNonLocal_;

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

        // accumulated plastic multiplier
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
        linearLemaitreDamage(const linearLemaitreDamage&);

        //- Disallow default bitwise assignment
        void operator=(const linearLemaitreDamage&);

        //calcualte damage
        void damageFunction
        (
            scalar triaxiality,    
           const scalar epsilonPEq,        
           const scalar sigmaEq,       
           const scalar damageOld,
           scalar& damage,
           const scalar DepsilonPeq,
           const scalar DLambda,
           const scalar sigmaHyd
        ) const;

        scalar curYieldStress
        (
            const scalar plasticMultOld,  
            const scalar DLambda
        ) const;

        //- Evaulate current value of the yield function
        scalar yieldFunction
        (
            const scalar qTrial,      
            const scalar DLambda, 
            const scalar plasticMultOld,
            scalar damage, 
            const scalar muBar           
        ) const;

        void newtonLoop
        (
            scalar& DLambda,             
            scalar& curSigmaY, 
            const scalar plasticMult,
            const scalar qTrial,       
            const scalar muBar,           
            const scalar maxMagDEpsilon,   
            scalar damage
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
    TypeName("linearLemaitreDamage");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearLemaitreDamage
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~linearLemaitreDamage();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return the implicit stiffness as a diagTensor
        virtual tmp<volDiagTensorField> impKdiagTensor() const;

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
