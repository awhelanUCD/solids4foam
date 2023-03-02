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
    GreenStrainLemaitreDamage

Description

    Lemaitre damage model coupled with large-strain elasto-plastic model.
    This class also give users the option to incorporate non-local gradient damage.

     - Neo-Hookean elasticity with Mises/J2 plasticity  (Simo & Hughes, Computational Inelasticity, 1998, Springer.).

     - Lemaitre damage model 

     - non local gradient damage (Doi: Peerlings et al.  10.1016/S0020-7683(01)00087-7 )

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    The plastic hardening law is defined in a user file (see plasticStrainVsYieldStress in tutorial case) 

    For Lemaitre damage parameters, the user can specify:
        - Damage denominator (s0), damage exponent (b), plastic strain threshold (epsilonD) and critical damage parameter (damageC)

    For non-local gradient damage model, the user can specify:
        - the characteristic length (l) (for local damage model lc=0)


SourceFiles

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef GreenStrainLemaitreDamage_H
#define GreenStrainLemaitreDamage_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class GreenStrainLemaitreDamage Declaration
\*---------------------------------------------------------------------------*/

class GreenStrainLemaitreDamage
:
    public mechanicalLaw
{
    // Private data

        const scalar relaxationFactor_;

        //- Reference yield stress
/*        const scalar y0_;*/
/*        // Isotropic hardening slope*/
/*        const scalar hIso_;*/

/*        // Hardening rate*/
/*        const scalar beta_;*/
/*        //- Asymptotic yield stress*/
/*        const scalar yInf_;*/


        //- Characteristic length associated with non-local damage
        dimensionedScalar charLength_;

        //- Plastic strain threshold for damage growth
        dimensionedScalar epsilonD_;

        //- Initial density
        const dimensionedScalar rho_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Lemaitre exponent
        dimensionedScalar b_;

        //- Lemaitre denominator
        dimensionedScalar s0_;

        //- Critical damage parameter
        const scalar damageC_;

        const Switch debug_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Relative deformation gradient volField
        volTensorField* relFPtr_;

        //- Relative deformation gradient surfaceField
        surfaceTensorField* relFfPtr_;

        //- Jacobian of the toal deformation gradient volField
        volScalarField* JPtr_;

        //- Jacobian of the toal deformation gradient surfaceField
        surfaceScalarField* JfPtr_;

        //- Total deformation gradient volField
        volTensorField* FPtr_;

        //- Total deformation gradient surfaceField
        surfaceTensorField* FfPtr_;

        //- Table of post-yield stress versus plastic strain
        //  For large strain total Lagrangian models, it expects
        //  the true stress and true plastic strain.
        interpolationTable<scalar> stressPlasticStrainSeries_;

        volScalarField plasticMult_;

        //- Damage
        volScalarField damage_;
 
        //- Non-local damage
        volScalarField damageNonLocal_;

        //- Gradient of the non-local damage
        volVectorField gradDamageNonLocal_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigma)/3.0 == -p,
        //  where p is the hydrostatic pressure
        volScalarField sigmaHyd_;

        //- Switch to enable solution of a pressure equation, which ensures a
        //  smooth sigmaHyd field
        Switch smoothPressure_;

        //- Cauchy yield stress
        volScalarField sigmaY_;

        //- Incremental change in sigmaY
        volScalarField DSigmaY_;

        //- Total plastic strain
        volSymmTensorField epsilonP_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        volSymmTensorField bEbarTrial_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        volSymmTensorField bEbar_;

        //- Equivalent plastic strain increment
        volScalarField DEpsilonPEq_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField DLambda_;

        //- Equivalent plastic strain
        volScalarField epsilonPEq_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;

        //- An iterative procedure is used when the yield stress is a nonlinear
        //  function of plastic strain
       // const bool nonLinearPlasticity_;

        //- Update bEbar consistently with the assumption that det(bEbar) == 1
        //  defaults to off
        const Switch updateBEbarConsistent_;

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
        GreenStrainLemaitreDamage(const GreenStrainLemaitreDamage&);

        //- Disallow default bitwise assignment
        void operator=(const GreenStrainLemaitreDamage&);

        //- Make the relF field
        void makeRelF();

        //- Return a reference to the relF field
        volTensorField& relF();

        //- Make the J field
        void makeJ();

        //- Return a reference to the J field
        volScalarField& J();

        //- Make the F field
        void makeF();

        //- Return a reference to the F field
        volTensorField& F();

        //- Calculate the damage value
        void damageFunction
        (
           const scalar triaxiality,    
           const scalar epsilonPEq,        
           const scalar sigmaEq,      
           scalar& damage,
           const scalar damageOld,
           const scalar& damageNonLocal,
           const scalar DLambda
        ) const;

        //- Return the current Kirchhoff yield stress
        scalar curYieldStress
        (
            const scalar plasticMultOld,    // Current equivalent plastic strain
            const scalar DLambda, 
            const scalar J                // Current Jacobian
  
        ) const;

        //- Evaulate current value of the yield function
        scalar yieldFunction
        (
            const scalar plasticMultOld,    // Old equivalent plastic strain
            const scalar magSTrial,        // Deviatoric trial stress magnitude
            const scalar DLambda,          // Plastic multiplier
            const scalar muBar,            // Scaled shear modulus
            const scalar J,                 // Current Jacobian
            const scalar damage
        ) const;

        //- Iteratively calculate plastic multiplier increment (DLambda)
        //  and current yield stress using Newton's method
        void newtonLoop
        (
            scalar& DLambda,               // Plastic multiplier
            scalar& curSigmaY,             // Current yield stress
            const scalar plasticMultOld,    // Old equivalent plastic strain
            const scalar magSTrial,        // Deviatoric trial stress magnitude
            const scalar muBar,            // Scaled shear modulus
            const scalar J,                // Current Jacobian
            const scalar maxMagDEpsilon,    // Max strain increment magnitude
            const scalar damage
        ) const;

        //- Calcualte Ibar such that det(bEbar) == 1
        tmp<volScalarField> Ibar
        (
            const volSymmTensorField& devBEbar
        );

        //- Calcualte Ibar such that det(bEbar) == 1
        tmp<surfaceScalarField> Ibar
        (
            const surfaceSymmTensorField& devBEbar
        );

public:

    //- Runtime type information
    TypeName("GreenStrainLemaitreDamage");

    // Static data members


    // Constructors

        //- Construct from dictionary
        GreenStrainLemaitreDamage
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~GreenStrainLemaitreDamage();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

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
