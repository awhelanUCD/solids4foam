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
    linearGTN

Description
    Small strain  Gurson-Tveergard-Needleman (GTN) micro-mechanical damage model. This class also provides users the option to include the non-local porosity f.

    Small strain elasto-plastic model with GTN yield equation uses algorithmic approach similar to that outlined in Achouri et Al. (Doi: 10.1016/j.commatsci.2013.01.035).

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    The plastic hardening law is defined in a user file (see plasticStrainVsYieldStress in tutorial case) 

    For GTN parameters, the user can specify:
        - the yield equation paramters q1 q2 and q3
        - void nucleation paramters epsilonN, sN, fN
        - void coalescence parameters f0, fC, fU and fF;
        - shear parameter kW
    
    For non-local porosity the user can specify the characteristic length (lc) 

SourceFiles
    linearGTN.C

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef linearGTN_H
#define linearGTN_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"
//#include "fadiff.h"
//using namespace fadbad;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class linearGTN Declaration
\*---------------------------------------------------------------------------*/

class linearGTN
:
    public mechanicalLaw
{
    // Private data

        const scalar relaxationFactor_;

        //Switch to include coalescene of voids effects
        Switch includeCoalescence_;

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

        const Switch debug_;

 
        // Parameters for GTN yield equation
        dimensionedScalar q1_;
        dimensionedScalar q2_;
        dimensionedScalar q3_;

        // Parameters for the nucleation of voids
        dimensionedScalar fN_;
        dimensionedScalar epsilonN_;
        dimensionedScalar sN_;

        // Initial porosity
        dimensionedScalar f0_;

        // Parameters for void coalescence
        dimensionedScalar fC_;        
        dimensionedScalar fU_;
        dimensionedScalar fF_;

        //- Table of post-yield stress versus plastic strain
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Switch to enable solution of pressure Laplacian equation
        //  This can help quell oscillations in the hydrostatic stress
        const Switch solvePressureEqn_;

        //- Optional: it is possible to scale the amount of smoothing in the
        //  pressure equation with this coefficient
        scalar pressureSmoothingCoeff_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigma)/3.0 == -p,
        //  where p is the hydrostatic pressure
        volScalarField sigmaHyd_;

        //- Yield stress
        volScalarField sigmaY_;

        volScalarField f_;

        volScalarField fStar_;

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

        volScalarField hydroEpsilonP_;
        volScalarField DHydroEpsilonP_;
        volScalarField qEpsilonP_;
        volScalarField DqEpsilonP_;

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
        linearGTN(const linearGTN&);

        //- Disallow default bitwise assignment
        void operator=(const linearGTN&);

        // update porosity
        void calculatef
        (
            scalar& DHyd,    
            scalar& f,   
            scalar& fOld,   
            scalar& fStar,    
            scalar& fStarOld,
            scalar epsilonPEqOld, 
            scalar& DEpsilonPEq    
        ) const;

        //- Iteratively calculate deviatoric, spherical and equivalent plastic strain
        void newtonLoop
        (            
            scalar& DHyd,   
            scalar& DqEpsilonP, 
            scalar& fStar,       
            scalar& fStarOld,
            scalar& sigmaY,
            scalar epsilonPEqOld, 
            scalar& DEpsilonPEq,    
            scalar& pTrial,
            scalar& qTrial
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
    TypeName("linearGTN");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearGTN
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~linearGTN();


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
