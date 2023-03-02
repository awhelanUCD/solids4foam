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

\*---------------------------------------------------------------------------*/


#include <iostream>
#include<stdio.h>
#include "linearGTN.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
//#include "fadiff.h"
#include <Eigen/Dense>
using namespace Eigen;
//#include "fadone.H"
using namespace std;

//using namespace fadbad;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearGTN, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearGTN, linGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar linearGTN::LoopTol_ = 1e-4;

    // Maximum number of iterations for Newton loop
    label linearGTN::MaxNewtonIter_ = 1000;

    // finiteDiff is the delta for finite difference differentiation
    scalar linearGTN::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar linearGTN::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::linearGTN::calculatef
(
    scalar& DHyd, // Increment of spherical strain
    scalar& f, // porosity
    scalar& fOld, // Old porosity
    scalar& fStar, // effective porosity        
    scalar& fStarOld, // old value of the effective porosity
    scalar epsilonPEqOld, // old equivalent plastic strain 
    scalar& DEpsilonPEq  // equivalent plastic strain increment  
) const
{
    scalar A = (fN_.value()/(sN_.value()*sqrt(2*3.14)))
               *pow(2.71828,-0.5*pow((epsilonPEqOld+DEpsilonPEq-epsilonN_.value())/sN_.value(), 2.0));

    f=fOld+(1-fStar)*DHyd+A*DEpsilonPEq;

    fStar=f;
    
    if (f >= fC_.value() && includeCoalescence_)
    {
       fStar = fC_.value()+(f-fC_.value())*(1/q1_.value()-fC_.value())/(fF_.value()-fC_.value());
    }

}
void Foam::linearGTN::newtonLoop
(
    scalar& DHyd, // spherical plastic strain increment
    scalar& DqEpsilonP, // Deviatoric plastic strain increment
    scalar& fStar, // effective porosity
    scalar& fStarOld, //old effective porosity
    scalar& sigmaY, // yield stress
    scalar epsilonPEqOld, //old equivalent plastic strain value
    scalar& DEpsilonPEq, // increment of the equivalent plastic strain    
    scalar& pTrial, // trial pressure
    scalar& qTrial // trial equivalent stress
) const
{
    // set initial values for the deviatoric, spherical and equivalent plastic strain
    DqEpsilonP=0;             
    DHyd=0;               
    DEpsilonPEq=0; 

        //Newton loop to solve for speherical, deviatoric and equivalent plastic strain increments
        for (int j=0; j<=1000; j++) 
        {

            //update pressure and equivalent stress
            scalar p = pTrial - K_.value()*DHyd;
            scalar q = qTrial - 3*mu_.value()*DqEpsilonP;

            //update yield stress
            sigmaY = stressPlasticStrainSeries_(epsilonPEqOld+DEpsilonPEq);

            // initialse fcosh and fsinh - these are components of the GTN yield equation
            scalar fCosh = cosh(3*q2_.value()*p/(2*sigmaY));
            scalar fSinh = sinh(3*q2_.value()*p/(2*sigmaY));

            //derivaitive of the yield stress
            scalar dSigmaY = (stressPlasticStrainSeries_(epsilonPEqOld+DEpsilonPEq+1e-8) 
                            -stressPlasticStrainSeries_(epsilonPEqOld+DEpsilonPEq))/1e-8;


            //these derivatives are components of the derivatives which fill the Jacobian matrix                
            scalar df1dq = 2.0*(q/pow(sigmaY,2.0));
            scalar df1dp = 3*fStar*q1_.value()*q2_.value()*fSinh/sigmaY;
            scalar df1dSigmaY =- 2.0*pow(q,2.0)/pow(sigmaY,3.0) - df1dp*p/sigmaY;
            scalar df1dpdp = 9.0*q1_.value()*pow(q2_.value(),2.0)*fStar*fCosh/(2.0*pow(sigmaY,2.0));
            scalar df1dpdSigmaY = -9.0*q1_.value()*pow(q2_.value(),2.0)*
                                   fStar*p*fCosh/(2.0*pow(sigmaY,3.0))
                                  -3*fStar*q1_.value()*q2_.value()*fSinh/pow(sigmaY,2.0);
            scalar df1dqdSigmaY = -4.0*(q/pow(sigmaY,3.0));

            //Calcualte derivatives which constitute the Jacobian matrix
            scalar a1 = -3.0*mu_.value()*df1dq;
            scalar a2 = -K_.value()*df1dp;
            scalar a3 = df1dSigmaY*dSigmaY;

            scalar b1 = df1dp+3.0*DqEpsilonP*mu_.value()*2.0/pow(sigmaY,2.0);
            scalar b2 = -DqEpsilonP*df1dpdp*K_.value() - df1dq;
            scalar b3 = DqEpsilonP*df1dpdSigmaY*dSigmaY - DHyd*df1dqdSigmaY*dSigmaY;

            scalar c1 = -q+3.0*mu_.value()*DqEpsilonP;
            scalar c2 = -p+K_.value()*DHyd;
            scalar c3 = (1-fStar)*(sigmaY+dSigmaY*DEpsilonPEq);

            if (debug_)
            {
                Info << "a1= " << a1 << endl;
	        Info << "a2= " << a2 << endl;
	        Info << "a3= " << a3 << endl;

	        Info << "b1= " << b1 << endl;
	        Info << "b2= " << b2 << endl;
	        Info << "b3= " << b3 << endl;
 
	        Info << "c1= " << c1 << endl;
	        Info << "c2= " << c2 << endl;
	        Info << "c3= " << c3 << endl;
            }
            
            //calculate values of fucnctions f1, f2 and f3 
            scalar f1val = pow(q/sigmaY, 2.0) + 2*q1_.value()*fStar*fCosh - (1+q3_.value()*pow(fStar,2.0));
            scalar f2val = DqEpsilonP*df1dp-DHyd*df1dq;
            scalar f3val = (1-fStar)*sigmaY*DEpsilonPEq - q*DqEpsilonP-p*DHyd;

            //create 3x3 Jacobian matrix
            Matrix3d J;
            J << a1, a2, a3, 
                 b1, b2, b3,
                 c1, c2, c3;

            if (debug_){cout<< "The inverse of J is:" << J.inverse() << endl;Info<<""<<endl;}
            if (debug_){cout<< "J test:" << J.inverse()*J << endl; Info << "" << endl;
                        Info << "" << endl;}

            //create 3 dimensional vector
            Vector3d fx(f1val,f2val,f3val);

            //find inverse of matrix
            Matrix3d Jinv = J.inverse();

            //Gauss-Newton step
            Vector3d hgn = -Jinv*fx;

            // update values of deviatoric, spherical and equivalent plastic strain increments
            DqEpsilonP = DqEpsilonP + (hgn(0));
            DHyd = DHyd + (hgn(1));
            DEpsilonPEq = DEpsilonPEq + (hgn(2));
   
            if (debug_)
            {

                Info << "f1val= " << f1val << endl;
	        Info << "f2val= " << f2val << endl;
	        Info << "f3val= " << f3val << endl;

	        Info << "h0: " <<hgn(0) << endl;
	        Info  << "h1: " <<hgn(1) << endl;
	        Info  << "h2: " << hgn(2)<< endl;

     	        Info << "DqEpsilonP: " << DqEpsilonP << endl;
	        Info << "DHyd: " << DHyd << endl;
	        Info << "DEpsilonPEq: " << DEpsilonPEq<< endl;}

            if (j == 1000)
            {
                WarningIn("logStrainGTN::newtonLoop()")
                    << "Plasticity Newton loop not converging" << endl;
            }

            if ((mag(f1val) < LoopTol_) && (mag(f2val) < LoopTol_) && (mag(f3val) < LoopTol_)) {break;}
        
    }

}

void Foam::linearGTN::calculateHydrostaticStress
(
    volScalarField& sigmaHyd,
    const volScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        // Store previous iteration to allow relaxation, if needed
        sigmaHyd.storePrevIter();

        // Lookup the momentum equation inverse diagonal field
        const volScalarField* ADPtr = NULL;
        if (mesh().foundObject<volScalarField>("DEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DEqnA");
        }
        else if (mesh().foundObject<volScalarField>("DDEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DDEqnA");
        }
        else
        {
            FatalErrorIn
            (
                "void Foam::linearGTN::"
                "calculateHydrostaticStress\n"
                "(\n"
                "    volScalarField& sigmaHyd,\n"
                "    const volScalarField& trEpsilon\n"
                ")"
            )   << "Cannot find the DEqnA or DDEqnA field: this should be "
                << "stored in the solidModel" << abort(FatalError);
        }
        const volScalarField& AD = *ADPtr;

        // Pressure diffusivity field multiple by (4.0/3.0)*mu + K, which is
        // equivalent to (2*mu + lambda)
        // Note: we can scale this coefficient by pressureSmoothingCoeff to
        // provide greater smoothing
        const surfaceScalarField rDAf
        (
            "rDAf",
            pressureSmoothingCoeff_
           *fvc::interpolate
            (
                ((4.0/3.0)*mu_ + K_)/AD, "interpolate(grad(sigmaHyd))"
            )
        );

        // Solve pressure laplacian
        // Note: the the laplacian and div terms combine to produce a
        // third-order smoothing/dissipation term
        // If we only used the laplacian term then the smoothing/dissipation
        // would be second-order.
        // It would be interesting to see how this compares to the JST 2nd/4th
        // order dissipation term
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(1.0, sigmaHyd)
          - fvm::laplacian(rDAf, sigmaHyd, "laplacian(DA,sigmaHyd)")
          + fvc::div(rDAf*fvc::interpolate(fvc::grad(sigmaHyd)) & mesh().Sf())
         ==
            K_*trEpsilon
        );

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Relax the field
        sigmaHyd.relax();
    }
    else
    {
        // Directly calculate hydrostatic stress from displacement field
        sigmaHyd = K_*trEpsilon;
    }
}


void Foam::linearGTN::calculateHydrostaticStress
(
    surfaceScalarField& sigmaHyd,
    const surfaceScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        FatalErrorIn
        (
            "void Foam::linearGTN::calculateHydrostaticStress\n"
            "(\n"
            "    surfaceScalarField& sigmaHyd,\n"
            "    const surfaceScalarField& trEpsilon\n"
            ")"
        )   << "'solvePressureEqn' option only implemented for volField stress "
            << "calculation" << abort(FatalError);
    }
    else
    {
        // Directly calculate hydrostatic stress from displacement field
        sigmaHyd = K_*trEpsilon;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearGTN::linearGTN
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    includeCoalescence_
    (
        dict.lookupOrDefault<Switch>("includeCoalescence", false)
    ),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    debug_(dict.lookup("debug")),
    q1_("zero", dimless, 0.0),
    q2_("zero", dimless, 0.0),
    q3_("zero", dimless, 0.0),
    fN_("zero", dimless, 0.0),
    epsilonN_("zero", dimless, 0.0),
    sN_("zero", dimless, 0.0),
    f0_("zero", dimless, 0.0),
    fC_("zero", dimless, 0.0),
    fU_("zero", dimless, 0.0),
    fF_("zero", dimless, 0.0),
    stressPlasticStrainSeries_(dict),
    solvePressureEqn_(dict.lookup("solvePressureEqn")),
    pressureSmoothingCoeff_
    (
        dict.lookupOrDefault<scalar>("pressureSmoothingCoeff", 1.0)
    ),
    sigmaHyd_
    (
        IOobject
        (
            "sigmaHyd",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    ),
    f_
    (
        IOobject
        (
            "f",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    fStar_
    (
        IOobject
        (
            "fStar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    DSigmaY_
    (
        IOobject
        (
            "DSigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonPEq_
    (
        IOobject
        (
            "DEpsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    hydroEpsilonP_
    (
        IOobject
        (
            "HydroEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DHydroEpsilonP_
    (
        IOobject
        (
            "DHydroEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    qEpsilonP_
    (
        IOobject
        (
            "qEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DqEpsilonP_
    (
        IOobject
        (
            " DqEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    plasticN_
    (
        IOobject
        (
            "plasticN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    epsilon_.oldTime();
    epsilonP_.oldTime();
    epsilonPEq_.oldTime();
    plasticN_.oldTime();
    fStar_.oldTime();

    q1_ = dimensionedScalar(dict.lookup("q1"));
    q2_ = dimensionedScalar(dict.lookup("q2"));
    q3_ = dimensionedScalar(dict.lookup("q3"));
    fN_ = dimensionedScalar(dict.lookup("fN")); 
    epsilonN_ = dimensionedScalar(dict.lookup("epsilonN"));
    sN_ = dimensionedScalar(dict.lookup("sN"));
    f0_ = dimensionedScalar(dict.lookup("f0"));

    f_=f0_;
    fStar_=f0_;

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        E_ = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        nu_ = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E_/(2.0*(1.0 + nu_));

        // Set the bulk modulus
        K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        // Read shear modulus
        mu_ = dimensionedScalar(dict.lookup("mu"));

        // Read bulk modulus
        K_ = dimensionedScalar(dict.lookup("K"));

        // Calculate Young's modulus
        E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

        // Calculate Poisson's ratio
        nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));
    }
    else
    {
        FatalErrorIn
        (
            "linearGTN::linearGTN::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingCoeff: " << pressureSmoothingCoeff_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearGTN::~linearGTN()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearGTN::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "lawRho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::linearGTN::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial =
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL));

    // Calculate scaling factor
    const volScalarField scaleFactor = 1.0 - (2.0*mu_*DLambda_/magSTrial);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            scaleFactor*(4.0/3.0)*mu_ + K_
        )
    );
}


Foam::tmp<Foam::volDiagTensorField>
Foam::linearGTN::impKdiagTensor() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial =
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL));

    // Calculate scaling factor
    const volScalarField theta = 1.0 - (2.0*mu_*DLambda_/magSTrial);

    // Calculate N squared where N is the plastic return direction
    const volTensorField NsquaredTensor = plasticN_ & plasticN_;
    volDiagTensorField Nsquared
    (
        IOobject
        (
            "Nsquared",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedDiagTensor("zero", dimless, diagTensor::zero)
    );

    Nsquared.internalField() = diag(NsquaredTensor.internalField());

    forAll(Nsquared.boundaryField(), patchI)
    {
        Nsquared.boundaryField()[patchI] =
            diag(NsquaredTensor.boundaryField()[patchI]);
    }

    const diagTensor Idiag = diagTensor(1.0, 1.0, 1.0);

    return tmp<volDiagTensorField>
    (
        new volDiagTensorField
        (
            IOobject
            (
                "impKdiagTensor",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            K_*Idiag + mu_*0.5*theta*(Idiag*4.0/3.0 - 2.0*Nsquared)
            //K_*Idiag + mu_*theta*(Idiag*4.0/3.0)
            //K_*Idiag + mu_*(theta/theta)*(Idiag*4.0/3.0)
        )
    );
}


void Foam::linearGTN::correct(volSymmTensorField& sigma)
{
    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilon_ = epsilon_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilon_ = symm(gradD);
    }

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearGTN::"
                "correct(volSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        epsilon_.replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
          - (
                epsilonP_.component(symmTensor::XX)
              + epsilonP_.component(symmTensor::YY)
            )
        );
    }

    // trial elastic strain
    volSymmTensorField epsilonElTrial = epsilon_-epsilonP_.oldTime();
 
    // initialise elastic strain
    volSymmTensorField epsilonEl = epsilonElTrial;

    //Calculate trial Hydrostatic Strain
     volScalarField hydroEpsilonTrial = tr(epsilonElTrial);

    //Calculate trial Hydrostatic Stress
    volScalarField  sigmaHydTr = K_*tr(epsilonElTrial);

    // Calculate deviatoric strain
    volSymmTensorField e = dev(epsilonElTrial);

    // Calculate deviatoric trial stress
    volSymmTensorField sTrial = 2.0*mu_*e;
    
    // trial equivalent stress
    volScalarField qTrial = sqrt(3.0/2.0)*mag(sTrial);

    // Calculate the yield function
    const volScalarField fTrial = pow(qTrial/sigmaY_, 2.0) 
                                 + 2*q1_.value()*fStar_*cosh((3*sigmaHydTr*q2_.value())/(2*sigmaY_))
                                 - (1+pow(q2_.value()*fStar_,2.0));

    volScalarField q = qTrial;
    volSymmTensorField s = sTrial;

    // Store previous iteration for residual calculation
    DEpsilonP_.storePrevIter();

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    const scalarField& sigmaYI = sigmaY_.internalField();
    const scalarField& sigmaYOldI = sigmaY_.oldTime().internalField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().internalField();

    scalarField& pI=sigmaHyd_.internalField();
    scalarField& qTrialI = qTrial.internalField();
    symmTensorField& sigmaI = sigma.internalField();
    symmTensorField& sI = s.internalField();
    scalarField& qI = q.internalField();
    scalarField& fOldI = f_.oldTime().internalField();
    scalarField& fI = f_.internalField();
    scalarField& fStarI = fStar_.internalField();
    scalarField& fStarOldI = fStar_.oldTime().internalField();
    scalarField& pTrialI = sigmaHydTr.internalField();
    scalarField& DhydroEpsilonPI = DHydroEpsilonP_.internalField();
    scalarField& DqEpsilonPI = DqEpsilonP_.internalField();
    scalarField& DEpsilonPEqI = DEpsilonPEq_.internalField();

    // Calculate deviatoric, spherical and equivalent plastic strain
    forAll(fTrialI, cellI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[cellI]);
        if (magS > SMALL)
        {
            plasticNI[cellI] = (3*sTrialI[cellI])/(2*qTrialI[cellI]);
        }
        
        // Update variables when no plastic deformation
        if (fTrialI[cellI] < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
            plasticNI[cellI] = symmTensor::zero;

            pI = pTrialI[cellI];
            fI[cellI] = fOldI[cellI];
            fStarI[cellI] = fStarOldI[cellI];
            sI[cellI] = sTrialI[cellI];
            qI[cellI] = qTrialI[cellI];
            sigmaI[cellI] = pI[cellI]*I+sTrialI[cellI];
            DEpsilonPEqI[cellI] = 0;

        }
        else
        {
            scalar curSigmaY = sigmaYOldI[cellI]; // updated in loop below

            // Calculates deviatoric, spherical and equivalent plastic strain using Newtons's method
            newtonLoop
            (
                DhydroEpsilonPI[cellI],
                DqEpsilonPI[cellI],
                fStarI[cellI],
                fStarOldI[cellI],
                curSigmaY ,      
                epsilonPEqOldI[cellI],
                DEpsilonPEqI[cellI],
                pTrialI[cellI],
                qTrialI[cellI]
            );

            pI[cellI] = pTrialI[cellI] - K_.value()*DhydroEpsilonPI[cellI];
            qI[cellI] = qTrialI[cellI] - 3*mu_.value()*DqEpsilonPI[cellI];
            sigmaI[cellI] = pI[cellI]*I + (2.0/3.0)*plasticNI[cellI]*qI[cellI];
            DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];

            // update porosity
            calculatef
            (
                DhydroEpsilonPI[cellI],
                fI[cellI],
                fOldI[cellI], 
                fStarI[cellI],
                fStarOldI[cellI],
                epsilonPEqOldI[cellI]+DEpsilonPEqI[cellI],
                DEpsilonPEqI[cellI]
            );
        }

    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
        const scalarField& sigmaYOldP = sigmaY_.oldTime().boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEq_.oldTime().boundaryField()[patchI];

        scalarField& pP=sigmaHyd_.boundaryField()[patchI];
        scalarField& qTrialP = qTrial.boundaryField()[patchI];
        symmTensorField& sigmaP = sigma.boundaryField()[patchI];
        symmTensorField& sP = s.boundaryField()[patchI];
        scalarField& qP = q.boundaryField()[patchI];
        scalarField& fOldP = f_.oldTime().boundaryField()[patchI];
        scalarField& fP = f_.boundaryField()[patchI];
        scalarField& fStarP = fStar_.boundaryField()[patchI];
        scalarField& fStarOldP = fStar_.oldTime().boundaryField()[patchI];
        scalarField& pTrialP = sigmaHydTr.boundaryField()[patchI];

        scalarField& DhydroEpsilonPP = DHydroEpsilonP_.boundaryField()[patchI];
        scalarField& DqEpsilonPP = DqEpsilonP_.boundaryField()[patchI];
        scalarField& DEpsilonPEqP = DEpsilonPEq_.boundaryField()[patchI];

        scalarField hydEpsP = hydroEpsilonTrial.boundaryField()[patchI];
        symmTensorField  epsilonElP = epsilonElTrial.boundaryField()[patchI];


        forAll(fTrialP, faceI)
        {

            // Calculate direction plasticN
            const scalar magS = mag(sTrialP[faceI]);
            if (magS > SMALL)
            {
                plasticNP[faceI] = (3*sTrialP[faceI])/(2*qTrialP[faceI]);
            }

            // Update variables when no plastic deformation
            if (fTrialP[faceI] < SMALL)
            {
                // elasticity
                DSigmaYP[faceI] = 0.0;
                DLambdaP[faceI] = 0.0;
                plasticNP[faceI] = symmTensor::zero;

		pP = pTrialP[faceI];
		fP[faceI] = fOldP[faceI];
		fStarP[faceI] = fStarOldP[faceI];
		sP[faceI] = sTrialP[faceI];
		qP[faceI] = qTrialP[faceI];
		sigmaP[faceI] = pP[faceI]*I+sTrialP[faceI];
		DEpsilonPEqP[faceI] = 0;
            }
            else
            {
                scalar curSigmaY = sigmaYOldP[faceI] ;// updated in loop below

                // Calculate DEpsilonPEq and curSigmaY
                newtonLoop
                (
                    DhydroEpsilonPP[faceI],
                    DqEpsilonPP[faceI],
                    fStarP[faceI],
                    fStarOldP[faceI],
                    curSigmaY ,      
                    epsilonPEqOldP[faceI],
                    DEpsilonPEqP[faceI],
                    pTrialP[faceI],
                    qTrialP[faceI]
                );

                pP = pTrialP[faceI] - K_.value()*DhydroEpsilonPP[faceI];
                qP[faceI] = qTrialP[faceI] - 3*mu_.value()*DqEpsilonPP[faceI];
                sigmaP[faceI] = pP[faceI]*I + (2.0/3.0)*plasticNP[faceI]*qP[faceI];
                symmTensor depsilonP = (1/3)*DhydroEpsilonPP[faceI]*I + DqEpsilonPP[faceI]*plasticNP[faceI];
                epsilonElP[faceI] = epsilonElP[faceI] - depsilonP;

                // calculate porosity
                calculatef
                (
                    DhydroEpsilonPP[faceI],
                    fP[faceI],
                    fOldP[faceI],
                    fStarP[faceI],
                    fStarOldP[faceI],
                    epsilonPEqOldP[faceI] + DEpsilonPEqP[faceI],
                    DEpsilonPEqP[faceI]
                );

                // Update increment of yield stress
                DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                
            }

        }
    }

    DEpsilonP_ = (1/3)*DHydroEpsilonP_*I + DqEpsilonP_*plasticN_;
    f_ = fStar_;

}


void Foam::linearGTN::correct(surfaceSymmTensorField& sigma)
{
   notImplemented("wip");
}


Foam::scalar Foam::linearGTN::residual()
{
    // Calculate residual based on change in plastic strain increment

        return
            gMax
            (
                mag
                (
                    DEpsilonP_.internalField()
                  - DEpsilonP_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().internalField()));
    
}


void Foam::linearGTN::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
    epsilonP_ += DEpsilonP_;
    hydroEpsilonP_ += DHydroEpsilonP_;
    qEpsilonP_ += DqEpsilonP_;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.internalField(), celli)
    {
        if (DEpsilonPEq_.internalField()[celli] > SMALL)
        {
            activeYield_.internalField()[celli] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.internalField()[celli] = 0.0;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchi)
    {
        if (!activeYield_.boundaryField()[patchi].coupled())
        {
            forAll(activeYield_.boundaryField()[patchi], facei)
            {
                if (DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
                {
                    activeYield_.boundaryField()[patchi][facei] = 1.0;
                }
                else
                {
                    activeYield_.boundaryField()[patchi][facei] = 0.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    const int nTotalCells = returnReduce(mesh().nCells(), sumOp<int>());

    Info<< "    " << numCellsYielding << " cells ("
        << 100.0*scalar(numCellsYielding)/scalar(nTotalCells)
        << "% of the cells in this material) are actively yielding"
        << nl << endl;

    // Testing
    if (mesh().time().outputTime())
    {
        volScalarField epsilonPMag
        (
            IOobject
            (
                "epsilonPMag",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((2.0/3.0)*magSqr(dev(epsilonP_)))
        );

        epsilonPMag.write();
    }
}


Foam::scalar Foam::linearGTN::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the Euler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon_)));

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    const symmTensorField& plasticNI = plasticN_.internalField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();

    // Calculate error field
    const symmTensorField DEpsilonPErrorI =
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
       /(epsilonEqI + SMALL);

    // Max error
    const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    if (maxMagDEpsilonPErr > SMALL)
    {
        Info<< "    " << name() << ": max time integration error = "
            << maxMagDEpsilonPErr
            << endl;

        if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
        {
            WarningIn
            (
                "Foam::scalar Foam::linearGTN::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is over 50 times larger "
                << "than the desired value!\n    Consider starting the "
                << "simulation with a smaller initial time-step" << endl;
        }

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


// ************************************************************************* //
