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

#include "GreenStrainLemaitreDamage.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GreenStrainLemaitreDamage, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, GreenStrainLemaitreDamage, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar GreenStrainLemaitreDamage::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label GreenStrainLemaitreDamage::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar GreenStrainLemaitreDamage::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar GreenStrainLemaitreDamage::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GreenStrainLemaitreDamage::makeRelF()
{
    if (relFPtr_)
    {
        FatalErrorIn("void Foam::GreenStrainLemaitreDamage::makeRelF()")
            << "pointer already set" << abort(FatalError);
    }

    relFPtr_ =
        new volTensorField
        (
            IOobject
            (
                "lawRelF",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::volTensorField& Foam::GreenStrainLemaitreDamage::relF()
{
    if (!relFPtr_)
    {
        makeRelF();
    }

    return *relFPtr_;
}

void Foam::GreenStrainLemaitreDamage::makeJ()
{
    if (JPtr_)
    {
        FatalErrorIn("void Foam::GreenStrainLemaitreDamage::makeJ()")
            << "pointer already set" << abort(FatalError);
    }

    JPtr_ =
        new volScalarField
        (
            IOobject
            (
                "lawJ",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );

    // Store the old-time
    JPtr_->oldTime();
}


Foam::volScalarField& Foam::GreenStrainLemaitreDamage::J()
{
    if (!JPtr_)
    {
        makeJ();
    }

    return *JPtr_;
}

void Foam::GreenStrainLemaitreDamage::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::GreenStrainLemaitreDamage::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_ =
        new volTensorField
        (
            IOobject
            (
                "lawF",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );

    // Store the old-time
    FPtr_->oldTime();
}


Foam::volTensorField& Foam::GreenStrainLemaitreDamage::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}

void Foam::GreenStrainLemaitreDamage::damageFunction
(
    const scalar triaxiality, // Triaxiality   
    const scalar epsilonPEq, // Equivalent plastic strain       
    const scalar sigmaEq, //Equivalent stress     
    scalar& damage, // Damage
    const scalar damageOld, // Old damage value
    const scalar& damageNonLocal, // Non-local damage 
    const scalar DLambda // Plastic multiplier increment
) const
{
    // Initialise energy release rate Y
    scalar Y;

    scalar Ystar;

    if ( triaxiality> (-1.0/3.0)&& epsilonPEq > epsilonD_.value())
    {
        // Denominator for damage law
        const scalar denom = (2.0/3.0)*(1.0 + nu_.value()) + 3.0*(1.0 - 2.0* nu_.value())*pow(triaxiality, 2.0);

        // Energy release rate
        Y = -(pow(sigmaEq, 2.0)/(2.0*E_.value()))*denom;

        // Update damage
        Ystar = pow(-Y/s0_.value(), b_.value())*DLambda*sqrtTwoOverThree_;
        damage = damageOld + Ystar/(1-damageNonLocal);

        if (damage > damageC_) {damage = 0.99;}
    }
}

Foam::scalar Foam::GreenStrainLemaitreDamage::curYieldStress
(
    const scalar plasticMultOld,    // Current equivalent plastic strain
    const scalar DLambda,
    const scalar J                 // Current Jacobian
) const
{
    // We assume that the stress-strain curve was specifed as Cauchy stress vs
    // true strain, but we want the Kirchhoff (tau) yield stress,
    // so we multiply Cauchy stress by J as tauSigmaY = J*sigmaCauchySigmaY

    return J*(stressPlasticStrainSeries_(plasticMultOld + sqrtTwoOverThree_*DLambda));
}


Foam::scalar Foam::GreenStrainLemaitreDamage::yieldFunction
(
    const scalar plasticMultOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar,            // Scaled shear modulus
    const scalar J ,                // Current Jacobian
    const scalar damage
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN/(1-damage)) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda/(1-damage) - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current Kirchhoff yield stress which is typically a
    // function of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
       magSTrial - 2*muBar*DLambda/(1-damage)
      - 
          sqrtTwoOverThree_* curYieldStress
            (
                plasticMultOld , DLambda,
                J
            );
}


void Foam::GreenStrainLemaitreDamage::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar plasticMultOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar J,                // Current Jacobian
    const scalar maxMagDEpsilon,    // Max strain increment magnitude
    const scalar damage // Damage (non-local)
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction(plasticMultOld, magSTrial, DLambda, muBar, J, damage);
    scalar residual = 1.0;

if (debug_) {Info<<"ftrial: "<<fTrial<<endl;}
    do
    {

if (debug_) {Info <<"i:" <<i<<endl;}
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluations are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                plasticMultOld, magSTrial, DLambda + finiteDiff_, muBar, J, damage
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;

        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(plasticMultOld, magSTrial, DLambda, muBar,  J, damage);

        if (i == MaxNewtonIter_)
        {
            WarningIn("GreenStrainLemaitreDamage::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Note: we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY =
        curYieldStress
        (
            plasticMultOld,DLambda,  J
        )/J;
}


Foam::tmp<Foam::volScalarField> Foam::GreenStrainLemaitreDamage::Ibar
(
    const volSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<volScalarField> tIbar
    (
        new volScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& Ibar = tIbar();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
            scalarField& IbarP = Ibar.boundaryField()[patchI];
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}


Foam::tmp<Foam::surfaceScalarField> Foam::GreenStrainLemaitreDamage::Ibar
(
    const surfaceSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<surfaceScalarField> tIbar
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    surfaceScalarField& Ibar = tIbar();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
            scalarField& IbarP = Ibar.boundaryField()[patchI];
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::GreenStrainLemaitreDamage::GreenStrainLemaitreDamage
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
//    y0_(readScalar(dict.lookup("y0"))),
//    hIso_(readScalar(dict.lookup("hIso"))),
//    beta_(readScalar(dict.lookup("y0"))),
//    yInf_(readScalar(dict.lookup("hIso"))),
    charLength_("zero", dimLength, 0.0),
    epsilonD_("zero", dimless, 0.0),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    b_("zero", dimless, 0.0),
    s0_("zero", dimPressure, 0.0),
    damageC_(dict.lookupOrDefault<scalar>("damageC", 0.99)),
    debug_(dict.lookup("debug")),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    relFPtr_(NULL),
    JPtr_(NULL),
    FPtr_(NULL),
    stressPlasticStrainSeries_(dict),
    plasticMult_
    (
        IOobject
        (
            "plasticMult",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    damage_
    (
        IOobject
        (
            "damage",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    damageNonLocal_
    (
        IOobject
        (
            "damageNonLocal",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    gradDamageNonLocal_
    (
        IOobject
        (
            "grad(damageNonLocal)",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless/dimLength, vector::zero)
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
    smoothPressure_(dict.lookupOrDefault<Switch>("smoothPressure", true)),
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
    bEbarTrial_
    (
        IOobject
        (
            "bEbarTrial",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
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
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    ),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.oldTime();
    bEbar_.oldTime();
    damage_.oldTime();
    plasticMult_.oldTime();

    // Material parameters associated with damage evolution
    s0_ = dimensionedScalar(dict.lookup("s0"));
    b_ = dimensionedScalar(dict.lookup("b"));
    epsilonD_ = dimensionedScalar(dict.lookup("epsilonD"));
    charLength_ = dimensionedScalar(dict.lookup("charLength"));

    Info<< "    smoothPressure: " << smoothPressure_ << nl
        << "    updateBEbarConsistent: " << updateBEbarConsistent_ << endl;

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
        if (planeStress())
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));
    }
    else
    {
        FatalErrorIn
        (
            "GreenStrainLemaitreDamage::GreenStrainLemaitreDamage::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GreenStrainLemaitreDamage::~GreenStrainLemaitreDamage()
{
    deleteDemandDrivenData(relFPtr_);
    deleteDemandDrivenData(JPtr_);
    deleteDemandDrivenData(FPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GreenStrainLemaitreDamage::rho() const
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
Foam::GreenStrainLemaitreDamage::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = mu_*dev(bEbarTrial_);

    const volScalarField Ibar = tr(bEbarTrial_)/3.0;
    const volScalarField muBar = Ibar*mu_;

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial =
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL));

    // Calculate scaling factor
    const volScalarField scaleFactor = 1.0 - (2.0*muBar*DLambda_/magSTrial);

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


void Foam::GreenStrainLemaitreDamage::correct(volSymmTensorField& sigma)
{
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the relative deformation gradient
        relF() = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF() & F().oldTime();
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            F() = F().oldTime() + gradDD.T();

            // Update the relative deformation gradient
            relF() = F() & inv(F().oldTime());
        }
        else
        {
            // Lookup gradient of displacement
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");

            // Update the total deformation gradient
            F() = I + gradD.T();

            // Update the relative deformation gradient
            relF() = F() & inv(F().oldTime());
        }
    }
    else
    {
        FatalErrorIn
        (
            type() + "::correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // Update the Jacobian of the total deformation gradient
    J() = det(F());

    // Calculate the relative Jacobian
    const volScalarField relJ = J()/J().oldTime();

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    const volTensorField relFbar = pow(relJ, -1.0/3.0)*relF();

    // Update bE trial
    bEbarTrial_ = transform(relFbar, bEbar_.oldTime());

    // Calculate trial deviatoric stress
    const volSymmTensorField sTrial = mu_*dev(bEbarTrial_);

    const volScalarField Ibar = tr(bEbarTrial_)/3.0;

    const volScalarField muBar = Ibar*mu_;

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonP_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(bEbarTrial_.internalField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const volScalarField fTrial = mag(sTrial) - sqrtTwoOverThree_*J()*sigmaY_;

    // Take references to the internal fields for efficiency
    // Note: these references are also used to calculate the damage 
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    const scalarField& JI = J().internalField();
    const scalarField& sigmaYI = sigmaY_.internalField();

    const scalarField& plasticMultOldI = plasticMult_.oldTime().internalField();
    scalarField& damageI = damage_.internalField();
    const scalarField& damageOldI = damage_.oldTime().internalField();
    const scalarField& damageNonLocalI = damageNonLocal_.internalField();

    // Calculate DLambda and plasticN
    forAll(fTrialI, cellI)
    {

        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[cellI]);

        if (magS > SMALL)
        {
            plasticNI[cellI] = sTrialI[cellI]/(magS*(1-damageNonLocalI[cellI]));
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[cellI] < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
            damageI[cellI] = damageOldI[cellI];

        }
        else
        {
                 // Total equivalent plastic strain where t is start of time-step
                scalar curSigmaY = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambdaI[cellI],
                    curSigmaY,
                    plasticMultOldI[cellI],
                    magS,
                    mu_.value(),
                    JI[cellI],
                    maxMagBE,
                    damageNonLocalI[cellI]
                );

                // Update increment of yield stress
                DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
            
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
        const scalarField& JP = J().boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];

        const scalarField& plasticMultOldP = plasticMult_.oldTime().boundaryField()[patchI];
        scalarField& damageP = damage_.boundaryField()[patchI];
        const scalarField& damageOldP = damage_.oldTime().boundaryField()[patchI];
        const scalarField& damageNonLocalP = damageNonLocal_.boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {

            // Calculate direction plasticN
            const scalar magS = mag(sTrialP[faceI]);
            if (magS > SMALL)
            {
                plasticNP[faceI] = sTrialP[faceI]/(magS*(1-damageNonLocalP[faceI]));
            }

            // Calculate DEpsilonPEq
            if (fTrialP[faceI] < SMALL)
            {
                // elasticity
                DSigmaYP[faceI] = 0.0;
                DLambdaP[faceI] = 0.0;
                damageP[faceI]=damageOldP[faceI];
            }
            else
            {
                    scalar curSigmaY = 0.0; // updated in loop below

                    // Calculate DEpsilonPEq and curSigmaY
                    newtonLoop
                        (
                            DLambdaP[faceI],
                            curSigmaY,
                            plasticMultOldP[faceI],
                            magS,
                            mu_.value(),
                            JP[faceI],
                            maxMagBE,
                            damageNonLocalP[faceI]
                        );

                    // Update increment of yield stress
                    DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
            }
        }
    }

    // Update accumulated plastic multiplier
    plasticMult_=plasticMult_.oldTime()+sqrtTwoOverThree_*DLambda_;

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_/(1-damageNonLocal_);
    DEpsilonP_ = Ibar*DLambda_*plasticN_;
    DEpsilonP_.relax();

    epsilonPEq_ = epsilonPEq_.oldTime() + DEpsilonPEq_;

    // Calculate deviatoric stress
    const volSymmTensorField s = (1-damageNonLocal_)*(sTrial - 2*mu_*DEpsilonP_);


    // Update bEbar
    if (updateBEbarConsistent_)
    {
        const volSymmTensorField devBEbar = (s/mu_)/(1-damageNonLocal_);
        bEbar_ = devBEbar + this->Ibar(devBEbar)*I;
    }
    else
    {
        bEbar_ = (s/mu_)/(1-damageNonLocal_) + Ibar*I;
    }

    // Update hydrostatic stress (negative of pressure)
    if (smoothPressure_)
    {
        // Calculate the hydrostatic pressure by solving a Laplace equation;
        // this ensures smoothness of the field and quells oscillations

        // Lookup the momentum equation inverse diagonal field
        const volScalarField AD = mesh().lookupObject<volScalarField>("DEqnA");

        // Pressure diffusivity field
        // Note: (4.0/3.0)*mu + K == 2*mu + lambda
        const surfaceScalarField rDAf
        (
            "rDAf",
            fvc::interpolate
            (
                ((4.0/3.0)*mu_ + K_)/AD, "interpolate(grad(sigmaHyd))"
            )
        );

        const dimensionedScalar one("one", dimless, 1.0);
        //const dimensionedScalar fac(dict().lookup("smoothFactor"));

        // Construct the pressure equation
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(one, sigmaHyd_)
          - fvm::laplacian(rDAf, sigmaHyd_, "laplacian(DDD,DD)")
         ==
            0.5*K_*(pow(J(), 2.0) - 1.0)
          - fvc::div(rDAf*fvc::interpolate(fvc::grad(sigmaHyd_)) & mesh().Sf())
        );

        // Store the pressue field to allow under-relaxation
        sigmaHyd_.storePrevIter();

        // Under-relax the linear system
        sigmaHydEqn.relax(0.7);

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Under-relax the pressure field
        sigmaHyd_.relax(0.2);
    }
    else
    {
        // Calculate the hydrostatic pressure directly from the displacement
        // field
        sigmaHyd_ = 0.5*K_*(pow(J(), 2.0) - 1.0);
    }

    sigmaHyd_=(1-damageNonLocal_)*sigmaHyd_;

    // Update the Cauchy stress
    sigma = (1.0/J())*(sigmaHyd_*I + s);



// Calculate damage field


    // Set fields needed to calculate damage
    const volSymmTensorField tau = J()*sigma;
    const volScalarField sigmaEq =
            sqrt(3.0/2.0)*mag(dev(tau/(1-damageNonLocal_)));
    const volScalarField  triaxiality = sigmaHyd_/(sqrt(3.0/2.0)*mag(dev(tau)));//sigmaHyd_(1.0/3.0)*tr(tau)

    // Take refernces to internal fields
    const scalarField& sigmaEqI = sigmaEq.internalField();   
    const scalarField& triaxialityI = triaxiality.internalField();
    const scalarField& DepsilonPeqI = DEpsilonPEq_.internalField();
    const scalarField& epsilonPEqI = epsilonPEq_.internalField();

    forAll(sigmaEqI, cellI)
    {
        damageFunction
        (
            triaxialityI[cellI],
            epsilonPEqI[cellI],
            sigmaEqI[cellI],
            damageI[cellI],
            damageOldI[cellI],
            damageNonLocalI[cellI],
            DLambdaI[cellI]
        ); 
    }
   
     // Take refernces to boundary fields
     forAll(damage_.boundaryField(), patchI)
     {

         const scalarField& sigmaEqP = sigmaEq.boundaryField()[patchI];  
         const scalarField& triaxialityP = triaxiality.boundaryField()[patchI];
         const scalarField& damageNonLocalP = damageNonLocal_.boundaryField()[patchI];
         scalarField& damageP = damage_.boundaryField()[patchI];
         const scalarField& damageOldP = damage_.oldTime().boundaryField()[patchI];
         const scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
         const scalarField plasticMultOldP = plasticMult_.boundaryField()[patchI];
         const scalarField epsilonPEqP = epsilonPEq_.boundaryField()[patchI];

         forAll(damageNonLocalP, faceI)
         {
             damageFunction
             (
                 triaxialityP[faceI],
                 epsilonPEqP[faceI],
                 sigmaEqP[faceI],
                 damageP[faceI],
                 damageOldP[faceI],
                 damageNonLocalP[faceI],
                 DLambdaP[faceI]
             ); 

         }
     }

    // Calculate non-local damage
    fvScalarMatrix damageEqn
    {
        fvm::Sp(1.0, damageNonLocal_)
      - fvm::laplacian(pow(charLength_, 2.0), damageNonLocal_)
     == damage_
    };

    damageEqn.solve();

    gradDamageNonLocal_ = fvc::grad(damageNonLocal_);

}


void Foam::GreenStrainLemaitreDamage::correct(surfaceSymmTensorField& sigma)
{

   notImplemented("wip");

}


Foam::scalar Foam::GreenStrainLemaitreDamage::residual()
{
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


void Foam::GreenStrainLemaitreDamage::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    //epsilonPEq_ += DEpsilonPEq_;
   // epsilonPEqf_ += DEpsilonPEqf_;
    epsilonP_ += DEpsilonP_;

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

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


Foam::scalar Foam::GreenStrainLemaitreDamage::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Update the total deformatio gradient
    F() = relF() & F().oldTime();

    // Calculate the total true (Hencky) strain
    const volSymmTensorField epsilon = 0.5*log(symm(F().T() & F()));

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));

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
                "Foam::scalar Foam::GreenStrainLemaitreDamage::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is lover 50 times larger "
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
