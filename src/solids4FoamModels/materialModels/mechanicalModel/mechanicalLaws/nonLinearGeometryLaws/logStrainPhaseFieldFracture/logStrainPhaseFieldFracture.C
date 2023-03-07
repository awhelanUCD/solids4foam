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

#include "logStrainPhaseFieldFracture.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(logStrainPhaseFieldFracture, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, logStrainPhaseFieldFracture, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar logStrainPhaseFieldFracture::LoopTol_ = 1e-4;

    // Maximum number of iterations for Newton loop
    label logStrainPhaseFieldFracture::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar logStrainPhaseFieldFracture::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar logStrainPhaseFieldFracture::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::symmTensor Foam::logStrainPhaseFieldFracture::logm(const symmTensor& T)
{
    // Finds the matrix log of the tensor T

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    Eigen::Vector3d LogEVals;
    LogEVals(0) = Foam::log(EVals(0));
    LogEVals(1) = Foam::log(EVals(1));
    LogEVals(2) = Foam::log(EVals(2));

    Eigen::Matrix3d D = LogEVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Eigen::Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}
Foam::symmTensor Foam::logStrainPhaseFieldFracture::expm(const symmTensor& T)
{
    // Calculate exponential of a tensor

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    Eigen::Vector3d ExpEVals;
    ExpEVals(0) = Foam::exp(EVals(0));
    ExpEVals(1) = Foam::exp(EVals(1));
    ExpEVals(2) = Foam::exp(EVals(2));

    Eigen::Matrix3d D = ExpEVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Eigen::Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}
void Foam::logStrainPhaseFieldFracture::makeRelF()
{
    if (relFPtr_)
    {
        FatalErrorIn("void Foam::logStrainPhaseFieldFracture::makeRelF()")
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


Foam::volTensorField& Foam::logStrainPhaseFieldFracture::relF()
{
    if (!relFPtr_)
    {
        makeRelF();
    }

    return *relFPtr_;
}

void Foam::logStrainPhaseFieldFracture::makeJ()
{
    if (JPtr_)
    {
        FatalErrorIn("void Foam::logStrainPhaseFieldFracture::makeJ()")
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
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );

    // Store the old-time
    JPtr_->oldTime();
}


Foam::volScalarField& Foam::logStrainPhaseFieldFracture::J()
{
    if (!JPtr_)
    {
        makeJ();
    }

    return *JPtr_;
}

void Foam::logStrainPhaseFieldFracture::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::logStrainPhaseFieldFracture::makeF()")
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


Foam::volTensorField& Foam::logStrainPhaseFieldFracture::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}

void Foam::logStrainPhaseFieldFracture::smallStrainReturnMap
(
    const symmTensor& Ee,  //elastic strain	            
    scalar& sigmaY,  //yield stress
    scalar& DSigmaY, // increment of yield stress         
    scalar& DLambda, //plastic strain increment                               
    const scalar& epsilonPEqOld,  //old value of the equivalent plastic strain
    scalar& activeYield, //defined if cell is undergoing active yielding
    const scalar& gc //degredation fucntion
) const
{
    sigmaY = stressPlasticStrainSeries_(epsilonPEqOld);
    scalar qTrial=sqrt(3.0/2.0)*mag(gc*2*mu_.value()*(dev(Ee)));
    scalar fTrial=qTrial-gc*sigmaY;
    DLambda=0;

    int i=0;
    if (fTrial>0.0)
    {
        //newton loop to calculate plastic strain increment
        do
        {
            fTrial = qTrial - gc*3*mu_.value()*DLambda 
                            - gc*(stressPlasticStrainSeries_(epsilonPEqOld+DLambda));

            scalar dSigmaY = (stressPlasticStrainSeries_(epsilonPEqOld+DLambda+1e-8)
                             -stressPlasticStrainSeries_(epsilonPEqOld+DLambda))/1e-8;

            scalar dfTrial = -gc*3*mu_.value() - gc*dSigmaY;
            DLambda -= fTrial/dfTrial; 

        if (i == MaxNewtonIter_)
        {
            WarningIn("neoHookeanElasticMisesPlasticDamage::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(fTrial) > LoopTol_) && ++i < MaxNewtonIter_);

    activeYield = 1.0;

    }
    else
    {
        DLambda = 0;
        activeYield = 0.0;
    }

    sigmaY = stressPlasticStrainSeries_(epsilonPEqOld+DLambda);
}


Foam::tmp<Foam::volScalarField> Foam::logStrainPhaseFieldFracture::Ibar
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


Foam::tmp<Foam::surfaceScalarField> Foam::logStrainPhaseFieldFracture::Ibar
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
Foam::logStrainPhaseFieldFracture::logStrainPhaseFieldFracture
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    charLength_("zero", dimLength, 0.0),
    Gc_("zero", dimPressure, 0.0),
    d1_("zero", dimless, 0.0),
    d2_("zero", dimless, 0.0),
    d3_("zero", dimless, 0.0),
    betaE_(dict.lookupOrDefault<scalar>("betaE", 1.0)),
    betaP_(dict.lookupOrDefault<scalar>("betaP", 1.0)),
    w0_("zero", dimPressure, 0.0),
    monolithic_
    (
        dict.lookupOrDefault<Switch>("monolithic", true)
    ),
    triaxiality_
    (
        dict.lookupOrDefault<Switch>("triaxiality", false)
    ),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    debug_(dict.lookup("debug")),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    relFPtr_(NULL),
    JPtr_(NULL),
    FPtr_(NULL),
    stressPlasticStrainSeries_(dict),
    Dpf_
    (
        IOobject
        (
            "Dpf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("1", dimless, 1)
    ),
    one_
    (
        IOobject
        (
            "one",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    ),
    zero_
    (
        IOobject
        (
            "zero",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    gradDpf_
    (
        IOobject
        (
            "grad(Dpf)",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless/dimLength, vector::zero)
    ),
    strainEnergy_
    (
        IOobject
        (
            "strainEnergy",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0)
    ),
    HStrainEnergy_
    (
        IOobject
        (
            "HStrainEnergy",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0)
    ),
    plasticStrainEnergy_
    (
        IOobject
        (
            "plasticStrainEnergy",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0)
    ),
    devStrainEnergy_
    (
        IOobject
        (
            "devStrainEnergy",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0)
    ),
    volStrainEnergy_
    (
        IOobject
        (
            "volStrainEnergy",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0)
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
    tau_
    (
        IOobject
        (
            "tau",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
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
            "initialYieldStress", dimPressure, 0.0
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
    epsilonEl_
    (
        IOobject
        (
            "epsilonEl",
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
            IOobject::AUTO_WRITE
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
    //nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    ),
  // Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.oldTime();
    bEbar_.oldTime();
    epsilonEl_.oldTime();
    HStrainEnergy_.oldTime();
    plasticStrainEnergy_.oldTime();
  
    if (triaxiality_)
    {
        d1_ = dimensionedScalar(dict.lookup("d1"));
        d2_ = dimensionedScalar(dict.lookup("d2"));
        d3_ = dimensionedScalar(dict.lookup("d3"));
    }

    charLength_ = dimensionedScalar(dict.lookup("charLength"));
    Gc_ = dimensionedScalar(dict.lookup("Gc"));
    w0_ = dimensionedScalar(dict.lookup("w0"));

    Dpf_=1.0;

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
            "logStrainPhaseFieldFracture::logStrainPhaseFieldFracture::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }



    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::logStrainPhaseFieldFracture::~logStrainPhaseFieldFracture()
{
    deleteDemandDrivenData(relFPtr_);
    deleteDemandDrivenData(JPtr_);
    deleteDemandDrivenData(FPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::logStrainPhaseFieldFracture::rho() const
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
Foam::logStrainPhaseFieldFracture::impK() const
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

void Foam::logStrainPhaseFieldFracture::calcPhase()
{

    if (!monolithic_) 
    {
        Info << "Updating Phase Field" << endl;
    }

    //set deviatoric elastic strain
    volSymmTensorField devEpsilonEl = dev(epsilonEl_);
     
    //calculate plastic strain energy
    plasticStrainEnergy_ = plasticStrainEnergy_.oldTime() + sigmaY_*DLambda_;

    //initialise effective plastic strain energy
    volScalarField effPlasticStrainEnergy = plasticStrainEnergy_;

    //loop to incorporate triaxiality effects
    if (triaxiality_)
    {
        volScalarField magS = mag(pow(one_-Dpf_, 2.0)*2.0*mu_*devEpsilonEl);
        volScalarField triaxParam = d1_ + d2_*pow(2.71828, d3_*sigmaHyd_/magS);
        plasticStrainEnergy_ = plasticStrainEnergy_.oldTime() + sigmaY_*DLambda_/triaxParam;
    }

    //set devatoric and volumetric elastic strain energy
    volStrainEnergy_ = 0.5*K_*pow(tr(epsilonEl_),2.0);
    devStrainEnergy_ = mu_*(devEpsilonEl && devEpsilonEl);

    //initialise internal fields
    symmTensorField& epsilonElI = epsilonEl_.internalField();
    scalarField& strainEnergyI = strainEnergy_.internalField();
    scalarField volStrainEnergyI = volStrainEnergy_.internalField();
    scalarField devStrainEnergyI = devStrainEnergy_.internalField();
    scalarField& HStrainEnergyI = HStrainEnergy_.internalField();
    scalarField& plasticStrainEnergyI = plasticStrainEnergy_.internalField();
    scalarField& effPlasticStrainEnergyI = effPlasticStrainEnergy.internalField();
    scalarField& HStrainEnergyOldI = HStrainEnergy_.oldTime().internalField();

    //set history strain energy field and effective plastic strain energy
    forAll(HStrainEnergyI, cellI)
    { 
        if (tr(epsilonElI[cellI]) >= 0.0)
        {
            strainEnergyI[cellI] = volStrainEnergyI[cellI] + devStrainEnergyI[cellI];
        }
        else 
        {
            strainEnergyI[cellI] = devStrainEnergyI[cellI];
        }

        effPlasticStrainEnergyI[cellI] = effPlasticStrainEnergyI[cellI] - w0_.value();

        if (plasticStrainEnergyI[cellI] < w0_.value())
        {
            effPlasticStrainEnergyI[cellI] = 0;
        }

        HStrainEnergyI[cellI] = max(strainEnergyI[cellI],HStrainEnergyOldI[cellI]);
    }

    //set history strain energy field and effective plastic strain energy
    forAll(HStrainEnergy_.boundaryField(), patchI)
    {
        //initialise boundary fields
        symmTensorField& epsilonElP = epsilonEl_.boundaryField()[patchI];
        scalarField& strainEnergyP = strainEnergy_.boundaryField()[patchI];
        scalarField volStrainEnergyP = volStrainEnergy_.boundaryField()[patchI];
        scalarField devStrainEnergyP = devStrainEnergy_.boundaryField()[patchI];
        scalarField& HStrainEnergyP = HStrainEnergy_.boundaryField()[patchI];
        scalarField& plasticStrainEnergyP = plasticStrainEnergy_.boundaryField()[patchI];
        scalarField& effPlasticStrainEnergyP = effPlasticStrainEnergy.boundaryField()[patchI];
        scalarField& HStrainEnergyOldP = HStrainEnergy_.oldTime().boundaryField()[patchI];

        forAll(HStrainEnergyP, faceI)
        {

            if (tr(epsilonElP[faceI]) >= 0.0)
            {
                strainEnergyP[faceI] = volStrainEnergyP[faceI] + devStrainEnergyP[faceI];
            }
            else 
            {
                strainEnergyP[faceI] = devStrainEnergyP[faceI];
            }

            effPlasticStrainEnergyP[faceI] = effPlasticStrainEnergyP[faceI] - w0_.value();

            if (plasticStrainEnergyP[faceI] < w0_.value())
            {
                effPlasticStrainEnergyP[faceI] = 0;
            }

            HStrainEnergyP[faceI] = max(strainEnergyP[faceI],HStrainEnergyOldP[faceI]);
        }
    }

    volScalarField elasticCoeff = (4.0*charLength_*HStrainEnergy_)/(Gc_);
    volScalarField plasticCoeff = (4.0*charLength_*effPlasticStrainEnergy)/(Gc_);
    volScalarField coeff = elasticCoeff + plasticCoeff;
    dimensionedScalar diml = charLength_/charLength_.value(); //this is used to ensure OpenFOAM doesn't give a dimension error

    //calculates phase field variable 
    fvScalarMatrix DpfEqn
    (
        -fvm::laplacian(4.0*pow(charLength_,2.0), Dpf_)
        +fvm::Sp(1.0, Dpf_)
        +fvm::Sp(coeff/diml, Dpf_)
        ==   
        one_
    );

    DpfEqn.solve();

    gradDpf_ = fvc::grad(Dpf_);
  

}

void Foam::logStrainPhaseFieldFracture::correct(volSymmTensorField& sigma)
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
    DEpsilonP_.storePrevIter();

    volScalarField gc_ = pow(Dpf_, 2.0);//set crack degredation function

    // Take references to the internal fields
    const tensorField& FI = F().internalField();
    symmTensorField epsilonElOldI = epsilonEl_.oldTime().internalField();
    symmTensorField& epsilonElI = epsilonEl_.internalField();
    symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    symmTensorField& tauI = tau_.internalField(); 
    symmTensorField& sigmaI = sigma.internalField();
    tensorField  relFI = relF().internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& pI = sigmaHyd_.internalField();
    scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().internalField();
    scalarField& DEpsilonPEqI = DEpsilonPEq_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    scalarField& sigmaYI = sigmaY_.internalField();
    scalarField& activeYieldI = activeYield_.internalField();
    scalarField& gcI = gc_.internalField();

    //initialise variables for plasticity model
    symmTensor BeOld = I; //old value of the left Cachy-Green tensor
    symmTensor Be = I; //left Cachy-Green tensor
    symmTensor Ee = symmTensor::zero;  //elastic strain tensor
    tensor relf=I; //relative deformation gradient
    symmTensor plasticN=symmTensor::zero; //normal to the plastic yielding

    // Calculate DLambda and plasticN
    forAll(relFI, cellI)
    {
        //pre-processing step
        relf = relFI[cellI];
        BeOld = expm(2.0*epsilonElOldI[cellI]);
        Be = symm(relf & BeOld & relf.T());
        Ee = 0.5*logm(Be);

        smallStrainReturnMap
        (
            Ee,
            sigmaYI[cellI],
            DSigmaYI[cellI],
            DLambdaI[cellI],
            epsilonPEqOldI[cellI],
            activeYieldI[cellI],
            gcI[cellI]
        );

        if (mag(Ee) > 0)
        {
            plasticN = dev(Ee)/(mag(dev(Ee)));
        }

        DEpsilonPI[cellI] = sqrt(3.0/2.0)*DLambdaI[cellI]*plasticN ;
        epsilonElI[cellI] = Ee-DEpsilonPI[cellI];

        if (tr(epsilonElI[cellI]) >= 0)
        {   
            pI[cellI] = gcI[cellI]*K_.value()*tr(epsilonElI[cellI]);
            tauI[cellI] = gcI[cellI]*2*mu_.value()*dev(epsilonElI[cellI]) + pI[cellI]*I;
        }
        else 
        {    
            pI[cellI] = K_.value()*tr(epsilonElI[cellI]);
            tauI[cellI] = gcI[cellI]*2*mu_.value()*dev(epsilonElI[cellI]) + pI[cellI]*I;
        }

        sigmaI[cellI] = (1/det(FI[cellI]))*tauI[cellI];
        DEpsilonPEqI[cellI] = DLambdaI[cellI];


    }

    forAll(F().boundaryField(), patchI)
    {
        // Take references to boundary fields
        const tensorField& FP = F().boundaryField()[patchI];
        symmTensorField epsilonElOldP = epsilonEl_.oldTime().boundaryField()[patchI];
        symmTensorField& epsilonElP = epsilonEl_.boundaryField()[patchI];
        symmTensorField& DEpsilonPP = DEpsilonP_.boundaryField()[patchI];
        symmTensorField& tauP = tau_.boundaryField()[patchI]; 
        symmTensorField& sigmaP = sigma.boundaryField()[patchI];
        tensorField  relFP = relF().boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& pP = sigmaHyd_.boundaryField()[patchI];
        scalarField& epsilonPEqOldP = epsilonPEq_.oldTime().boundaryField()[patchI];
        scalarField& DEpsilonPEqP = DEpsilonPEq_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
        scalarField& activeYieldP = activeYield_.boundaryField()[patchI];
        scalarField& gcP = gc_.boundaryField()[patchI];

        // Calculate DLambda and plasticN
        forAll(FP, faceI)
        {
            //pre-processing step
            relf = relFP[faceI];
            BeOld = expm(2.0*epsilonElOldP[faceI]);
            Be = symm(relf & BeOld & relf.T());
            Ee = 0.5*logm(Be);

            smallStrainReturnMap
            (
                Ee,
                sigmaYP[faceI],
                DSigmaYP[faceI],
                DLambdaP[faceI],
                epsilonPEqOldP[faceI],
                activeYieldP[faceI],
                gcP[faceI]
            );

            //post processing step
            if (mag(Ee)>0)
            {
                plasticN = dev(Ee)/(mag(dev(Ee)));
            }

            DEpsilonPP[faceI] = sqrt(3.0/2.0)*DLambdaP[faceI]*plasticN ;
            epsilonElP[faceI] = Ee-DEpsilonPP[faceI];

            if (tr(epsilonElP[faceI]) >= 0)
            {   
                pP[faceI] = gcP[faceI]*K_.value()*tr(epsilonElP[faceI]);
                tauP[faceI] = gcP[faceI]*2*mu_.value()*dev(epsilonElP[faceI]) + pP[faceI]*I;
            }
            else 
            {    
                pP[faceI] = K_.value()*tr(epsilonElP[faceI]);
                tauP[faceI] = gcP[faceI]*2*mu_.value()*dev(epsilonElP[faceI]) + pP[faceI]*I;
            }

            sigmaP[faceI] = (1/det(FP[faceI]))* tauP[faceI];
            DEpsilonPEqP[faceI] = DLambdaP[faceI];

        }
    }

    if (monolithic_)
    {
        calcPhase();
    }


}


void Foam::logStrainPhaseFieldFracture::correct(surfaceSymmTensorField& sigma)
{

   notImplemented("wip");
}


Foam::scalar Foam::logStrainPhaseFieldFracture::residual()
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


void Foam::logStrainPhaseFieldFracture::updateTotalFields()
{
    if (!monolithic_)
    {
        calcPhase();
    }

    Info<< nl << "Updating total accumulated fields" << endl;
    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
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


Foam::scalar Foam::logStrainPhaseFieldFracture::newDeltaT()
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
                "Foam::scalar Foam::logStrainPhaseFieldFracture::newDeltaT()"
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
