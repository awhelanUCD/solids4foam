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

#include "linearLemaitreDamage.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearLemaitreDamage, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearLemaitreDamage, linGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar linearLemaitreDamage::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label linearLemaitreDamage::MaxNewtonIter_ = 1000;

    // finiteDiff is the delta for finite difference differentiation
    scalar linearLemaitreDamage::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar linearLemaitreDamage::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);


} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearLemaitreDamage::damageFunction
(
    const scalar triaxiality, // triaxiality value  
    const scalar epsilonPEq, // equivalent plastic strain
    const scalar sigmaEq, // equivalent stress     
    const scalar damageOld, // old damage value
    scalar& damage, // damage
    const scalar DepsilonPeq, // equivalent plastic strain increment
    const scalar DLambda, // plastic multiplier increment
    const scalar sigmaHyd // hydrostatic stress
) const
{
    // initialise total energy release rate
    scalar Y;
    
    scalar Ystar;

    if (triaxiality> (-1.0/3.0) && epsilonPEq > epsilonD_.value())
    {
        // denominator for lemaitre damage evolution equation 
        // (not used here as an alternative mathamatically equivalent way of calculating Y is used)
        // const scalar denom = (2.0/3.0)*(1.0 + nu_.value())+ 3.0*(1.0 - 2.0* nu_.value())*pow(triaxiality, 2.0);

        // energy release rate Y
        Y = -(pow(sigmaEq, 2.0)/(6*mu_.value()) + (pow(sigmaHyd,2.0))/(2*K_.value()));

        Ystar = pow(-Y/s0_.value(), b_.value())*DLambda;

        // update damage
        damage = damageOld + Ystar/(1-damage);

        if (damage > damageC_) {damage = 0.99;}

     }

}

Foam::scalar Foam::linearLemaitreDamage::curYieldStress
(
    const scalar plasticMultOld, // old value of the accumulated plastic multiplier
    const scalar DLambda    // plastic multiplier increment
) const
{

    return  stressPlasticStrainSeries_(plasticMultOld+DLambda);

}

Foam::scalar Foam::linearLemaitreDamage::yieldFunction
(
    const scalar qTrial, // Deviatoric trial stress magnitude
    const scalar DLambda, // plastic multiplier 
    const scalar plasticMultOld, // old value of the accumulated plastic multiplier
    scalar damage, // damage value
    const scalar muBar // Scaled shear modulus
) const
{
    // Evaluate current yield function

    return 
        qTrial - 3*muBar*DLambda/(1-damage)
        - curYieldStress(plasticMultOld, DLambda);
}



void Foam::linearLemaitreDamage::newtonLoop
(
    scalar& DLambda, // Plastic multiplier
    scalar& curSigmaY, // Current yield stress
    const scalar plasticMultOld, // old value of the accumulated plastic multiplier
    const scalar qTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar maxMagDEpsilon,    // Max strain increment magnitude
const scalar damage
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    scalar fTrial=0;

    int i = 0;


    fTrial = yieldFunction(qTrial, DLambda, plasticMultOld, damage, muBar);

    scalar residual = 1.0;

    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluations are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =   
            yieldFunction
            (
                 qTrial,  DLambda + finiteDiff_, plasticMultOld,damage,  muBar
            );

        // Numerical derivative of fTrial
        scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;
   
        // Update DLambda
        residual = fTrial/fTrialDerivative;

        DLambda =(1-relaxationFactor_)*DLambda+relaxationFactor_*(DLambda-residual);

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial =   yieldFunction(qTrial, DLambda, plasticMultOld, damage, muBar);


        if (i == MaxNewtonIter_)
        {
            WarningIn("linearLemaitreDamage::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    curSigmaY=curYieldStress(plasticMultOld, DLambda);

}


void Foam::linearLemaitreDamage::calculateHydrostaticStress
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
                "void Foam::linearLemaitreDamage::"
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


void Foam::linearLemaitreDamage::calculateHydrostaticStress
(
    surfaceScalarField& sigmaHyd,
    const surfaceScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        FatalErrorIn
        (
            "void Foam::linearLemaitreDamage::calculateHydrostaticStress\n"
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

// Construct from dictionarylinearElasticMisesPlas
Foam::linearLemaitreDamage::linearLemaitreDamage
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    charLength_("zero", dimLength, 0.0),
    rho_(dict.lookup("rho")),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    stressPlasticStrainSeries_(dict),
    b_("zero", dimless, 0.0),
    s0_("zero", dimPressure, 0.0),
    epsilonD_("zero", dimless, 0.0),
    damageC_
    (
        dict.lookupOrDefault<scalar>("damageC", 0.99)
    ),
    solvePressureEqn_(dict.lookup("solvePressureEqn")),
    debug_(dict.lookup("debug")),
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
    plasticN_.oldTime();
    damage_.oldTime();

    // read paramters for damage model
    charLength_ = dimensionedScalar(dict.lookup("charLength"));
    s0_ = dimensionedScalar(dict.lookup("s0"));
    b_ = dimensionedScalar(dict.lookup("b"));
    epsilonD_ = dimensionedScalar(dict.lookup("epsilonD"));

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
            "linearLemaitreDamage::linearLemaitreDamage::()"
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

Foam::linearLemaitreDamage::~linearLemaitreDamage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearLemaitreDamage::rho() const
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
Foam::linearLemaitreDamage::impK() const
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
Foam::linearLemaitreDamage::impKdiagTensor() const
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


void Foam::linearLemaitreDamage::correct(volSymmTensorField& sigma)
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
                "void Foam::linearLemaitreDamage::"
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

    // Calculate deviatoric strain
    volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // trial equivalent stress
    volScalarField qTrial = sqrt(3.0/2.0)*mag(sTrial);

    // Calculate the yield function
    const volScalarField fTrial =qTrial - sigmaY_;

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilon_.internalField())), SMALL);

    // Store previous iteration for residual calculation
    DEpsilonP_.storePrevIter();


    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    scalarField& qTrialI = qTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    const scalarField& sigmaYI = sigmaY_.internalField();
    scalarField& damageOldNonLocalI = damageNonLocal_.oldTime().internalField();
    scalarField& damageNonLocalI = damageNonLocal_.internalField();

    scalarField& plasticMultOldI = plasticMult_.oldTime().internalField();

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
            plasticNI[cellI] = symmTensor::zero;
            damageNonLocalI[cellI] = damageOldNonLocalI[cellI];
                           
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
                    qTrialI[cellI],
                    mu_.value(),
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
        const scalarField& qTrialP = qTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
        scalarField& damageOldNonLocalP = damageNonLocal_.oldTime().boundaryField()[patchI];
        scalarField& damageNonLocalP = damageNonLocal_.boundaryField()[patchI];

        scalarField& plasticMultOldP = plasticMult_.oldTime().boundaryField()[patchI]; 

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
                plasticNP[faceI] = symmTensor::zero;
                damageNonLocalP[faceI] = damageOldNonLocalP[faceI];
               
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
                        qTrialP[faceI],
                        mu_.value(),
                        maxMagBE,
                        damageNonLocalP[faceI]
                    );

                    // Update increment of yield stress
                    DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];

            }
       }
   }

    plasticMult_ = plasticMult_.oldTime() + DLambda_;

    // Update DEpsilonPEq
    DEpsilonPEq_ = DLambda_/(1-damageNonLocal_);

    // Update DEpsilonP
    DEpsilonP_ = sqrt(3.0/2.0)*DLambda_*plasticN_;

    // Update total plastic strain
    epsilonP_ = epsilonP_.oldTime() + DEpsilonP_;

    // Update equivalent total plastic strain
    epsilonPEq_ = epsilonPEq_.oldTime() + DEpsilonPEq_;

    e=epsilon_-epsilonP_;

    // Calculate deviatoric stress
    const volSymmTensorField s = 2*mu_*dev(e)*(1-damageNonLocal_);

    const volScalarField trEpsilon = tr(epsilon_)*(1-damageNonLocal_);
    calculateHydrostaticStress(sigmaHyd_, trEpsilon);

    // Update the stress
    sigma = sigmaHyd_*I + s;

//****Damage model calculations***//

    volScalarField sigmaEq =
            sqrt((3.0/2.0)*magSqr(dev(sigma/(1-damageNonLocal_))));
    const volScalarField  triaxiality = sigmaHyd_/sigmaEq;
    const volScalarField effSigmaHyd_= sigmaHyd_/(1-damageNonLocal_);

    const scalarField& sigmaEqI = sigmaEq.internalField();
    const scalarField& effSigmaHydI = effSigmaHyd_.internalField();
    const scalarField& epsilonPEqI = epsilonPEq_.internalField();     
    const scalarField& triaxialityI = triaxiality.internalField();
    const scalarField& DEpsilonPeqI=DEpsilonPEq_.internalField();
    scalarField& damageLocalI = damage_.internalField();
    scalarField& damageOldLocalI = damage_.oldTime().internalField();
     
    forAll(damageLocalI, cellI)
    {

        damageFunction
        (
            triaxialityI[cellI],
            epsilonPEqI[cellI],
            sigmaEqI[cellI],
            damageOldLocalI[cellI],
            damageLocalI[cellI],
            DEpsilonPeqI[cellI],
            DLambdaI[cellI],
            effSigmaHydI[cellI]
        );
 
    }
   
    forAll(damage_.boundaryField(), patchI)
    {
        const scalarField& sigmaEqP = sigmaEq.boundaryField()[patchI];
        const scalarField& effSigmaHydP = effSigmaHyd_.boundaryField()[patchI];
        const scalarField& epsilonPEqP = epsilonPEq_.boundaryField()[patchI];     
        const scalarField& triaxialityP = triaxiality.boundaryField()[patchI];
        const scalarField& DLambdaP = DLambda_.boundaryField()[patchI];

        const scalarField& DEpsilonPeqP=DEpsilonPEq_.boundaryField()[patchI];
        const scalarField plasticMultP=plasticMult_.boundaryField()[patchI];
        scalarField& damageLocalP = damage_.boundaryField()[patchI];
        const scalarField damageOldLocalP=damage_.oldTime().boundaryField()[patchI];

        forAll(damageLocalP, faceI)
        {

            damageFunction
            (
                triaxialityP[faceI],
                epsilonPEqP[faceI],
                sigmaEqP[faceI],
                damageOldLocalP[faceI],
                damageLocalP[faceI],
                DEpsilonPeqP[faceI],
                DLambdaP[faceI],
                effSigmaHydP[faceI]
            ); 

        }
    }

    // calculate non-local damage
    fvScalarMatrix damageEqn
    {
        fvm::Sp(1.0, damageNonLocal_)
      - fvm::laplacian(pow(charLength_, 2.0), damageNonLocal_)
     == damage_
    };

    damageEqn.solve();
    gradDamageNonLocal_ = fvc::grad(damageNonLocal_);

//damageNonLocal_=damage_;

}


void Foam::linearLemaitreDamage::correct(surfaceSymmTensorField& sigma)
{
    notImplemented("wip");

}


Foam::scalar Foam::linearLemaitreDamage::residual()
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


void Foam::linearLemaitreDamage::updateTotalFields()
{

    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;

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


Foam::scalar Foam::linearLemaitreDamage::newDeltaT()
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
                "Foam::scalar Foam::linearLemaitreDamage::newDeltaT()"
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
