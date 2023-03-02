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

#include "linearElasticPlasticPhaseFracture.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticPlasticPhaseFracture, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElasticPlasticPhaseFracture, linGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar linearElasticPlasticPhaseFracture::LoopTol_ = 1e-4;

    // Maximum number of iterations for Newton loop
    label linearElasticPlasticPhaseFracture::MaxNewtonIter_ = 1000;

    // finiteDiff is the delta for finite difference differentiation
    scalar linearElasticPlasticPhaseFracture::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar linearElasticPlasticPhaseFracture::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);


} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearElasticPlasticPhaseFracture::newtonLoop
(
    scalar& DLambda, // plastic mulitplier            
    scalar& curSigmaY, // current yield stress
    const scalar& plasticMult, // accumulated plastic multiplier  
    const scalar& magSTrial, // magnitude of the trial deviatoric stress       
    const scalar& Dpf // phase field variable
) const
{
    // initialise and set yield stress
    scalar sigmaY = stressPlasticStrainSeries_(plasticMult+DLambda);

    //set equivalent stress
    scalar qTrial = sqrt(3.0/2.0)*magSTrial;
    
    // initialise variables for Newton loop
    scalar fTrial, dfTrial;
    scalar DLambdaPrev = 0.0;

    //initialise crack degredation function
    scalar gd = pow(Dpf,2.0); 

    //Newton loop
    for (int iteration = 0; iteration < 200; iteration++)
    {

        sigmaY = stressPlasticStrainSeries_(plasticMult+DLambda);

        fTrial = qTrial - gd*3*mu_.value()*DLambda-gd*sigmaY;

        scalar DSigmaY = (
                           stressPlasticStrainSeries_(plasticMult+DLambda+1e-8)
                          -stressPlasticStrainSeries_(plasticMult+DLambda)
                         )/1e-8;

        dfTrial = -gd*3*mu_.value()-gd*DSigmaY;
        DLambda -= fTrial/dfTrial;


        scalar residual = (DLambdaPrev-DLambda)/(DLambda);

        if (mag(residual) < LoopTol_)
        {
            break;
        }

        DLambdaPrev = DLambda;
     
    }

    curSigmaY=sigmaY;

}


void Foam::linearElasticPlasticPhaseFracture::calculateHydrostaticStress
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
                "void Foam::linearElasticPlasticPhaseFracture::"
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


void Foam::linearElasticPlasticPhaseFracture::calculateHydrostaticStress
(
    surfaceScalarField& sigmaHyd,
    const surfaceScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        FatalErrorIn
        (
            "void Foam::linearElasticPlasticPhaseFracture::calculateHydrostaticStress\n"
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
Foam::linearElasticPlasticPhaseFracture::linearElasticPlasticPhaseFracture
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    charLength_("zero", dimLength, 0.0),
    Gc_("zero", dimPressure, 0.0),
    w0_("zero", dimPressure, 0.0),
    residualStiff_   
    (
        dict.lookupOrDefault<scalar>("residualStiffness", 0.0)
    ),
    stressPlasticStrainSeries_(dict),
    solvePressureEqn_(dict.lookup("solvePressureEqn")),
    debug_(dict.lookup("debug")),
    monolithic_
    (
        dict.lookupOrDefault<Switch>("monolithic", true)
    ),
    pressureSmoothingCoeff_
    (
        dict.lookupOrDefault<scalar>("pressureSmoothingCoeff", 1.0)
    ),
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
        dimensionedScalar("zero", dimless,0.0),
        zeroGradientFvPatchScalarField::typeName
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
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
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
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
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
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
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
    magS_
    (
        IOobject
        (
            "magS",
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
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
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
    plasticMult_
    (
        IOobject
        (
            "plasticMult",
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
    plasticN_.oldTime();
    plasticMult_.oldTime();
    plasticStrainEnergy_.oldTime();

    // parameters for phase-field model
    charLength_ = dimensionedScalar(dict.lookup("charLength"));
    //viscousD_ = dimensionedScalar(dict.lookup("viscousD"));
    Gc_ = dimensionedScalar(dict.lookup("Gc"));
    w0_ = dimensionedScalar(dict.lookup("w0"));
    Dpf_=1.0;  

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
            "linearElasticPlasticPhaseFracture::linearElasticPlasticPhaseFracture::()"
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

Foam::linearElasticPlasticPhaseFracture::~linearElasticPlasticPhaseFracture()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElasticPlasticPhaseFracture::rho() const
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
Foam::linearElasticPlasticPhaseFracture::impK() const
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
Foam::linearElasticPlasticPhaseFracture::impKdiagTensor() const
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

void Foam::linearElasticPlasticPhaseFracture::calcPhase()
{

    if (!monolithic_){Info<<"Updating Phase Field"<<endl;}

    // set fields needed to calculate strain energy contributions
    volSymmTensorField epsilonEl_ = epsilon_ - (epsilonP_.oldTime() + DEpsilonP_);
    volScalarField trEpsilon = tr(epsilonEl_);
    volScalarField sigmaY1 = sigmaY_ + DSigmaY_;

    // set initial values for strain energy contributions
    plasticStrainEnergy_ = plasticStrainEnergy_.oldTime() + DLambda_*(sigmaY1);
    volScalarField effPlasticStrainEnergy = plasticStrainEnergy_;
    strainEnergy_ = ((0.5*K_*trEpsilon*trEpsilon) + (mu_*dev(epsilonEl_)&&dev(epsilonEl_)));

    // references to internal fields
    scalarField& strainEnergyI = strainEnergy_.internalField();
    scalarField& HStrainEnergyI = HStrainEnergy_.internalField();
    scalarField& HStrainEnergyOldI = HStrainEnergy_.oldTime().internalField();
    scalarField& plasticStrainEnergyI = plasticStrainEnergy_.internalField();
    scalarField& effPlasticStrainEnergyI = effPlasticStrainEnergy.internalField();

    // Set elastic history strain energy variable and effective plastic strain energy
    forAll(HStrainEnergyI, cellI)
    { 
         HStrainEnergyI[cellI] = max(strainEnergyI[cellI],HStrainEnergyOldI[cellI]);
         effPlasticStrainEnergyI[cellI] = effPlasticStrainEnergyI[cellI] - w0_.value();

         if (plasticStrainEnergyI[cellI] < w0_.value())
         {
             effPlasticStrainEnergyI[cellI] = 0;
         }
    }

    // Set elastic history strain energy variable and effective plastic strain energy
    forAll(HStrainEnergy_.boundaryField(), patchI)
    {
        // references to boundary fields
        scalarField& strainEnergyP = strainEnergy_.boundaryField()[patchI];
        scalarField& HStrainEnergyP = HStrainEnergy_.boundaryField()[patchI];
        scalarField& HStrainEnergyOldP = HStrainEnergy_.oldTime().boundaryField()[patchI];
        scalarField& plasticStrainEnergyP = plasticStrainEnergy_.boundaryField()[patchI];
        scalarField& effPlasticStrainEnergyP = effPlasticStrainEnergy.boundaryField()[patchI];

        forAll(HStrainEnergyP, faceI)
        {
            HStrainEnergyP[faceI] = max(strainEnergyP[faceI],HStrainEnergyOldP[faceI]);
            effPlasticStrainEnergyP[faceI] = effPlasticStrainEnergyP[faceI] - w0_.value();

            if (plasticStrainEnergyP[faceI] < w0_.value())
            {
                effPlasticStrainEnergyP[faceI] = 0;
            }
        }
    }

    // Set variables which drive the crack growth
    volScalarField elasticCoeff = (4.0*charLength_*HStrainEnergy_)/(Gc_);
    volScalarField plasticCoeff = (4.0*charLength_*effPlasticStrainEnergy)/(Gc_);
    volScalarField coeff = elasticCoeff + plasticCoeff;

    // this variable is used to ensure that OpenFOAM doesn't give a dimension error
    dimensionedScalar diml = charLength_/charLength_.value();

    volScalarField totalStrainEnergy = HStrainEnergy_ + plasticStrainEnergy_;

    fvScalarMatrix DpfEqn
    (
       -fvm::laplacian(4.0*pow(charLength_,2.0), Dpf_)
       +fvm::Sp(1.0, Dpf_)
       +fvm::Sp(coeff/diml, Dpf_)
       ==   
       one_
    );

    DpfEqn.solve();

}

void Foam::linearElasticPlasticPhaseFracture::correct(volSymmTensorField& sigma)
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
                "void Foam::linearElasticPlasticPhaseFracture::"
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
    const volSymmTensorField e = dev(epsilon_);

    volScalarField gd = pow(Dpf_ ,2.0);
  
    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = gd*2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Calculate the yield function
    const volScalarField fTrial = sqrt(3.0/2.0)*mag(sTrial) - gd*sigmaY_;

    // Store previous iteration for residual calculation
    DEpsilonP_.storePrevIter();

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const scalarField& DpfI = Dpf_.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    scalarField& plasticMultOldI = plasticMult_.oldTime().internalField();
    const scalarField& sigmaYI = sigmaY_.internalField();

    // Calculate DLambda and plasticN
    forAll(fTrialI, cellI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[cellI]);
        if (magS > SMALL)
        {
            plasticNI[cellI] = sqrt(3.0/2.0)*sTrialI[cellI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[cellI] < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
            plasticNI[cellI] = symmTensor::zero;
        }
        else 
        { 

            DLambdaI[cellI]=0;

            scalar curSigmaY = 0.0; // updated in loop below

            // Calculates DLambda using Newtons's method
            newtonLoop
            (
                DLambdaI[cellI],
                curSigmaY,
                plasticMultOldI[cellI],
                magS,
                DpfI[cellI]
            );

            // Update increment of yield stress
             DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
            
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& DpfP = Dpf_.boundaryField()[patchI];
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        scalarField& plasticMultOldP = plasticMult_.oldTime().boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Calculate direction plasticN
            const scalar magS = mag(sTrialP[faceI]);
            if (magS > SMALL)
            {
                plasticNP[faceI] = sqrt(3.0/2.0)*sTrialP[faceI]/magS;
            }

            // Calculate DEpsilonPEq
            if (fTrialP[faceI] < SMALL)
            {
                // elasticity
                DSigmaYP[faceI] = 0.0;
                DLambdaP[faceI] = 0.0;
                plasticNP[faceI] = symmTensor::zero;
            }
            else
            {
  
                DLambdaP[faceI]=0;
                scalar curSigmaY = 0.0; // updated in loop below

                // Calculate DLambda and curSigmaY
                newtonLoop
                (
                    DLambdaP[faceI],
                    curSigmaY,
                    plasticMultOldP[faceI],
                    magS,
                    DpfP[faceI]
                );

                // Update increment of yield stress
                DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];

            }
        }
    }

    plasticMult_=plasticMult_.oldTime()+DLambda_;

    // Update DEpsilonPEq
    DEpsilonPEq_ = DLambda_;

    // Update DEpsilonP
    DEpsilonP_ = DLambda_*plasticN_;

    // Update total plastic strain
    epsilonP_ = epsilonP_.oldTime() + DEpsilonP_;

    // Update equivalent total plastic strain
    epsilonPEq_ = epsilonPEq_.oldTime() + DEpsilonPEq_;

    // Calculate deviatoric stress
    const volSymmTensorField s = (sTrial - gd*2*mu_*DEpsilonP_);

    // Calculate the hydrostatic pressure
    const volScalarField trEpsilon = tr(epsilon_);
    calculateHydrostaticStress(sigmaHyd_, trEpsilon);

    // Update the stress
    sigma = gd*sigmaHyd_*I + s;

    magS_=mag(s);//CHECK THIS

    const volScalarField sigmaEq =
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
          + dimensionedScalar("SMALL", dimPressure, SMALL);

     if (monolithic_)
     {
         calcPhase();
     }

}


void Foam::linearElasticPlasticPhaseFracture::correct(surfaceSymmTensorField& sigma)
{
    notImplemented
    (
        "linearElasticPlasticPhaseFracture::correct(surfaceSymmTensorField& sigma)"
    );
}


Foam::scalar Foam::linearElasticPlasticPhaseFracture::residual()
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


void Foam::linearElasticPlasticPhaseFracture::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;
//    sigmaYf_ += DSigmaYf_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
//    epsilonPEq_ += DEpsilonPEq_;
//    epsilonPEqf_ += DEpsilonPEqf_;
//    epsilonP_ += DEpsilonP_;
//    epsilonPf_ += DEpsilonPf_;

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


Foam::scalar Foam::linearElasticPlasticPhaseFracture::newDeltaT()
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
                "Foam::scalar Foam::linearElasticPlasticPhaseFracture::newDeltaT()"
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
