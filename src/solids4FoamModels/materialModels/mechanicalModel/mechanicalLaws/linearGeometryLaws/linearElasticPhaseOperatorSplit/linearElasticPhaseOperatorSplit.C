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

#include "linearElasticPhaseOperatorSplit.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticPhaseOperatorSplit, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElasticPhaseOperatorSplit, linGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::linearElasticPhaseOperatorSplit::makeSigma0f() const
{
    if (sigma0fPtr_)
    {
        FatalErrorIn("void Foam::linearElasticPhaseOperatorSplit::makeSigma0f() const")
            << "pointer already set" << abort(FatalError);
    }

    sigma0fPtr_ =
        new surfaceSymmTensorField
        (
            "sigma0f",
            fvc::interpolate(sigma0_)
        );
}


const Foam::surfaceSymmTensorField& Foam::linearElasticPhaseOperatorSplit::sigma0f() const
{
    if (!sigma0fPtr_)
    {
        makeSigma0f();
    }

    return *sigma0fPtr_;
}


void Foam::linearElasticPhaseOperatorSplit::decompose
(
  const symmTensor& epsilon,
  symmTensor& epsilonPositive,
  symmTensor& epsilonNegative
)
{

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        epsilon.xx(), epsilon.xy(), epsilon.xz(),
        epsilon.xy(), epsilon.yy(), epsilon.yz(),
        epsilon.xz(), epsilon.yz(), epsilon.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    //Initialise Vectors which will hold positive and negative Eigenvalues respectively
    Eigen::Vector3d EPosVals=EVals;
    Eigen::Vector3d ENegVals=EVals;

    //Set positive and negative Eignvalue vectors
    for (int i=0;i<3;i++)
    {
        if (EVals(i)>0)
        {
            ENegVals(i)=0.0;
        }
        else
        {
            EPosVals(i)=0.0;
        }
    }
    
    Eigen::Matrix3d P = EPosVals.asDiagonal();
    Eigen::Matrix3d N = ENegVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();
    Eigen::Matrix3d resultPos = EVecs*P*EVecs.inverse();
    Eigen::Matrix3d resultNeg = EVecs*N*EVecs.inverse();

    // set positive elastic strain tensor
    epsilonPositive.xx()=resultPos(0,0);
    epsilonPositive.xy()=resultPos(0,1);
    epsilonPositive.xz()=resultPos(0,2);
    epsilonPositive.yy()=resultPos(1,1);
    epsilonPositive.yz()=resultPos(1,2);
    epsilonPositive.zz()=resultPos(2,2);

    // set negative elastic strain tensor
    epsilonNegative.xx()=resultNeg(0,0);
    epsilonNegative.xy()=resultNeg(0,1);
    epsilonNegative.xz()=resultNeg(0,2);
    epsilonNegative.yy()=resultNeg(1,1);
    epsilonNegative.yz()=resultNeg(1,2);
    epsilonNegative.zz()=resultNeg(2,2);

}


Foam::scalar Foam::linearElasticPhaseOperatorSplit::calcStrainEnergy
(
  const symmTensor& epsilon
)
{

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        epsilon.xx(), epsilon.xy(), epsilon.xz(),
        epsilon.xy(), epsilon.yy(), epsilon.yz(),
        epsilon.xz(), epsilon.yz(), epsilon.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    //Initialise Vectors which will hold positive Eigenvalues 
    Eigen::Vector3d EPosVals=EVals;

    //Set positive Eignvalue vectors

    for (int i=0;i<3;i++)
    {
        if (EVals(i)<0)
        {
            EPosVals(i)=0.0;
        }
    } 

    // Initialse positive strain tensor
    symmTensor epsilonPositive;

    // Eigen::Matrix3d D = EVals.asDiagonal();
    Eigen::Matrix3d P = EPosVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();
    Eigen::Matrix3d resultPos = EVecs*P*EVecs.inverse();

    // Set positive stres tensor
    epsilonPositive.xx()=resultPos(0,0);
    epsilonPositive.xy()=resultPos(0,1);
    epsilonPositive.xz()=resultPos(0,2);
    epsilonPositive.yy()=resultPos(1,1);
    epsilonPositive.yz()=resultPos(1,2);
    epsilonPositive.zz()=resultPos(2,2);

    scalar trEpsilon=tr(epsilon);

    if (trEpsilon>0.0)
    {
        return
            (lambda_.value()*pow(tr(epsilon),2.0))/2.0+mu_.value()*(epsilonPositive&&epsilonPositive);
    }  
    else
    {
        return mu_.value()*(epsilonPositive&&epsilonPositive);
    }
}

void Foam::linearElasticPhaseOperatorSplit::calculateHydrostaticStress
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
                "void Foam::linearElasticPhaseOperatorSplitMisesPlastic::"
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


void Foam::linearElasticPhaseOperatorSplit::calculateHydrostaticStress
(
    surfaceScalarField& sigmaHyd,
    const surfaceScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        FatalErrorIn
        (
            "void Foam::linearElasticPhaseOperatorSplit::calculateHydrostaticStress\n"
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
Foam::linearElasticPhaseOperatorSplit::linearElasticPhaseOperatorSplit
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    charLength_("zero", dimLength, 0.0),
    viscousD_("zero", dimTime, 0.0),
    Gc_("zero", dimPressure, 0.0),
    residualStiff_("zero", dimless, 0.0),
    rho_(dict.lookup("rho")),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    lambda_("lambda", dimPressure, 0.0),
    monolithic_
    (
        dict.lookupOrDefault<Switch>("monolithic", true)
    ),
    solvePressureEqn_
    (
        dict.lookupOrDefault<Switch>("solvePressureEqn", false)
    ),
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
    epsilonPositive_
    (
        IOobject
        (
            "epsilonPositive",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonNegative_
    (
        IOobject
        (
            "epsilonNegative",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    sigma0_
    (
        IOobject
        (
            "sigma0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dict.lookupOrDefault<dimensionedSymmTensor>
        (
            "sigma0",
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        )
    ),
    sigma0fPtr_(NULL)
{

    // material parameters associated with the phase field variable 
    charLength_ = dimensionedScalar(dict.lookup("charLength"));
    viscousD_ = dimensionedScalar(dict.lookup("viscousD"));
    Gc_ = dimensionedScalar(dict.lookup("Gc"));
    residualStiff_ = dimensionedScalar(dict.lookup("residualStiff"));

    // Force storage of strain old time
    epsilon_.oldTime();
    HStrainEnergy_.oldTime();
    strainEnergy_.oldTime();

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
        if (nu_.value() < 0.5)
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_.value() = GREAT;
        }
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
            "linearElasticPhaseOperatorSplitMisesPlastic::linearElasticPhaseOperatorSplitMisesPlastic::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Set first Lame parameter
    if (nu_.value() < 0.5)
    {
        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    }
    else
    {
        lambda_.value() = GREAT;
    }

    // Check for physical Poisson's ratio
    if (nu_.value() < -1.0 || nu_.value() > 0.5)
    {
        FatalErrorIn
        (
            "Foam::linearElasticPhaseOperatorSplit::linearElasticPhaseOperatorSplit\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unphysical Poisson's ratio: nu should be >= -1.0 and <= 0.5"
            << abort(FatalError);
    }

    // Check for incompressibility or quasi-incompressibility
    if (nu_.value() > 0.49 && !solvePressureEqn_)
    {
        WarningIn(type() + "::" + type())
            << "Poisson's ratio is greater than 0.49: "
            << "consider setting 'solvePressureEqn' to 'yes'!" << endl;
    }

    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingCoeff: " << pressureSmoothingCoeff_
            << endl;
    }

    if (gMax(mag(sigma0_)()) > SMALL)
    {
        Info<< "Reading sigma0 initial/residual stress field" << endl;
    }

    if (gMax(mag(sigma0_)()) > SMALL)
    {
        Info<< "Reading sigma0 initial/residual stress field" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticPhaseOperatorSplit::~linearElasticPhaseOperatorSplit()
{
    deleteDemandDrivenData(sigma0fPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElasticPhaseOperatorSplit::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::linearElasticPhaseOperatorSplit::impK() const
{
    if (nu_.value() == 0.5)
    {
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
                mesh(),
                2.0*mu_
            )
        );
    }
    else
    {
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
                mesh(),
                2.0*mu_ + lambda_
            )
        );
    }

}


const Foam::dimensionedScalar& Foam::linearElasticPhaseOperatorSplit::mu() const
{
    return mu_;
}


const Foam::dimensionedScalar& Foam::linearElasticPhaseOperatorSplit::K() const
{
    return K_;
}


const Foam::dimensionedScalar& Foam::linearElasticPhaseOperatorSplit::E() const
{
    return E_;
}


const Foam::dimensionedScalar& Foam::linearElasticPhaseOperatorSplit::nu() const
{
    return nu_;
}


const Foam::dimensionedScalar& Foam::linearElasticPhaseOperatorSplit::lambda() const
{
    return lambda_;
}

void Foam::linearElasticPhaseOperatorSplit::calcPhase()
{

    if (monolithic_==false)
    {
        Info<<"Updating Phase Field"<<endl;Info<<endl;
    }

    // take references to internal fields
    symmTensorField& epsilonI = epsilon_.internalField();
    scalarField& strainEnergyI = strainEnergy_.internalField();
    scalarField& HStrainEnergyI = HStrainEnergy_.internalField();
    scalarField& HStrainEnergyOldI = HStrainEnergy_.oldTime().internalField();

    // Calculate history strain energy

    forAll(HStrainEnergyI, cellI)
    { 
         strainEnergyI[cellI] = calcStrainEnergy(epsilonI[cellI]);
         HStrainEnergyI[cellI] = max(strainEnergyI[cellI],HStrainEnergyOldI[cellI]);
    }

    forAll(HStrainEnergy_.boundaryField(), patchI)
    {
        // take references to boundary fields
        symmTensorField& epsilonP = epsilon_.boundaryField()[patchI];
        scalarField& strainEnergyP = strainEnergy_.boundaryField()[patchI];
        scalarField& HStrainEnergyP = HStrainEnergy_.boundaryField()[patchI];
        scalarField& HStrainEnergyOldP = HStrainEnergy_.oldTime().boundaryField()[patchI];

        forAll(HStrainEnergyP, faceI)
        {
            strainEnergyP[faceI] = calcStrainEnergy(epsilonP[faceI]);
            HStrainEnergyP[faceI] = max(strainEnergyP[faceI],HStrainEnergyOldP[faceI]);
        }
    }

    // This variable is to ensure that OpenFOAM doesn't give a dimension error  
    dimensionedScalar diml=charLength_/charLength_.value();

    // Calculate the phase field variable d
    fvScalarMatrix DEqn
    (
        -fvm::laplacian(pow(charLength_, 2.0), Dpf_)
        +fvm::Sp(1.0, Dpf_)
       // +fvm::ddt(charLength_.value()*viscousD_/Gc_.value(),Dpf_)
        ==   
        (2.0*charLength_*HStrainEnergy_/(Gc_*diml))
        -fvm::Sp((2.0*charLength_*HStrainEnergy_/(Gc_*diml)), Dpf_)
    );

    DEqn.solve();
}

void Foam::linearElasticPhaseOperatorSplit::correct(volSymmTensorField& sigma)
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
                "void Foam::linearElasticPhaseOperatorSplitMisesPlastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        epsilon_.replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
        );
    }

    // take references to the internal fields
    symmTensorField& epsilonI =epsilon_.internalField();
    symmTensorField& epsilonPositiveI =epsilonPositive_.internalField();
    symmTensorField& epsilonNegativeI =epsilonNegative_.internalField();

    // Decompose elastic strain tensor into it's positive and negative components
    forAll(epsilonI, cellI)
    {
          decompose
          (
           epsilonI[cellI],
           epsilonPositiveI[cellI],
           epsilonNegativeI[cellI]
          );
    }

    forAll(epsilon_.boundaryField(), patchI)
    {
        // take references to the boundary fields
        symmTensorField& epsilonP = epsilon_.boundaryField()[patchI];
        symmTensorField& epsilonPositiveP = epsilonPositive_.boundaryField()[patchI];
        symmTensorField& epsilonNegativeP = epsilonNegative_.boundaryField()[patchI];

        // Decompose elastic strain tensor into it's positive and negative components
        forAll(epsilonP, faceI)
        {
            decompose
            (
                epsilonP[faceI],
                epsilonPositiveP[faceI],
                epsilonNegativeP[faceI]
            );
        }
    }

    // Hooke's law 

    const volScalarField trEpsilonPos = tr(epsilonPositive_);
    const volScalarField trEpsilonNeg = tr(epsilonNegative_);

    // Calculate positive and negative stress tensor's
    volSymmTensorField sigmaPos = 2.0*mu_*epsilonPositive_ + lambda_*trEpsilonPos*I + sigma0_;
    volSymmTensorField sigmaNeg = 2.0*mu_*epsilonNegative_ + lambda_*trEpsilonNeg*I + sigma0_;

    sigma = sigmaPos + sigmaNeg;

    if (monolithic_)
    {
        calcPhase();
    }


}


void Foam::linearElasticPhaseOperatorSplit::correct(surfaceSymmTensorField& sigma)
{
    notImplemented
    (
        "linearElasticPlasticPhaseFracture::correct(surfaceSymmTensorField& sigma)"
    );
}


// ************************************************************************* //
