/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DSMCCloud.H"
#include "BinaryCollisionModel.H"
#include "WallInteractionModel.H"
#include "InflowBoundaryModel.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"
#include "polyMeshTetDecomposition.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] = typename ParcelType::constantProperties(molDict);
    }
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(typename DSMCCloud<ParcelType>, *this, iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::initialise()
{
    Info<< nl << "Initialising particles" << endl;

    const dictionary& moleculeAbundancesDict
    (
        particleProperties_.subDict("initialAbundances")
    );

    List<word> molecules(moleculeAbundancesDict.toc());

    Field<scalar> moleculeAbundancies(molecules.size());

    forAll(molecules, i)
    {
        moleculeAbundancies[i] = moleculeAbundancesDict.lookup<scalar>(molecules[i]);
    }

    label mostAbundantType(findMax(moleculeAbundancies));

    const typename ParcelType::constantProperties& cP = constProps(mostAbundantType);

    forAll(mesh_.cells(), celli)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh_,
            celli
        );

        forAll(cellTets, tetI)
        {
            const tetIndices& cellTetIs = cellTets[tetI];
            tetPointRef tet = cellTetIs.tet(mesh_);
            scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName(molecules[i]);

                label typeId(findIndex(typeIdList_, moleculeName));

                if (typeId == -1)
                {
                    FatalErrorInFunction
                        << "typeId " << moleculeName << "not defined." << nl
                        << abort(FatalError);
                }

                const typename ParcelType::constantProperties& cP = constProps(typeId);

                // Calculate the number of particles required
                scalar particlesRequired = rhoN_[celli] * moleculeAbundancies[i] / nParticle_;

                particlesRequired =
                        (particlesRequired - label(particlesRequired)) > rndGen_.scalar01() ?
                            label(particlesRequired) + 1 : label(particlesRequired);

                for (label pI = 0; pI < particlesRequired; pI++)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U = equipartitionLinearVelocity(T_[celli], cP.mass());

                    scalar Ei = equipartitionInternalEnergy(
                        T_[celli],
                        cP.internalDegreesOfFreedom()
                    );

                    U += U_[celli];

                    addNewParcel(p, celli, U, Ei, typeId);
                }
            }
        }
        // Initialise the sigmaTcRMax_ field to the product of the cross section of
        // the most abundant species and the most probable thermal speed (Bird,
        // p222-223)



        sigmaTcRMax_[celli] = cP.sigmaT()*maxwellianMostProbableSpeed
        (
            T_[celli],
            cP.mass()
        );
    }

    sigmaTcRMax_.correctBoundaryConditions();
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::collisions()
{
    if (!binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List<DynamicList<label>> subCells(8);

    scalar deltaT = mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    forAll(cellOccupancy_, celli)
    {
        const DynamicList<ParcelType*>& cellParcels(cellOccupancy_[celli]);

        label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            const point& cC = mesh_.cellCentres()[celli];

            forAll(cellParcels, i)
            {
                const ParcelType& p = *cellParcels[i];
                vector relPos = p.position() - cC;

                label subCell =
                    pos0(relPos.x()) + 2*pos0(relPos.y()) + 4*pos0(relPos.z());

                subCells[subCell].append(i);
                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = sigmaTcRMax_[celli];

            scalar selectedPairs =
                collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC - 1)*nParticle_*sigmaTcRMax*deltaT
               /mesh_.cellVolumes()[celli];

            label nCandidates(selectedPairs);
            collisionSelectionRemainder_[celli] = selectedPairs - nCandidates;
            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen_.sampleAB<label>(0, nC);

                // Declare the second collision candidate
                label candidateQ = -1;

                List<label> subCellPs = subCells[whichSubCell[candidateP]];
                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        candidateQ = subCellPs[rndGen_.sampleAB<label>(0, nSC)];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        candidateQ = rndGen_.sampleAB<label>(0, nC);
                    } while (candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = rndGen_.sampleAB<label>(0, nC);

                // // Select a possible second collision candidate
                // label candidateQ = rndGen_.sampleAB<label>(0, nC);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = rndGen_.sampleAB<label>(0, nC);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                ParcelType& parcelP = *cellParcels[candidateP];
                ParcelType& parcelQ = *cellParcels[candidateQ];

                scalar sigmaTcR = binaryCollision().sigmaTcR
                (
                    parcelP,
                    parcelQ
                );

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this

                if (sigmaTcR > sigmaTcRMax_[celli])
                {
                    sigmaTcRMax_[celli] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > rndGen_.scalar01())
                {
                    binaryCollision().collide
                    (
                        parcelP,
                        parcelQ
                    );

                    collisions++;
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    sigmaTcRMax_.correctBoundaryConditions();

    if (collisionCandidates)
    {
        Info<< "    Collisions                      = "
            << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl
            << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::calculateFields()
{
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const label celli,
    const vector& U,
    const scalar Ei,
    const label typeId
)
{
    this->addParticle(new ParcelType(mesh_, position, celli, U, Ei, typeId));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DSMCCloud<ParcelType>::DSMCCloud
(
    const word& cloudName,
    const fvMesh& mesh
)
:
    Cloud<ParcelType>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_
    (
        particleProperties_.template lookup<scalar>("nEquivalentParticles")
    ),
    cellOccupancy_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( dimensionSet(0, 3, -1, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    collisionSelectionRemainder_
    (
        IOobject
        (
            this->name() + ":collisionSelectionRemainder",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),
    constProps_(),
    rndGen_(label(971501) + 1526*Pstream::myProcNo()),
    rhoN_
    (
        IOobject
        (
            this->name() + "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( dimensionSet(0, -3, 0, 0, 0), vSmall)
    ),
    T_
    (
        volScalarField
        (
            IOobject
            (
                "T",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    U_
    (
        volVectorField
        (
            IOobject
            (
                "U",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    binaryCollisionModel_(),
    wallInteractionModel_(),
    inflowBoundaryModel_()
{
    clear();
    buildConstProps();
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DSMCCloud<ParcelType>::~DSMCCloud()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::evolve()
{
    typename ParcelType::trackingData td(*this);

    // Reset the data collection fields

    if (debug)
    {
        this->dumpParticlePositions();
    }

    // Insert new particles from the inflow boundary
    this->inflowBoundary().inflow();

    // Move the particles ballistically with their current velocities
    Cloud<ParcelType>::move(*this, td, mesh_.time().deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    // Calculate new velocities via stochastic collisions
    collisions();

    // Calculate the volume field data
    calculateFields();
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::info() const
{
    label nDSMCParticles = this->size();
    reduce(nDSMCParticles, sumOp<label>());

    scalar nMol = nDSMCParticles*nParticle_;

    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar internalEnergy = internalEnergyOfSystem();
    reduce(internalEnergy, sumOp<scalar>());

    Info<< "Cloud name: " << this->name() << nl
        << "    Number of dsmc particles        = "
        << nDSMCParticles
        << endl;

    if (nDSMCParticles)
    {
        Info<< "    Number of molecules             = "
            << nMol << nl
            << "    Mass in system                  = "
            << returnReduce(massInSystem(), sumOp<scalar>()) << nl
            << "    Average linear momentum         = "
            << linearMomentum/nMol << nl
            << "   |Average linear momentum|        = "
            << mag(linearMomentum)/nMol << nl
            << "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average internal energy         = "
            << internalEnergy/nMol << nl
            << "    Average total energy            = "
            << (internalEnergy + linearKineticEnergy)/nMol
            << endl;
    }
}


template<class ParcelType>
Foam::vector Foam::DSMCCloud<ParcelType>::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *rndGen_.sampleNormal<vector>();
}


template<class ParcelType>
Foam::scalar Foam::DSMCCloud<ParcelType>::equipartitionInternalEnergy
(
    scalar temperature,
    direction iDof
)
{
    scalar Ei = 0.0;

    if (iDof == 0)
    {
        return Ei;
    }
    else if (iDof == 2)
    {
        // Special case for iDof = 2, i.e. diatomics;
        Ei = -log(rndGen_.scalar01())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*iDof - 1;
        scalar energyRatio;
        scalar P = -1;

        do
        {
            energyRatio = 10*rndGen_.scalar01();
            P = pow((energyRatio/a), a)*exp(a - energyRatio);
        } while (P < rndGen_.scalar01());

        Ei = energyRatio*physicoChemical::k.value()*temperature;
    }

    return Ei;
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(typename DSMCCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<ParcelType>::autoMap(mapper);

    // Update the cell occupancy field
    cellOccupancy_.setSize(mesh_.nCells());
    buildCellOccupancy();

    // Update the inflow BCs
    this->inflowBoundary().autoMap(mapper);
}


// ************************************************************************* //
