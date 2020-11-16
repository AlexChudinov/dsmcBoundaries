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

Class
    Foam::DSMCCloud

Description
    Templated base class for dsmc cloud

SourceFiles
    DSMCCloudI.H
    DSMCCloud.C

\*---------------------------------------------------------------------------*/

#ifndef DSMCCloud_H
#define DSMCCloud_H

#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "barycentric.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class BinaryCollisionModel;

template<class CloudType>
class WallInteractionModel;

template<class CloudType>
class InflowBoundaryModel;

/*---------------------------------------------------------------------------*\
                      Class DSMCCloudName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(DSMCCloud);


/*---------------------------------------------------------------------------*\
                         Class DSMCCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DSMCCloud
:
    public Cloud<ParcelType>,
    public DSMCCloudName
{
    // Private Data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudName_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- A list of unique instances of molecule types in the
        //  simulation.  The position of an entry in the list maps to
        //  the label identifying the typeId, i.e. where typeIdList_ =
        //  (N2 O2 CO2) N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Number of real atoms/molecules represented by a parcel
        scalar nParticle_;

        //- A data structure holding which particles are in which cell
        List<DynamicList<ParcelType*>> cellOccupancy_;

        //- A field holding the value of (sigmaT * cR)max for each
        //  cell (see Bird p220). Initialised with the parcels,
        //  updated as required, and read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        volScalarField::Internal collisionSelectionRemainder_;

        //- Parcel constant properties - one for each type
        List<typename ParcelType::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        //gas number density
        volScalarField rhoN_;

        //temperature
        volScalarField T_;

        //velocity
        volVectorField U_;

        //square of velocity
        volScalarField U2_;

        //number of times of cloud was evolved
        uint64_t stepCounter_;

        // References to the cloud sub-models

            //- Binary collision model
            autoPtr<BinaryCollisionModel<DSMCCloud<ParcelType>>>
                binaryCollisionModel_;

            //- Wall interaction model
            autoPtr<WallInteractionModel<DSMCCloud<ParcelType>>>
                wallInteractionModel_;

            //- Inflow boundary model
            autoPtr<InflowBoundaryModel<DSMCCloud<ParcelType>>>
                inflowBoundaryModel_;


    // Private Member Functions

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Record which particles are in which cell
        void buildCellOccupancy();

        //- Initialise the system
        void initialise();

        //- Calculate collisions between molecules
        void collisions();

        //- Reset the data accumulation field values to zero
        void resetFields();

        //- Calculate the volume field data
        void calculateFields();


public:

    // Constructors

        //- Construct given name, mesh and initialisation dictionary.
        DSMCCloud
        (
            const word& cloudName,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        DSMCCloud(const DSMCCloud&) = delete;


    //- Destructor
    virtual ~DSMCCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            // References to the mesh and databases

                //- Return the cloud type
                inline const word& cloudName() const;

                //- Return references to the mesh
                inline const fvMesh& mesh() const;


            // References to the dsmc specific data

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return the idList
                inline const List<word>& typeIdList() const;

                //- Return the number of real particles represented by one
                //  parcel
                inline scalar nParticle() const;

                //- Return the cell occupancy addressing
                inline const List<DynamicList<ParcelType*>>&
                    cellOccupancy() const;

                inline List<DynamicList<ParcelType*>>&
                    cellOccupancy();

                //- Return the sigmaTcRMax field.  non-const access to allow
                // updating.
                inline volScalarField& sigmaTcRMax();

                //- Return the collision selection remainder field.  non-const
                // access to allow updating.
                inline scalarField& collisionSelectionRemainder();

                //- Return all of the constant properties
                inline const List<typename ParcelType::constantProperties>&
                    constProps() const;

                //- Return the constant properties of the given typeId
                inline const typename ParcelType::constantProperties&
                    constProps(label typeId) const;

                //- Return references to the random object
                inline Random& rndGen();

                //- Return macroscopic temperature
                inline const volScalarField::Boundary& boundaryT() const;

                //- Return macroscopic velocity
                inline const volVectorField::Boundary& boundaryU() const;

                //- Return gas number density
                inline const volScalarField::Boundary& boundaryRhoN() const;

            // Kinetic theory helper functions

                //- Generate a random velocity sampled from the Maxwellian speed
                // distribution
                vector equipartitionLinearVelocity
                (
                    scalar temperature,
                    scalar mass
                );

                //- Generate a random internal energy, sampled from the
                // equilibrium distribution (Bird eqn 11.22 and 11.23 and
                // adapting code from DSMC3.FOR)
                scalar equipartitionInternalEnergy
                (
                    scalar temperature,
                    direction internalDegreesOfFreedom
                );


                // From the Maxwellian distribution:
                //- Average particle speed
                inline scalar maxwellianAverageSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianAverageSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- RMS particle speed
                inline scalar maxwellianRMSSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianRMSSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- Most probable speed
                inline scalar maxwellianMostProbableSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianMostProbableSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;


            // Sub-models

                //- Return reference to binary elastic collision model
                inline const BinaryCollisionModel<DSMCCloud<ParcelType>>&
                    binaryCollision() const;

                //- Return non-const reference to binary elastic collision model
                inline BinaryCollisionModel<DSMCCloud<ParcelType>>&
                    binaryCollision();

                //- Return reference to wall interaction model
                inline const WallInteractionModel<DSMCCloud<ParcelType>>&
                    wallInteraction() const;

                //- Return non-const reference to wall interaction model
                inline WallInteractionModel<DSMCCloud<ParcelType>>&
                    wallInteraction();

                //- Return reference to wall interaction model
                inline const InflowBoundaryModel<DSMCCloud<ParcelType>>&
                    inflowBoundary() const;

                //- Return non-const reference to wall interaction model
                inline InflowBoundaryModel<DSMCCloud<ParcelType>>&
                    inflowBoundary();


        // Check

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total internal energy in the system
            inline scalar internalEnergyOfSystem() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;




        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const label celli,
                const vector& U,
                const scalar Ei,
                const label typeId
            );

            //- Evolve the cloud (move, collide)
            void evolve();

            //- Clear the Cloud
            inline void clear();


        // Mapping

            //- Remap the particles to the correct cells following mesh change
            virtual void autoMap(const mapPolyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DSMCCloud&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::DSMCCloud<ParcelType>::constProps(label typeId) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorInFunction
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}

template<class ParcelType>
inline const Foam::word& Foam::DSMCCloud<ParcelType>::cloudName() const
{
    return cloudName_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::DSMCCloud<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::DSMCCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::List<Foam::word>&
Foam::DSMCCloud<ParcelType>::typeIdList() const
{
    return typeIdList_;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloud<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline const Foam::List<Foam::DynamicList<ParcelType*>>&
Foam::DSMCCloud<ParcelType>::cellOccupancy() const
{
    return cellOccupancy_;
}

template<class ParcelType>
Foam::List<Foam::DynamicList<ParcelType *> > &
Foam::DSMCCloud<ParcelType>::cellOccupancy()
{
    return cellOccupancy_;
}


template<class ParcelType>
inline Foam::volScalarField& Foam::DSMCCloud<ParcelType>::sigmaTcRMax()
{
    return sigmaTcRMax_;
}


template<class ParcelType>
inline Foam::scalarField&
Foam::DSMCCloud<ParcelType>::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


template<class ParcelType>
inline const Foam::List<typename ParcelType::constantProperties>&
Foam::DSMCCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline Foam::Random& Foam::DSMCCloud<ParcelType>::rndGen()
{
    return rndGen_;
}

template<class ParcelType>
const Foam::volScalarField::Boundary &Foam::DSMCCloud<ParcelType>::boundaryT() const
{
    return T_.boundaryField();
}

template<class ParcelType>
const Foam::volVectorField::Boundary &Foam::DSMCCloud<ParcelType>::boundaryU() const
{
    return U_.boundaryField();
}

template<class ParcelType>
const Foam::volScalarField::Boundary &Foam::DSMCCloud<ParcelType>::boundaryRhoN() const
{
    return rhoN_.boundaryField();
}

template<class ParcelType>
inline const Foam::BinaryCollisionModel<Foam::DSMCCloud<ParcelType>>&
Foam::DSMCCloud<ParcelType>::binaryCollision() const
{
    return binaryCollisionModel_;
}


template<class ParcelType>
inline Foam::BinaryCollisionModel<Foam::DSMCCloud<ParcelType>>&
Foam::DSMCCloud<ParcelType>::binaryCollision()
{
    return binaryCollisionModel_();
}


template<class ParcelType>
inline const Foam::WallInteractionModel<Foam::DSMCCloud<ParcelType>>&
Foam::DSMCCloud<ParcelType>::wallInteraction() const
{
    return wallInteractionModel_;
}


template<class ParcelType>
inline Foam::WallInteractionModel<Foam::DSMCCloud<ParcelType>>&
Foam::DSMCCloud<ParcelType>::wallInteraction()
{
    return wallInteractionModel_();
}


template<class ParcelType>
inline const Foam::InflowBoundaryModel<Foam::DSMCCloud<ParcelType>>&
Foam::DSMCCloud<ParcelType>::inflowBoundary() const
{
    return inflowBoundaryModel_;
}


template<class ParcelType>
inline Foam::InflowBoundaryModel<Foam::DSMCCloud<ParcelType>>&
Foam::DSMCCloud<ParcelType>::inflowBoundary()
{
    return inflowBoundaryModel_();
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloud<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(typename DSMCCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysMass += cP.mass();
    }

    return nParticle_*sysMass;
}


template<class ParcelType>
inline Foam::vector Foam::DSMCCloud<ParcelType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(typename DSMCCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum += cP.mass()*p.U();
    }

    return nParticle_*linearMomentum;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCCloud<ParcelType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename DSMCCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U());
    }

    return nParticle_*linearKineticEnergy;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCCloud<ParcelType>::internalEnergyOfSystem() const
{
    scalar internalEnergy = 0.0;

    forAllConstIter(typename DSMCCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        internalEnergy += p.Ei();
    }

    return nParticle_*internalEnergy;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloud<ParcelType>::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


template<class ParcelType>
inline Foam::scalarField Foam::DSMCCloud<ParcelType>::maxwellianAverageSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
    return tfld();
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloud<ParcelType>::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField Foam::DSMCCloud<ParcelType>::maxwellianRMSSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(3.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField
Foam::DSMCCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(2.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}

template<class ParcelType>
inline void Foam::DSMCCloud<ParcelType>::clear()
{
    return IDLList<ParcelType>::clear();
}

// ************************************************************************* //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DSMCCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
