/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
inline const Foam::volScalarField&
Foam::DSMCCloud<ParcelType>::rhoN() const
{
    return rhoN_;
}

template<class ParcelType>
inline void Foam::DSMCCloud<ParcelType>::clear()
{
    return IDLList<ParcelType>::clear();
}


// ************************************************************************* //
