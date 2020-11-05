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

#include "DSMCParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::DSMCParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = cloud.pMesh();

    // For reduced-D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;

    while (td.keepParticle && !td.switchProcessor && p.stepFraction() < 1)
    {
        Utracking = U_;

        // Apply correction to velocity to constrain tracking for
        // reduced-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        const scalar f = 1 - p.stepFraction();
        p.trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::DSMCParcel<ParcelType>::hitPatch(TrackCloudType&, trackingData&)
{
    return false;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DSMCParcel<ParcelType>::hitProcessorPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DSMCParcel<ParcelType>::hitWallPatch
(
    TrackCloudType& cloud,
    trackingData&
)
{
    cloud.wallInteraction().correct(*this);
}


template<class ParcelType>
void Foam::DSMCParcel<ParcelType>::transformProperties
(
    const transformer& transform
)
{
    ParcelType::transformProperties(transform);
    U_ = transform.transform(U_);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "DSMCParcelIO.C"

// ************************************************************************* //
