#ifndef INPUTSTREAM_H
#define INPUTSTREAM_H

#include <vector>
#include <algorithm>

#include "InflowBoundaryModel.H"
#include "polyMesh.H"
#include "constants.H"
#include "tetIndices.H"

namespace Foam {

template<class CloudType>
class InputStream
        :
        public InflowBoundaryModel<CloudType>
{

    const scalar sqrtPi = sqrt(pi);

    //Name of dict with molecular abundances
    const word componentAbundancesDictName_ = "componentAbundances";

    const word moleculePropertiesDicrName_ = "moleculeProperties";

    const word massDictKey_ = "mass";

    // Private Data

    //- The indices of patches to introduce molecules across
    labelList patches_;

    //Names of the molecules
    wordList componentNames_;

    scalarList componentAbundances_;

    scalarList componentMasses_;

    labelList componentTypeId_;

    void checkTemperature() const;

    vector genMoleculeVelocity(
            Random& rndGen,
            const vector& U,
            double mass,
            double T);

public:

    //- Runtime type information
    TypeName("InputStream");


    // Constructors

    //- Construct from dictionary
    InputStream(
            const dictionary& dict,
            CloudType& cloud);


    //- Destructor
    virtual ~InputStream();


    // Member Functions

    // Mapping

    //- Remap the particles to the correct cells following mesh change
    virtual void autoMap(const mapPolyMesh&);

    //- Introduce particles
    virtual void inflow();
};

template<class CloudType>
void InputStream<CloudType>::checkTemperature() const
{
    const volScalarField::Boundary& boundaryT = this->owner().boundaryT();
    forAll(patches_, p){
        if(min(boundaryT[p]) < small){
            FatalErrorInFunction
                    << "Zero boundary temperature detected, check boundaryT "
                    << "condition." << nl
                    << nl << abort(FatalError);
        }
    }
}

template<class CloudType>
vector InputStream<CloudType>::genMoleculeVelocity(
        Random &rndGen,
        const vector &U,
        double mass,
        double T)
{
    //n - normal pointing inside
    CloudType& cloud = this->owner();
    scalar s = cloud.maxwellianMostProbableSpeed(T, mass) / sqrt(2.0);
    return s * vector(rndGen.scalarNormal(), rndGen.scalarNormal(), rndGen.scalarNormal())
            + U;
}

template<class CloudType>
InputStream<CloudType>::InputStream(const dictionary &dict, CloudType &cloud)
    :
      InflowBoundaryModel<CloudType>(dict, cloud, typeName)
{
    DynamicList<label> patches;
    forAll(cloud.mesh().boundaryMesh(), p){
        const polyPatch& patch = cloud.mesh().boundaryMesh()[p];
        if(isType<polyPatch>(patch)){
            patches.append(p);
        }
    }
    patches_.transfer(patches);

    const dictionary& componentAbundances = this->coeffDict().subDict(
                componentAbundancesDictName_);

    componentNames_ = componentAbundances.toc();
    componentAbundances_.setSize(componentNames_.size());
    componentMasses_.setSize(componentNames_.size());
    componentTypeId_.setSize(componentNames_.size());
    const dictionary& moleculeProperties = cloud.particleProperties().subDict(
                moleculePropertiesDicrName_);
    Info << "Create custom properties of molecules for molecular inflow: " << nl;
    forAll(componentNames_, n)
    {
        componentAbundances_[n] = componentAbundances.lookup<scalar>(componentNames_[n]);
        componentMasses_[n] = moleculeProperties.subDict(componentNames_[n])
                .lookup<scalar>(massDictKey_);
        componentTypeId_[n] = findIndex(cloud.typeIdList(), componentNames_[n]);
        Info << componentNames_[n] << " : "
             << componentTypeId_[n] << " : "
             << componentAbundances_[n] << " : "
             << componentMasses_[n] << nl;
    }

    checkTemperature();
}

template<class CloudType>
InputStream<CloudType>::~InputStream()
{

}

template<class CloudType>
void InputStream<CloudType>::inflow()
{
    CloudType& cloud = this->owner();

    List<DynamicList<typename CloudeType::particleType*>>& occupancy = cloud.cellOccupancy();

    const polyMesh & mesh = cloud.mesh();

    const scalar deltaT = mesh.time().deltaTValue();

    Random & rndGen = cloud.rndGen();

    label particlesInserted = 0;

    const volScalarField::Boundary& boundaryT = cloud.boundaryT();

    const volVectorField::Boundary& boundaryU = cloud.boundaryU();

    const volScalarField::Boundary& rhoN = cloud.boundaryRhoN();

    forAll(componentMasses_, componenti)
    {
        forAll(patches_, p)
        {
            const label patchi = patches_[p];

            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            forAll(patch, facei){
                //const face& f = patch[facei];

                label globalFaceIndex = facei + patch.start();

                label celli = mesh.faceOwner()[globalFaceIndex];

                forAll(occupancy[celli], parceli){
                    cloud.deleteParticle(*occupancy[celli][parceli]);
                }

                occupancy[celli].clear();

                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                        (
                            mesh,
                            celli
                            );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];
                    tetPointRef tet = cellTetIs.tet(mesh_);
                    scalar tetVolume = tet.mag();

                    scalar nP = tetVolume * rhoN[facei] / cloud.nParticle();

                    label numberOfParticles = rndGen.scalar01() > nP - label(nP) ?
                                label(nP) : label(nP) + 1;

                    for(label i = 0; i < numberOfParticles; ++i){

                        const point r0 = tet.randomPoint(rndGen);

                        const vector molecularVecolity = genMoleculeVelocity(
                                    rndGen,
                                    boundaryU[patchi][facei],
                                    componentMasses_[componenti],
                                    boundaryT[patchi][facei]);

                        scalar Ei = cloud.equipartitionInternalEnergy(
                                    boundaryT[patchi][facei],
                                    cloud.constProps(componentTypeId_[componenti]).internalDegreesOfFreedom());

                        cloud.addNewParcel(r0, celli, molecularVecolity, Ei, componentTypeId_[componenti]);

                        particlesInserted++;
                    }
                }
            }
        }
    }

    reduce(particlesInserted, sumOp<label>());

    Info<< "    Particles inserted              = "
        << particlesInserted << nl;
}

template<class CloudType>
void InputStream<CloudType>::autoMap(const mapPolyMesh&) {
    Info << "Just dummy yet" << nl;
}

}

#endif // INPUTSTREAM_H
