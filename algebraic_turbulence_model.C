#include "mixingLengthModelFINAL.H"
#include "fvOptions.H"
#include "bound.H"
namespace Foam
{
namespace RASModels
{
// --- Dummy implementation for k() so that the class is concrete ---
template<class BasicTurbulenceModel>
tmp<volScalarField> mixingLengthModelFINAL<BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>(
        new volScalarField
        (
            IOobject
            (
                "k",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("zero", dimensionSet(0,2,-2,0,0,0,0), 0.0)
        )
    );
}
template<class BasicTurbulenceModel>
void mixingLengthModelFINAL<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = TV_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}
template<class BasicTurbulenceModel>
mixingLengthModelFINAL<BasicTurbulenceModel>::mixingLengthModelFINAL
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    Lmix_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Lmix",
            this->coeffDict_,
            dimensionSet(0, 1, 0, 0, 0, 0, 0),
            0.01
        )
    ),   
    TV_
    (
        IOobject
        (
            IOobject::groupName("Turbulent_Viscosity", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )    
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
bool mixingLengthModelFINAL<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        
        Lmix_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}
template<class BasicTurbulenceModel>
void mixingLengthModelFINAL<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();
    tmp<volTensorField> GradientU = fvc::grad(U);
    const volTensorField& TRUEgradU = GradientU();
    const volSymmTensorField S(symm(TRUEgradU));
    const volScalarField SijSij(mag(S));
    this->TV_ = sqr(Lmix_) * sqrt(2.0) * SijSij;
    correctNut();
}
} // End namespace RASModels
} // End namespace Foam
#include "IncompressibleTurbulenceModel.H"
#include "transportModel.H"
#include "addToRunTimeSelectionTable.H"
namespace Foam {
namespace RASModels {
    typedef mixingLengthModelFINAL<IncompressibleTurbulenceModel<transportModel> > incompressiblemixingLengthModelFINALType;
    typedef RASModel<IncompressibleTurbulenceModel<transportModel> >          incompressibleRASModelType;
    addToRunTimeSelectionTable(incompressibleRASModelType, incompressiblemixingLengthModelFINALType, dictionary);
}
}

