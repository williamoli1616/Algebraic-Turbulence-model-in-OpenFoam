#ifndef mixingLengthModelFINAL_H
#define mixingLengthModelFINAL_H
#include "RASModel.H"
#include "eddyViscosity.H"
namespace Foam
{
namespace RASModels
{
template<class BasicTurbulenceModel>
class mixingLengthModelFINAL
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
        mixingLengthModelFINAL(const mixingLengthModelFINAL&) = delete;
        void operator=(const mixingLengthModelFINAL&) = delete;
protected:
            dimensionedScalar Lmix_;
            volScalarField TV_;
        virtual void correctNut();
public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;
    TypeName("mixingLengthModelFINAL");
    // Constructors
        mixingLengthModelFINAL
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );
    //- Destructor
    virtual ~mixingLengthModelFINAL() = default;
        virtual bool read();
        virtual tmp<volScalarField> TV() const
        {
            return TV_;
        }
        virtual tmp<volScalarField> k() const;
        virtual void correct();
};
} // End namespace RASModels
} // End namespace Foam
#endif
