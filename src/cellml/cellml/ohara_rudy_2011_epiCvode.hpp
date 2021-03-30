#ifdef CHASTE_CVODE
#ifndef CELLOHARA_RUDY_2011_EPIFROMCELLMLCVODE_HPP_
#define CELLOHARA_RUDY_2011_EPIFROMCELLMLCVODE_HPP_

//! @file
//!
//! This source file was generated from CellML by chaste_codegen version 0.6.3
//!
//! Model: ohara_rudy_2011_epi
//!
//! Processed by chaste_codegen: https://github.com/ModellingWebLab/chaste-codegen
//!     (translator: chaste_codegen, model type: AnalyticCvode)
//! on 2021-03-29 17:27:20
//!
//! <autogenerated>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractStimulusFunction.hpp"
#include "AbstractCvodeCell.hpp"

class Cellohara_rudy_2011_epiFromCellMLCvode : public AbstractCvodeCell
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCvodeCell >(*this);
        
    }

    //
    // Settable parameters and readable variables
    //

public:

    boost::shared_ptr<RegularStimulus> UseCellMLDefaultStimulus();
    double GetIntracellularCalciumConcentration();
    Cellohara_rudy_2011_epiFromCellMLCvode(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~Cellohara_rudy_2011_epiFromCellMLCvode();
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);
    void EvaluateYDerivatives(double var_chaste_interface__environment__time, const N_Vector rY, N_Vector rDY);
    N_Vector ComputeDerivedQuantities(double var_chaste_interface__environment__time, const N_Vector & rY);
    void EvaluateAnalyticJacobian(double var_chaste_interface__environment__time, N_Vector rY, N_Vector rDY, CHASTE_CVODE_DENSE_MATRIX rJacobian, N_Vector rTmp1, N_Vector rTmp2, N_Vector rTmp3);
};

// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Cellohara_rudy_2011_epiFromCellMLCvode)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const Cellohara_rudy_2011_epiFromCellMLCvode * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, Cellohara_rudy_2011_epiFromCellMLCvode * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)Cellohara_rudy_2011_epiFromCellMLCvode(p_solver, p_stimulus);
        }

    }

}

#endif // CELLOHARA_RUDY_2011_EPIFROMCELLMLCVODE_HPP_
#endif // CHASTE_CVODE