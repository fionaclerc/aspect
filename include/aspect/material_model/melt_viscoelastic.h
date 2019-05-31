/*
  Copyright (C) 2014 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_melt_viscoelastic_h
#define _aspect_material_model_melt_viscoelastic_h

#include <deal.II/base/function_lib.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>
//#include <aspect/material_model/viscoelastic.h>
#include <aspect/melt.h>
namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that applies "effective viscoelastic" viscosities to a ''base model''
     * chosen from any of the other available material models.
     * All other properties are derived from the base model.
     * @ingroup MaterialModels
     */

    /**
     * Additional output fields for the elastic shear modulus to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class ElasticAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        ElasticAdditionalOutputs(const unsigned int n_points);

        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * Elastic shear moduli at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> elastic_shear_moduli;
    };



    template <int dim>
    class MeltViscoelastic : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Initialize the base model at the beginning of the run.
         */
        virtual
        void initialize();

        /**
         * Update the base model and viscosity function at the beginning of
         * each timestep.
         */
       // virtual
       // void update();

	/**
	 * Function to compute viscoelasticity.
	 */
	virtual
	void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
//        virtual
//        void
//        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
//                  typename Interface<dim>::MaterialModelOutputs &out) const;
//		  const MaterialModel::MaterialModelInputs<dim> &inVE,
//                  MaterialModel::MaterialModelOutputs<dim> &outVE) const;
        /**
         * Method to declare parameters related to viscoelastic melt model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to viscoelastic melt model
         */
        virtual void
        parse_parameters (ParameterHandler &prm);

        /**
         * Method that indicates whether material is compressible. Viscoelastic melt
         * model is compressible if and only if base model is compressible.
	 * Make sure this returns false.
         */
        virtual bool is_compressible () const;

        /** 
         * Method that calculates reference viscosity (from melt model).
         */
        virtual double reference_viscosity () const;
	virtual double reference_darcy_coefficient () const;
	
	virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;


      private:

	

        /**
         * An enum to describe where the depth dependency of the viscosity is coming from.
         */
        enum ViscositySource
        {
          Function,
          File,
          List,
          None
        };

        /**
         * Currently chosen source for the viscosity.
         */
        ViscositySource viscosity_source;

        /**
         * Function to read depth-dependent lookup table and set up interpolating function
         * for File depth dependence method
         */
       // void
       // read_viscosity_file(const std::string &filename,
        //                    const MPI_Comm &comm);

        /**
         * Data structures to store depth and viscosity lookup tables as well as interpolating
         * function to calculate viscosity for File Depth dependence method
         */
        //std::unique_ptr< Functions::InterpolatedTensorProductGridData<1> > viscosity_file_function;

        /**
         * Function to calculate viscosity at depth using values provided as List input
         */
        //double
        //viscosity_from_list(const double &depth) const;

        /**
         * function to calculate viscosity at depth using values provided from File input
         */
        //double
        //viscosity_from_file(const double &depth) const;

        /**
         * Function to calculate depth-dependent multiplicative prefactor to be applied
         * to base model viscosity.
         */
        //double
        //calculate_depth_dependent_prefactor(const double &depth) const;

        /**
         * Values of depth specified by the user if using List depth dependence method
         */
        std::vector<double> depth_values;
        std::vector<double> viscosity_values;
        /**
         * Parsed function that specifies viscosity depth-dependence when using the Function
         * method.
         */
        //Functions::ParsedFunction<1> viscosity_function;

        /**
         * Pointer to the material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim> > base_model;

	// The following comes from the viscoelastic.h file.

	/**
         * Reference temperature for thermal expansion. All components use
         * the same reference_T.
         */
        double reference_T;
	double reference_permeability;
	double eta_f;
        /**
         * Enumeration for selecting which viscosity averaging scheme to use.
         */
        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

        /**
         * Used for calculating average elastic shear modulus and viscosity
         */
        double calculate_average_vector (const std::vector<double> &composition,
                                         const std::vector<double> &parameter_values,
                                         const MaterialUtilities::CompositionalAveragingOperation &average_type) const;



	double calculate_average_viscoelastic_viscosity (const double average_viscosity,
                                                 const double average_elastic_shear_modulus,
                                                 const double dte) const;


  	/**
         * Vector for field densities, read from parameter file.
         */
        std::vector<double> densities;

        /**
         * Vector for field viscosities, read from parameter file.
         */
        std::vector<double> viscosities;

        /**
         * Vector for field thermal expnsivities, read from parameter file.
         */
        std::vector<double> thermal_expansivities;

        /**
         * Vector for field thermal conductivities, read from parameter file.
         */
        std::vector<double> thermal_conductivities;

        /**
         * Vector for field specific heats, read from parameter file.
         */
        std::vector<double> specific_heats;

        /**
         * Vector for field elastic shear moduli, read from parameter file.
         */
        std::vector<double> elastic_shear_moduli;

        /**
         * Bool indicating whether to use a fixed material time scale in the
         * viscoelastic rheology for all time steps (if true) or to use the
         * actual (variable) advection time step of the model (if false). Read
         * from parameter file.
         */
        bool use_fixed_elastic_time_step;

        /**
         * Bool indicating whether to use a stress averaging scheme to account
         * for differences between the numerical and fixed elastic time step
         * (if true). When set to false, the viscoelastic stresses are not
         * modified to account for differences between the viscoelastic time
         * step and the numerical time step. Read from parameter file.
         */
        bool use_stress_averaging;

	 /**
         * Double for fixed elastic time step value, read from parameter file
         */
        double fixed_elastic_time_step;

    };
  }
}

#endif
