/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/material_model/melt_viscoelastic.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <array>

#include <utility>
#include <limits>
//#include <aspect/material_model/viscoelastic.h>

namespace aspect
{
  namespace MaterialModel
  {

// From viscoelastic model

	namespace
    	{
      	  std::vector<std::string> make_elastic_additional_outputs_names()
      	  {
        	std::vector<std::string> names;
        	names.emplace_back("elastic_shear_modulus");
        	return names;
      	  }
    	}

	template <int dim>
    	  ElasticAdditionalOutputs<dim>::ElasticAdditionalOutputs (const unsigned int n_points)
      	  :
      	  NamedAdditionalMaterialOutputs<dim>(make_elastic_additional_outputs_names()),
      	  elastic_shear_moduli(n_points, numbers::signaling_nan<double>())
    	  {}

    	template <int dim>
    	std::vector<double>
    	ElasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    	{
      	  AssertIndexRange (idx, 1);
      	  switch (idx)
            {
          	case 0:
            	  return elastic_shear_moduli;

          	default:
            	  AssertThrow(false, ExcInternalError());
            }
      	// we will never get here, so just return something
      	return elastic_shear_moduli;
    	}






	template <int dim>
    	double
    	MeltViscoelastic<dim>::
    	reference_darcy_coefficient () const
    	{
      	// 0.01 = 1% melt
      	return reference_permeability * std::pow(0.01,3.0) / eta_f;
    	}







// There before

    template <int dim>
    void
    MeltViscoelastic<dim>::initialize()
    {
      base_model->initialize();
    }




    template <int dim>
    double
    MeltViscoelastic<dim>::
    calculate_average_vector (const std::vector<double> &composition,
                              const std::vector<double> &parameter_values,
                              const MaterialUtilities::CompositionalAveragingOperation &average_type) const
    {
      // Store which components to exclude during volume fraction computation.
      ComponentMask composition_mask(this->n_compositional_fields(), true);
      // assign compositional fields associated with viscoelastic stress a value of 0
      // assume these fields are listed first
      for (unsigned int i=0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
        composition_mask.set(i, false);
      const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition, composition_mask);
	// The problem is that these are not the same size ?
        Assert(volume_fractions.size() == parameter_values.size(),
               ExcMessage ("The volume fractions and parameter values vectors used for averaging "
                           "have to have the same length! inside calculate_average_vector function"));
      const double averaged_vector = MaterialUtilities::average_value(volume_fractions, parameter_values, average_type);
      return averaged_vector;
    }

    template <int dim>
    double
    MeltViscoelastic<dim>::
    calculate_average_viscoelastic_viscosity (const double average_viscosity,
                                              const double average_elastic_shear_modulus,
                                              const double dte) const
    {
      return ( average_viscosity * dte ) / ( dte + ( average_viscosity / average_elastic_shear_modulus ) );
    }


// This is probably useless ?
//    template <int dim>
//    void
//    MeltViscoelastic<dim>::update()
//    {
//      base_model->update();
//    }

// Function to compute viscoelasticity.
    template <int dim>
    void
    MeltViscoelastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {

      base_model->evaluate(in,out);
//	cout << "Entering evaluate in MeltViscoelastic";

      // Create the structure for the elastic force terms that are needed to compute the
      // right-hand side of the Stokes system
      MaterialModel::ElasticOutputs<dim>
      *force_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim> >();


      // Store which components to exclude during volume fraction computation.
      ComponentMask composition_mask(this->n_compositional_fields(),true);
      // assign compositional fields associated with viscoelastic stress a value of 0
      // assume these fields are listed first
      for (unsigned int i=0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
        composition_mask.set(i,false);

      // The elastic time step (dte) is equal to the numerical time step if the time step number
      // is greater than 0 and the parameter 'use_fixed_elastic_time_step' is set to false.
      // On the first (0) time step the elastic time step is always equal to the value
      // specified in 'fixed_elastic_time_step', which is also used in all subsequent time
      // steps if 'use_fixed_elastic_time_step' is set to true.
      //
      // We also use this parameter when we are still *before* the first time step,
      // i.e., if the time step number is numbers::invalid_unsigned_int.
      const double dte = ( ( this->get_timestep_number() > 0 &&
                             this->get_timestep_number() != numbers::invalid_unsigned_int &&
                             use_fixed_elastic_time_step == false )
                           ?
                           this->get_timestep()
                           :
                           fixed_elastic_time_step);

        for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = in.temperature[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition, composition_mask);

	  //cout << "Volume fraction size: ";
	  //cout << volume_fractions.size();
	  //cout << "Specific heat size: ";
	  //cout << specific_heats.size();
	  Assert(volume_fractions.size() == specific_heats.size(),
               ExcMessage ("The volume fractions and specific heats vectors used for averaging "
                           "have to have the same length!"));
          out.specific_heat[i] = MaterialUtilities::average_value(volume_fractions, specific_heats, MaterialUtilities::arithmetic);

          // Arithmetic averaging of thermal conductivities
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
        Assert(volume_fractions.size() == thermal_conductivities.size(),
               ExcMessage ("The volume fractions and thermal conductivities vectors used for averaging "
                           "have to have the same length!"));
          out.thermal_conductivities[i] = MaterialUtilities::average_value(volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);

          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // not strictly correct if thermal expansivities are different, since we are interpreting
              // these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }
          out.densities[i] = density;

	        Assert(volume_fractions.size() == thermal_expansivities.size(),
               ExcMessage ("The volume fractions and thermal expansivities vectors used for averaging "
                           "have to have the same length!"));
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value(volume_fractions, thermal_expansivities, MaterialUtilities::arithmetic);

          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          // (here we use an incompressible medium)
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          // Average viscosity
          const double average_viscosity = calculate_average_vector(composition,
                                                                    viscosities,
                                                                    viscosity_averaging);

          // Average elastic shear modulus
          const double average_elastic_shear_modulus = calculate_average_vector(composition,
                                                                                elastic_shear_moduli,
                                                                                viscosity_averaging);

          // Average viscoelastic (e.g., effective) viscosity (equation 28 in Moresi et al., 2003, J. Comp. Phys.)
          out.viscosities[i] = calculate_average_viscoelastic_viscosity(average_viscosity,
                                                                        average_elastic_shear_modulus,
                                                                        dte);

// THIS is part I added (copied from melt_global.cc). Not sure if correct conditional statement for this script.
          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));

              // calculate viscosity based on local melt
              out.viscosities[i] *= (1.0 - porosity);
            }
//END
          //out.viscosities[i] *= in.porosity[i];
          // Fill the material properties that are part of the elastic additional outputs
          if (ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<ElasticAdditionalOutputs<dim> >())
            {
              elastic_out->elastic_shear_moduli[i] = average_elastic_shear_modulus;
            }

          // Fill elastic force outputs (assumed to be zero during intial time step)
          if (force_out)
            {
              force_out->elastic_force[i] = 0.;
            }

	}


      // Viscoelasticity section
      if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 && in.strain_rate.size() > 0)
        {
          // Get old (previous time step) velocity gradients
          std::vector<Point<dim> > quadrature_positions(in.position.size());
          for (unsigned int i=0; i < in.position.size(); ++i)
            quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

          // FEValues requires a quadrature and we provide the default quadrature
          // as we only need to evaluate the solution and gradients.
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   Quadrature<dim>(quadrature_positions),
                                   update_gradients);

          fe_values.reinit (in.current_cell);
          std::vector<Tensor<2,dim> > old_velocity_gradients (quadrature_positions.size(), Tensor<2,dim>());
          fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_old_solution(),
                                                                                         old_velocity_gradients);

          MaterialModel::ElasticOutputs<dim>
          *force_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim> >();


        for (unsigned int i=0; i < in.position.size(); ++i)
            {
              // Get old stresses from compositional fields
              SymmetricTensor<2,dim> stress_old;
              for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];

              // Calculate the rotated stresses
              // Rotation (vorticity) tensor (equation 25 in Moresi et al., 2003, J. Comp. Phys.)
              const Tensor<2,dim> rotation = 0.5 * ( old_velocity_gradients[i] - transpose(old_velocity_gradients[i]) );

              // Recalculate average elastic shear modulus
              const std::vector<double> composition = in.composition[i];
              const double average_elastic_shear_modulus = calculate_average_vector(composition,
                                                                                    elastic_shear_moduli,
                                                                                    viscosity_averaging);

              // Average viscoelastic viscosity
              const double average_viscoelastic_viscosity = out.viscosities[i];

              // Calculate the current (new) viscoelastic stress, which is a function of the material
              // properties (viscoelastic viscosity, shear modulus), elastic time step size, strain rate,
              // vorticity and prior (inherited) viscoelastic stresses (see equation 29 in Moresi et al.,
              // 2003, J. Comp. Phys.)
              SymmetricTensor<2,dim> stress_new = ( 2. * average_viscoelastic_viscosity * deviator(in.strain_rate[i]) ) +
                                                  ( ( average_viscoelastic_viscosity / ( average_elastic_shear_modulus * dte ) ) * stress_old ) +
                                                  ( ( average_viscoelastic_viscosity / average_elastic_shear_modulus ) *
                                                    ( symmetrize(rotation * Tensor<2,dim>(stress_old) ) - symmetrize(Tensor<2,dim>(stress_old) * rotation) ) );




              // Stress averaging scheme to account for difference betweed fixed elastic time step
              // and numerical time step (see equation 32 in Moresi et al., 2003, J. Comp. Phys.)
              const double dt = this->get_timestep();
              if (use_fixed_elastic_time_step == true && use_stress_averaging == true)
                {
                  stress_new = ( ( 1. - ( dt / dte ) ) * stress_old ) + ( ( dt / dte ) * stress_new ) ;
                }

              // Fill reaction terms
              for (unsigned int j = 0; j < SymmetricTensor<2,dim>::n_independent_components ; ++j)
                out.reaction_terms[i][j] = -stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)]
                                           + stress_new[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)];

              // Fill elastic force outputs (See equation 30 in Moresi et al., 2003, J. Comp. Phys.)
              if (force_out)
                {
                  force_out->elastic_force[i] = -1. * ( ( average_viscoelastic_viscosity / ( average_elastic_shear_modulus * dte  ) ) * stress_old );
                }

	   }
	}

    }



// porosity correction ? done above.
//    template <int dim>
//    void
//    MeltViscoelastic<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
//                                  typename Interface<dim>::MaterialModelOutputs &out) const
//    MeltViscoelastic<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
//                                  typename Interface<dim>::MaterialModelOutputs &out,
//				  const MaterialModel::MaterialModelInputs<dim> &inVE,
//				  MaterialModel::MaterialModelOutputs<dim> &outVE) const
//    {
      
// Should this go in the evaluate_VE function ?
//      base_model->evaluate(in,out);

//      if (in.strain_rate.size())
//        {
          // Scale the base model viscosity value by the depth dependent prefactor
//          for (unsigned int i=0; i < out.viscosities.size(); ++i)
//            {
    //          const double depth = this->get_geometry_model().depth(in.position[i]);
    //          out.viscosities[i] *= calculate_depth_dependent_prefactor( depth );
    //            out.viscosities[i] *= (1 - out.porosity[i]);
//		out.viscosities[i] = MeltViscoelastic<dim>::calculate_average_viscoelastic_viscosity(out.viscosities[i],100,100);
//		out.viscosities[i] = Viscoelastic<dim>::calculate_average_viscoelastic_viscosity(out.viscosities[i],100,100);
//		out.viscosities[i] = Viscoelastic<dim>::is_compressible();
//            }
//        }
//    }

// Leave the thermodynamic parameters/viscosities/densities to the melt (base) model.

    template <int dim>
    void
    MeltViscoelastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt viscoelastic model");
        {
          prm.declare_entry("Base model","melt global",
			    Patterns::Selection("melt simple|melt global"),
//                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by a depth "
                            "dependent viscosity. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.declare_entry ("Depth dependence method", "None",
                             Patterns::Selection("None"),
                             "Method that is used to specify how the viscosity should vary with depth. ");
	  prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $\\text{K}$.");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $Pa s$");
          prm.declare_entry ("Thermal expansivities", "4.e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $1/K$");
          prm.declare_entry ("Specific heats", "1250.",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heats $C_p$ for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $J /kg /K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $W/m/K$ ");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition "),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Elastic shear moduli", "75.0e9",
                             Patterns::List(Patterns::Double(0)),
                             "List of elastic shear moduli, $G$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The default value of 75 GPa is representative of mantle rocks. Units: Pa.");
	  prm.declare_entry ("Use fixed elastic time step", "unspecified",
                             Patterns::Selection("true|false|unspecified"),
                             "Select whether the material time scale in the viscoelastic constitutive "
                             "relationship uses the regular numerical time step or a separate fixed "
                             "elastic time step throughout the model run. The fixed elastic time step "
                             "is always used during the initial time step. If a fixed elastic time "
                             "step is used throughout the model run, a stress averaging scheme can be "
                             "applied to account for differences with the numerical time step. An "
                             "alternative approach is to limit the maximum time step size so that it "
                             "is equal to the elastic time step. The default value of this parameter is "
                             "'unspecified', which throws an exception during runtime. In order for "
                             "the model to run the user must select 'true' or 'false'.");
          prm.declare_entry ("Fixed elastic time step", "1.e3",
                             Patterns::Double (0),
                             "The fixed elastic time step $dte$. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Use stress averaging","false",
                             Patterns::Bool (),
                             "Whether to apply a stress averaging scheme to account for differences "
                             "between the fixed elastic time step and numerical time step. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltViscoelastic<dim>::parse_parameters (ParameterHandler &prm)
    {

      // Get the number of fields for composition-dependent material properties
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt viscoelastic model");
        {
	  reference_T = prm.get_double ("Reference temperature");

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          // Parse viscoelastic properties
          densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                              n_fields,
                                                              "Densities");
          viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosities"))),
                                                                n_fields,
                                                                "Viscosities");
          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");
          thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                          n_fields,
                                                                          "Thermal expansivities");
          specific_heats = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Specific heats"))),
                                                                   n_fields,
                                                                   "Specific heats");
          AssertThrow( prm.get("Base model") != "melt viscoelastic",
                       ExcMessage("You may not use ``melt viscoelastic'' as the base model for "
                                  "a melt-viscoelastic model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          if (  prm.get("Depth dependence method") == "None" )
            viscosity_source = None;
          else
            {
              AssertThrow(false, ExcMessage("Unknown method for depth dependence."));
            }


	//  Assign different shear modulus for porosity field...
	  elastic_shear_moduli = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Elastic shear moduli"))),
                                                                         n_fields,
                                                                         "Elastic shear moduli");

	  if (prm.get ("Use fixed elastic time step") == "true")
            use_fixed_elastic_time_step = true;
          else if (prm.get ("Use fixed elastic time step") == "false")
            use_fixed_elastic_time_step = false;
          else
            AssertThrow(false, ExcMessage("'Use fixed elastic time step' must be set to 'true' or 'false'"));

          use_stress_averaging = prm.get_bool ("Use stress averaging");
          if (use_stress_averaging)
            AssertThrow(use_fixed_elastic_time_step == true,
                        ExcMessage("Stress averaging can only be used if 'Use fixed elastic time step' is set to true'"));

          fixed_elastic_time_step = prm.get_double ("Fixed elastic time step");
          AssertThrow(fixed_elastic_time_step > 0,
                      ExcMessage("The fixed elastic time step must be greater than zero"));

          if (this->convert_output_to_years())
            fixed_elastic_time_step *= year_in_seconds;

          AssertThrow(this->get_parameters().enable_elasticity == true,
                      ExcMessage ("Material model Viscoelastic only works if 'Enable elasticity' is set to true"));

	// Make sure this doesn't conflict with porosity and peridotite fields from the melt model.
	// Check whether the compositional fields representing the viscoelastic
          // stress tensor are both named correctly and listed in the right order.
          if (dim == 2)
            {
              AssertThrow(this->introspection().compositional_index_for_name("stress_xx") == 0,
                          ExcMessage("Material model Viscoelastic only works if the first "
                                     "compositional field is called stress_xx."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_yy") == 1,
                          ExcMessage("Material model Viscoelastic only works if the second "
                                     "compositional field is called stress_yy."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_xy") == 2,
                          ExcMessage("Material model Viscoelastic only works if the third "
                                     "compositional field is called stress_xy."));
            }
          else if (dim == 3)
            {
              AssertThrow(this->introspection().compositional_index_for_name("stress_xx") == 0,
                          ExcMessage("Material model Viscoelastic only works if the first "
                                     "compositional field is called stress_xx."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_yy") == 1,
                          ExcMessage("Material model Viscoelastic only works if the second "
                                     "compositional field is called stress_yy."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_zz") == 2,
                          ExcMessage("Material model Viscoelastic only works if the third "
                                     "compositional field is called stress_zz."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_xy") == 3,
                          ExcMessage("Material model Viscoelastic only works if the fourth "
                                     "compositional field is called stress_xy."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_xz") == 4,
                          ExcMessage("Material model Viscoelastic only works if the fifth "
                                     "compositional field is called stress_xz."));
              AssertThrow(this->introspection().compositional_index_for_name("stress_yz") == 5,
                          ExcMessage("Material model Viscoelastic only works if the sixth "
                                     "compositional field is called stress_yz."));
            }
          else
            ExcNotImplemented();


          // Currently, it only makes sense to use this material model when the nonlinear solver
          // scheme does a single Advection iteration and at minimum one Stokes iteration. More
          // than one nonlinear Advection iteration will produce an unrealistic build-up of
          // viscoelastic stress, which are tracked through compositional fields.
          AssertThrow((this->get_parameters().nonlinear_solver ==
                       Parameters<dim>::NonlinearSolver::single_Advection_single_Stokes
                       ||
                       this->get_parameters().nonlinear_solver ==
                       Parameters<dim>::NonlinearSolver::single_Advection_iterated_Stokes),
                      ExcMessage("The material model will only work with the nonlinear "
                                 "solver schemes 'single Advection, single Stokes' and "
                                 "'single Advection, iterated Stokes'"));

          // Functionality to average the additional RHS terms over the cell is not implemented.
          // This enforces that the variable 'Material averaging' is set to 'none'.
          AssertThrow((this->get_parameters().material_averaging ==
                       MaterialModel::MaterialAveraging::none),
                      ExcMessage("The viscoelastic material model cannot be used with "
                                 "material averaging. The variable 'Material averaging' "
                                 "in the 'Material model' subsection must be set to 'none'."));
	//if (this->get_parameters().include_melt_transport)
	//	AssertThrow(false,ExcMessage("Include melt transport true."));

	}
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for depth dependent, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();

      // Declare dependencies on solution variables (from viscoelastic code). Is this necessary ?
      //this->model_dependence.viscosity = NonlinearDependence::compositional_fields;
      //this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      //this->model_dependence.compressibility = NonlinearDependence::none;
      //this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      //this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }

    template <int dim>
    bool
    MeltViscoelastic<dim>::
    is_compressible () const
    {
      //return base_model->is_compressible();
	return false;
    }

    template <int dim>
    double
    MeltViscoelastic<dim>::
    reference_viscosity() const
    {
      /* The viscosity returned by the base model is normalized by the base model's
       * reference viscosity and then scaled by the depth-dependency,
       * so the reference viscosity returned by the depth-dependent
       * model should be representative of the product of the depth dependent contribution
       * and the base model contribution to the total viscosity */
     // const double mean_depth = 0.5*this->get_geometry_model().maximal_depth();
      return base_model->reference_viscosity();
    }



    template <int dim>
    void
    MeltViscoelastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<ElasticAdditionalOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::ElasticAdditionalOutputs<dim>> (n_points));
        }
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class ElasticAdditionalOutputs<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

    ASPECT_REGISTER_MATERIAL_MODEL(MeltViscoelastic,
                                   "melt viscoelastic",
                                   "The ``melt viscoelastic'' Material model applies a depth-dependent scaling "
                                   "to any of the other available material models. In other words, it "
                                   "is a ``compositing material model''."
                                   "\n\n"
                                   "Parameters related to the depth dependent model are read from a subsection "
                                   "``Material model/Melt viscoelastic model''. "
                                   "The user must specify a ``Base model'' from which material properties are "
                                   "derived. Currently the depth dependent model only allows depth dependence of "
                                   "viscosity - other material properties are taken from the ``Base model''. "
                                   "Viscosity $\\eta$ at depth $z$ is calculated according to:"
                                   "\\begin{equation}"
                                   "\\eta(z,p,T,X,...) = \\eta(z) \\eta_b(p,T,X,..)/\\eta_{rb}"
                                   "\\end{equation}"
                                   "where $\\eta(z)$ is the depth-dependence specified by the depth dependent "
                                   "model, $\\eta_b(p,T,X,...)$ is the viscosity calculated from the base model, "
                                   "and $\\eta_{rb}$ is the reference viscosity of the ``Base model''. "
                                   "In addition to the specification of the ``Base model'', the user must specify "
                                   "the method to be used to calculate the depth-dependent viscosity $\\eta(z)$ as "
                                   "``Material model/Melt viscoelastic model/Depth dependence method'', which can be "
                                   "chosen among ``None|Function|File|List''. Each method and the associated parameters "
                                   "are as follows:"
                                   "\n"
                                   "\n"
                                   "``Function'': read a user-specified parsed function from the input file in a "
                                   "subsection ``Material model/Melt viscoelastic model/Viscosity depth function''. "
                                   "By default, this function is uniformly equal to 1.0e21. Specifying a function "
                                   "that returns a value less than or equal to 0.0 anywhere in the model domain will "
                                   "produce an error. "
                                   "\n"
                                   "\n"
                                   "``File'': read a user-specified file containing viscosity values at specified "
                                   "depths. The file containing depth-dependent viscosities is read from a "
                                   "directory specified by the user as "
                                   "``Material model/Melt viscoelastic model/Data directory'', from a file with name "
                                   "specified as ``Material model/Melt viscoelastic model/Viscosity depth file''. "
                                   "The format of this file is ascii text and contains two columns with one header line:"
                                   "\n"
                                   "\n"
                                   "example Viscosity depth file:\\\\"
                                   "Depth (m)    Viscosity (Pa-s)\\\\"
                                   "0.0000000e+00     1.0000000e+21\\\\"
                                   "6.7000000e+05     1.0000000e+22\\\\"
                                   "\n"
                                   "\n"
                                   "Viscosity is interpolated from this file using linear interpolation. "
                                   "``None'': no depth-dependence. Viscosity is taken directly from ``Base model''"
                                   "\n"
                                   "\n"
                                   "``List:'': read a comma-separated list of depth values corresponding to the maximum "
                                   "depths of layers having constant depth-dependence $\\eta(z)$. The layers must be "
                                   "specified in order of increasing depth, and the last layer in the list must have a depth "
                                   "greater than or equal to the maximal depth of the model. The list of layer depths is "
                                   "specified as ``Material model/Melt viscoelastic model/Depth list'' and the corresponding "
                                   "list of layer viscosities is specified as "
                                   "``Material model/Melt viscoelastic model/Viscosity list''"
                                  )
  }
}
