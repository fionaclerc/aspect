#include <aspect/melt.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{

  template <int dim>
  class TestMeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:

      virtual double reference_darcy_coefficient () const
      {
        return 1.0;
      }


      virtual bool
      viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        if ((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
          return true;
        return false;
      }

      virtual bool
      density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual double reference_density () const
      {
        return 1.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {


        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
        *force = out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();
        const double elastic_pore_modulus = 2.e3;
        const double tau = 1.e4;
        const double dte = this->get_timestep();
        const double xi = 1.0;
        enum CompactionPressureFormulation
        {
          strain_rate,
          fluid_pressure,
          compaction_pressure
        } p_c_formulation;

        //p_c_formulation = compaction_pressure;
        p_c_formulation = fluid_pressure;

        std::vector<double> fluid_pressures(in.position.size());
        std::vector<double> compaction_pressures(in.position.size());

        double p_c_scale = 0;
        if (in.current_cell.state() == IteratorState::valid)
          {
            std::vector<Point<dim> > quadrature_positions(in.position.size());
            for (unsigned int i=0; i < in.position.size(); ++i)
                quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

            FEValues<dim> fe_values (this->get_mapping(),
                                     this->get_fe(),
                                     Quadrature<dim>(quadrature_positions),
                                     update_values | update_gradients);
            fe_values.reinit (in.current_cell);
            const FEValuesExtractors::Scalar extractor_pressure = this->introspection().variable("fluid pressure").extractor_scalar();
            fe_values[extractor_pressure].get_function_values (this->get_solution(),
                                                               fluid_pressures);

            const FEValuesExtractors::Scalar extractor_pressure_compaction = this->introspection().variable("compaction pressure").extractor_scalar();
            fe_values[extractor_pressure].get_function_values (this->get_solution(),
                                                               compaction_pressures);

            p_c_scale = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>(this->get_material_model()).p_c_scale(in,
                                 out,
                                 this->get_melt_handler(),
                                 true);

          }

        double xi_ve = xi;
        if (this->get_timestep_number() > 0)
          xi_ve = 1./(1./xi + 1./(elastic_pore_modulus * dte));


        double p_c_err = 0.0;
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const double porosity = in.composition[i][porosity_idx];
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            const double p_c_old = in.composition[i][this->introspection().compositional_index_for_name("p_c_field")];


            double p_c_new = 0.;
            switch (p_c_formulation)
              {
                case fluid_pressure:
                {
                  p_c_new = (1. - porosity) * (in.pressure[i] - fluid_pressures[i]);
                  break;
                }
                case compaction_pressure:
                {
                  p_c_new = compaction_pressures[i] * p_c_scale; 
                  break;
                }
                case strain_rate:
                {
                  p_c_new = xi_ve * ((porosity - 1.0) * trace(in.strain_rate[i]) +  p_c_old/(elastic_pore_modulus * dte));
                  break;
                }
              }
             std::cout << p_c_new;
             std::cout << " \n";

            // get VE compaction viscosity now
            const double viscosity_VE = 1.0;
            out.viscosities[i] = viscosity_VE;
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0;

            out.reaction_terms[i][porosity_idx] = 0;
            out.reaction_terms[i][this->introspection().compositional_index_for_name("p_c_field")] = -p_c_old + p_c_new;
            double p_c_force;
            double p_c_force_sol_RHS;

            const double phi = 0.01 + 0.2 * exp(-20 * pow(x + 2. * z, 2.));

            if (this->get_timestep_number() > 0)
              {
                p_c_force = -1.0 * p_c_old/(elastic_pore_modulus * dte);
                p_c_force_sol_RHS = - (1 - phi) * tau * ( sin(x + z)  * cos(this->get_time()*tau) )/(elastic_pore_modulus );
              }
            p_c_err += p_c_force + p_c_force_sol_RHS;

//**********
// copy and paste here
            if (force)
              {
                force->rhs_u[i][0] = cos(z) + cos(x * z) * z ;// -1. * cos(x + z) * sin(tau * this->get_time());
                force->rhs_u[i][1] = sin(x) + x*cos(x * z) ;//-1. * cos(x + z) * sin(tau * this->get_time());
                force->rhs_p[i] = 2. * sin(x + z) *sin( tau * this->get_time()) - sin(x * z) * z * z - sin(x * z) * x * x;
                force->rhs_melt_pc[i] = p_c_force + p_c_force_sol_RHS -  ((1 -phi) * (sin(x + z) * sin(this->get_time()*tau))) / xi;
              }
//***********
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != nullptr)
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              const double x = in.position[i](0);
              const double z = in.position[i](1);

              melt_out->compaction_viscosities[i] = xi_ve; // xi
              melt_out->fluid_viscosities[i] = 1.0;
              melt_out->permeabilities[i] = 1.0; // K_D
              melt_out->fluid_density_gradients[i] = 0.0;
              melt_out->fluid_densities[i] = 0.5;
            }

      }

  };



  template <int dim>
  class RefFunction : public Function<dim> //, public ::aspect::SimulatorAccess<dim>
  {
    public:
      RefFunction () : Function<dim>(2*dim+3+3) {}
      virtual void vector_value (const Point< dim >   &p,
                                 Vector< double >   &values) const
      {
        double x = p(0);
        double z = p(1);
        double tt = this->get_time();
        double dt = 1.e-5;
        double tau = 1.e4;
        double p_c_0 = sin(x + z) * sin((tt - dt) * tau);
        const double phi = 0.01 + 0.2 * exp(-20 * pow(x + 2. * z, 2.));

//**********
// copy and paste here (add "out.")
        values[0] = cos(z); // u
        values[1] = sin(x); // v
        values[2] = - (sin(x + z) * sin(tt*tau)) + sin(x * z); // p_f
        values[3] = (1 -phi) * (sin(x + z)  * sin(tt*tau)); // p_c
        if (phi <= 0.0)
          {
            values[4] = cos(z);
            values[5] = sin(x);
          } else {
            values[4] = ((cos(x + z) * sin(tt*tau)) - z * cos(x * z))/phi  + cos(z); // uf_x
            values[5] =  ((cos(x + z) * sin(tt*tau)) - x * cos(x * z))/phi  + sin(x); // uf_y
          }
        values[6] = sin(x * z); // p_s
        values[7] = 0;
        values[8] = phi; // phi
        values[9] = (1 - phi) * p_c_0; // p_c_old
//**********
      }
  };



  /**
    * A postprocessor that evaluates the accuracy of the solution
    * by using the L2 norm.
    */
  template <int dim>
  class ConvergenceMeltPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Generate graphical output from the current solution.
       */
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics);

  };

  template <int dim>
  std::pair<std::string,std::string>
  ConvergenceMeltPostprocessor<dim>::execute (TableHandler &statistics)
  {
    
    RefFunction<dim> ref_func;
    const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +2);

    const unsigned int n_total_comp = this->introspection().n_components;

    Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_u_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c_field (this->get_triangulation().n_active_cells());

    ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                        n_total_comp);
    ComponentSelectFunction<dim> comp_p_f(dim, n_total_comp);
    ComponentSelectFunction<dim> comp_p_c(dim+1, n_total_comp);
    ComponentSelectFunction<dim> comp_u_f(std::pair<unsigned int, unsigned int>(dim+2,dim+2+
                                                                                dim),
                                          n_total_comp);
    ComponentSelectFunction<dim> comp_p(dim+2+dim, n_total_comp);
    ComponentSelectFunction<dim> comp_porosity(dim+2+dim+2, n_total_comp);
    ComponentSelectFunction<dim> comp_p_c_field(dim+2+dim+3, n_total_comp);

    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_f);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_c,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_c);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_porosity,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_porosity);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_c_field,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_c_field);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u_f);

    const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u, VectorTools::L2_norm);
    const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L2_norm);
    const double p_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_f, VectorTools::L2_norm);
    const double p_c_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_c, VectorTools::L2_norm);
    const double phi_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_porosity, VectorTools::L2_norm);
    const double p_c_field_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_c_field, VectorTools::L2_norm);
    const double u_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u_f, VectorTools::L2_norm);


    std::ostringstream os;
    os << std::scientific
       << " h = " << this->get_triangulation().begin_active()->diameter()
       << " ndofs= " << this->get_solution().size()
       << " u_L2= " << u_l2
       << " p_L2= "  << p_l2
       << " p_f_L2= " << p_f_l2
       << " p_c_L2= " << p_c_l2
       << " phi_L2= " << phi_l2
       << " p_c_field_L2= " << p_c_field_l2
       << " u_f_L2= " << u_f_l2
       ;



    std::string filename ="output_errors_" + Utilities::int_to_string(this->get_timestep_number(), 5);
    std::ofstream file (filename.c_str());

    file << u_l2 << " "
         << p_l2 << " "
         << p_f_l2 << " "
         << p_c_l2 << " "
         << phi_l2 << " "
         << p_c_field_l2 << " "
         << u_f_l2 << std::endl;
  

    return std::make_pair("Errors", os.str());
  }


  template <int dim>
  class PressureBdry:
    public BoundaryFluidPressure::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const dealii::types::boundary_id,
        const typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
        const std::vector<Tensor<1,dim> > &normal_vectors,
        std::vector<double> &output
      ) const
      {
        const double tau = 1.e4;
        for (unsigned int q=0; q<output.size(); ++q)
          {
            const double x = material_model_inputs.position[q][0];
            const double z = material_model_inputs.position[q][1];
            Tensor<1,dim> gradient;
//**********
// copy and paste here (add "out.")
            gradient[0] = -1. * cos(x + z) * sin(tau * this->get_time()) + cos(x * z) * z;
            gradient[1] = -1. * cos(x + z) * sin(tau * this->get_time())+ cos(x * z) * x;
//**********
            output[q] = gradient * normal_vectors[q];
          }
      }



  };

}

// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_MATERIAL_MODEL(TestMeltMaterial,
                                 "test melt material",
                                 "")

  ASPECT_REGISTER_POSTPROCESSOR(ConvergenceMeltPostprocessor,
                                "melt error calculation",
                                "A postprocessor that compares the numerical solution to the analytical solution")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}
