#include "Wave2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
  class EigenModeSolution : public Function<Wave2D::dim>
  {
  public:
    explicit EigenModeSolution(const double omega)
      : Function<Wave2D::dim>(1)
      , omega(omega)
    {}

    double
    value(const Point<Wave2D::dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]) *
             std::cos(omega * this->get_time());
    }

  private:
    const double omega;
  };

  class EigenModeVelocity : public Function<Wave2D::dim>
  {
  public:
    explicit EigenModeVelocity(const double omega)
      : Function<Wave2D::dim>(1)
      , omega(omega)
    {}

    double
    value(const Point<Wave2D::dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return -omega * std::sin(numbers::PI * p[0]) *
             std::sin(numbers::PI * p[1]) *
             std::sin(omega * this->get_time());
    }

  private:
    const double omega;
  };

  std::string
  sanitize_for_path(std::string s)
  {
    for (char &c : s)
      if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || c == '-' || c == '_' || c == '.'))
        c = '_';
    return s;
  }
} // namespace

Wave2D::Wave2D(const std::string &mesh_file_name_)
  : Wave2D(mesh_file_name_, SolverOptions())
{}

Wave2D::Wave2D(const std::string   &mesh_file_name_,
               const SolverOptions &options_)
  : mesh_file_name(mesh_file_name_)
  , options(options_)
  , dof_handler(mesh)
{}

void
Wave2D::setup()
{
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh);

    std::ifstream grid_in_file(mesh_file_name);
    if (!grid_in_file)
      throw std::runtime_error("Cannot open mesh file: " + mesh_file_name);

    grid_in.read_msh(grid_in_file);
  }

  {
    fe = std::make_unique<FE_SimplexP<dim>>(options.fe_degree);

    const unsigned int quad_order = std::max(2u, options.fe_degree + 1);
    quadrature                  = std::make_unique<QGaussSimplex<dim>>(quad_order);
  }

  {
    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);
  }

  {
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    stiffness_matrix.reinit(sparsity_pattern);
    mass_matrix_dirichlet.reinit(sparsity_pattern);
    mass_lumped.reinit(dof_handler.n_dofs());
  }

  boundary_dofs = DoFTools::extract_boundary_dofs(dof_handler);
}

void
Wave2D::assemble_matrices()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients | update_JxW_values);

  FullMatrix<double>                   cell_mass(dofs_per_cell, dofs_per_cell);
  FullMatrix<double>                   cell_stiffness(dofs_per_cell,
                                                      dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  mass_matrix      = 0.0;
  stiffness_matrix = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_mass      = 0.0;
      cell_stiffness = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              cell_mass(i, j) += fe_values.shape_value(i, q) *
                                 fe_values.shape_value(j, q) *
                                 fe_values.JxW(q);

              cell_stiffness(i, j) += fe_values.shape_grad(i, q) *
                                      fe_values.shape_grad(j, q) *
                                      fe_values.JxW(q);
            }

      cell->get_dof_indices(dof_indices);
      mass_matrix.add(dof_indices, cell_mass);
      stiffness_matrix.add(dof_indices, cell_stiffness);
    }
}

void
Wave2D::compute_lumped_mass()
{
  mass_lumped = 0.0;
  for (unsigned int i = 0; i < mass_matrix.m(); ++i)
    {
      double s = 0.0;
      for (SparseMatrix<double>::const_iterator it = mass_matrix.begin(i);
           it != mass_matrix.end(i);
           ++it)
        s += it->value();

      mass_lumped[i] = s;
    }

  for (const auto i : boundary_dofs)
    mass_lumped[i] = 1.0;
}

void
Wave2D::enforce_homogeneous_dirichlet_rows(SparseMatrix<double> &matrix) const
{
  for (const auto i : boundary_dofs)
    {
      for (auto it = matrix.begin(i); it != matrix.end(i); ++it)
        it->value() = 0.0;
      matrix.set(i, i, 1.0);
    }
}

void
Wave2D::apply_inverse_mass(Vector<double>       &dst,
                           const Vector<double> &rhs) const
{
  dst.reinit(rhs.size());

  if (options.mass_type == MassType::lumped)
    {
      for (unsigned int i = 0; i < dst.size(); ++i)
        dst[i] = rhs[i] / mass_lumped[i];
    }
  else
    {
      if (!mass_solver_initialized)
        throw std::runtime_error(
          "Internal error: consistent mass solver not initialized.");
      mass_solver.vmult(dst, rhs);
    }

  for (const auto i : boundary_dofs)
    dst[i] = 0.0;
}

void
Wave2D::apply_dirichlet(Vector<double> &u, const double time) const
{
  (void)time;
  for (const auto i : boundary_dofs)
    u[i] = 0.0;
}

std::string
Wave2D::method_tag() const
{
  std::string time_tag = "cd";
  if (options.time_scheme == TimeScheme::newmark)
    time_tag = "newmark";

  const std::string mass_tag =
    (options.mass_type == MassType::lumped ? "lumped" : "consistent");

  return time_tag + "-" + mass_tag + "-p" + std::to_string(options.fe_degree);
}

const std::filesystem::path &
Wave2D::output_directory() const
{
  return run_output_dir;
}

void
Wave2D::output(const Vector<double> &u,
               const unsigned int    step,
               const double          time) const
{
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, u, "u");
  data_out.build_patches();

  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string output_file_name = "output-" + mesh_path.stem().string() +
                                       "-" + method_tag() + "-" +
                                       std::to_string(step) + ".vtu";
  const std::filesystem::path output_path = run_output_dir / output_file_name;
  std::ofstream               output_file(output_path.string());
  data_out.write_vtu(output_file);

  std::cout << "  step " << step << " | t = " << time << " -> "
            << output_path.string() << std::endl;
}

double
Wave2D::compute_energy(const Vector<double> &u, const Vector<double> &v) const
{
  Vector<double> Mv(v.size());
  Vector<double> Ku(u.size());
  mass_matrix.vmult(Mv, v);
  stiffness_matrix.vmult(Ku, u);

  const double kinetic   = 0.5 * (v * Mv);
  const double potential = 0.5 * (u * Ku);
  return kinetic + potential;
}

double
Wave2D::compute_L2_error(const Vector<double> &u_h,
                         const Function<dim>  &u_exact) const
{
  Vector<double> difference(mesh.n_active_cells());

  FE_SimplexP<dim>   fe_map(options.fe_degree);
  MappingFE<dim>     mapping(fe_map);
  QGaussSimplex<dim> quad_err(std::max(3u, options.fe_degree + 2));

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    u_h,
                                    u_exact,
                                    difference,
                                    quad_err,
                                    VectorTools::L2_norm);

  return VectorTools::compute_global_error(mesh,
                                           difference,
                                           VectorTools::L2_norm);
}

double
Wave2D::compute_H1_error(const Vector<double> &u_h,
                         const Function<dim>  &u_exact) const
{
  Vector<double> difference(mesh.n_active_cells());

  FE_SimplexP<dim>   fe_map(options.fe_degree);
  MappingFE<dim>     mapping(fe_map);
  QGaussSimplex<dim> quad_err(std::max(3u, options.fe_degree + 2));

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    u_h,
                                    u_exact,
                                    difference,
                                    quad_err,
                                    VectorTools::H1_seminorm);

  return VectorTools::compute_global_error(mesh,
                                           difference,
                                           VectorTools::H1_seminorm);
}

double
Wave2D::compute_min_cell_diameter() const
{
  double h_min = std::numeric_limits<double>::infinity();
  for (const auto &cell : mesh.active_cell_iterators())
    h_min = std::min(h_min, cell->diameter());
  return h_min;
}

void
Wave2D::run(const double       dt,
            const double       T,
            const unsigned int output_every,
            const double       omega,
            const bool         compute_error_each_step)
{
  if (dt <= 0.0)
    throw std::runtime_error("dt must be positive.");
  if (T < 0.0)
    throw std::runtime_error("T must be non-negative.");

  Timer timer_total;
  timer_total.start();

  setup();
  assemble_matrices();
  compute_lumped_mass();

  if (options.mass_type == MassType::consistent)
    {
      mass_matrix_dirichlet.copy_from(mass_matrix);
      enforce_homogeneous_dirichlet_rows(mass_matrix_dirichlet);
      mass_solver.initialize(mass_matrix_dirichlet);
      mass_solver_initialized = true;
    }

  const double h_min = compute_min_cell_diameter();

  std::cout << "===============================================" << std::endl;
  std::cout << "  Mode: eigenmode" << std::endl;
  std::string time_scheme_name = "central_difference";
  if (options.time_scheme == TimeScheme::newmark)
    time_scheme_name = "newmark";
  std::cout << "  Time scheme: " << time_scheme_name << std::endl;
  std::cout << "  Mass type: "
            << (options.mass_type == MassType::lumped ? "lumped"
                                                      : "consistent")
            << std::endl;
  std::cout << "  FE degree: " << options.fe_degree << std::endl;
  std::cout << "  DoFs:  " << dof_handler.n_dofs() << std::endl;
  std::cout << "  h_min: " << h_min << std::endl;
  std::cout << "  dt:    " << dt << " | T: " << T << std::endl;
  std::cout << "  CFL (dt/h_min): " << dt / h_min << std::endl;
  std::cout << "===============================================" << std::endl;

  auto u_exact_ptr = std::make_unique<EigenModeSolution>(omega);
  auto v_exact_ptr = std::make_unique<EigenModeVelocity>(omega);

  const unsigned int n_steps = static_cast<unsigned int>(std::ceil(T / dt));

  const std::filesystem::path mesh_path(mesh_file_name);
  std::ostringstream          cfg_name;
  cfg_name << "mesh-" << mesh_path.stem().string() << "-mode-eigenmode"
           << "-time-" << method_tag().substr(0, method_tag().find('-'))
           << "-mass-" << (options.mass_type == MassType::lumped ? "lumped" :
                           "consistent")
           << "-p" << options.fe_degree << "-dt-" << std::setprecision(8) << dt
           << "-T-" << T << "-out-" << output_every << "-omega-" << omega
           << "-errstep-" << (compute_error_each_step ? 1 : 0);

  run_output_dir =
    std::filesystem::path(PROJECT_SOURCE_DIR) / "result" /
    sanitize_for_path(cfg_name.str());
  std::filesystem::create_directories(run_output_dir);
  std::cout << "  Output directory: " << run_output_dir.string() << std::endl;

  std::ofstream energy_file(
    (run_output_dir / ("energy-" + method_tag() + ".csv")).string());
  energy_file << "step,time,energy\n";

  std::ofstream error_file((run_output_dir / ("error-" + method_tag() + ".csv"))
                             .string());
  error_file << "step,time,L2_error,H1_error\n";

  Vector<double> u_n(dof_handler.n_dofs());
  u_exact_ptr->set_time(0.0);
  VectorTools::interpolate(dof_handler, *u_exact_ptr, u_n);
  apply_dirichlet(u_n, 0.0);

  Vector<double> v_n(dof_handler.n_dofs());
  v_exact_ptr->set_time(0.0);
  VectorTools::interpolate(dof_handler, *v_exact_ptr, v_n);
  for (const auto i : boundary_dofs)
    v_n[i] = 0.0;

  Vector<double> Ku(dof_handler.n_dofs());
  stiffness_matrix.vmult(Ku, u_n);
  Ku *= -1.0;

  Vector<double> a_n(dof_handler.n_dofs());
  apply_inverse_mass(a_n, Ku);

  Timer timer_stepping;
  timer_stepping.start();

  if (options.time_scheme == TimeScheme::central_difference)
    {
      Vector<double> u_nm1(dof_handler.n_dofs());
      Vector<double> u_np1(dof_handler.n_dofs());

      u_nm1 = u_n;
      u_n.add(dt, v_n);
      u_n.add(0.5 * dt * dt, a_n);
      apply_dirichlet(u_n, dt);

      auto compute_velocity = [&](const Vector<double> &u_curr,
                                  const Vector<double> &u_prev) {
        Vector<double> v(u_curr.size());
        v = u_curr;
        v.add(-1.0, u_prev);
        v *= (1.0 / dt);
        for (const auto idx : boundary_dofs)
          v[idx] = 0.0;
        return v;
      };

      if (output_every > 0)
        output(u_n, 1, dt);

      {
        const Vector<double> v_step = compute_velocity(u_n, u_nm1);
        energy_file << 1 << "," << dt << "," << compute_energy(u_n, v_step)
                    << "\n";
      }

      u_exact_ptr->set_time(dt);
      if (compute_error_each_step)
        error_file << 1 << "," << dt << ","
                   << compute_L2_error(u_n, *u_exact_ptr) << ","
                   << compute_H1_error(u_n, *u_exact_ptr) << "\n";

      Vector<double> rhs(dof_handler.n_dofs());
      for (unsigned int step = 1; step < n_steps; ++step)
        {
          const double time_n = step * dt;

          stiffness_matrix.vmult(rhs, u_n);
          rhs *= -1.0;

          Vector<double> acc(rhs.size());
          apply_inverse_mass(acc, rhs);

          u_np1 = u_n;
          u_np1 *= 2.0;
          u_np1.add(-1.0, u_nm1);
          u_np1.add(dt * dt, acc);
          apply_dirichlet(u_np1, time_n + dt);

          u_nm1 = u_n;
          u_n   = u_np1;

          const unsigned int out_step = step + 1;
          const double       time_np1 = time_n + dt;

          if (output_every > 0 && (out_step % output_every == 0))
            output(u_n, out_step, time_np1);

          const Vector<double> v_step = compute_velocity(u_n, u_nm1);
          energy_file << out_step << "," << time_np1 << ","
                      << compute_energy(u_n, v_step) << "\n";

          if (compute_error_each_step)
            {
              u_exact_ptr->set_time(time_np1);
              error_file << out_step << "," << time_np1 << ","
                         << compute_L2_error(u_n, *u_exact_ptr) << ","
                         << compute_H1_error(u_n, *u_exact_ptr) << "\n";
            }
        }

      v_n = compute_velocity(u_n, u_nm1);
    }
  else
    {
      constexpr double beta  = 0.25;
      constexpr double gamma = 0.5;

      SparseMatrix<double> system_matrix(sparsity_pattern);
      system_matrix.copy_from(stiffness_matrix);
      system_matrix *= (beta * dt * dt);

      if (options.mass_type == MassType::consistent)
        system_matrix.add(1.0, mass_matrix);
      else
        for (unsigned int i = 0; i < mass_lumped.size(); ++i)
          system_matrix.add(i, i, mass_lumped[i]);

      enforce_homogeneous_dirichlet_rows(system_matrix);

      SparseDirectUMFPACK newmark_solver;
      newmark_solver.initialize(system_matrix);

      Vector<double> u_np1(dof_handler.n_dofs());
      Vector<double> u_pred(dof_handler.n_dofs());
      Vector<double> v_pred(dof_handler.n_dofs());
      Vector<double> rhs(dof_handler.n_dofs());

      if (output_every > 0)
        output(u_n, 0, 0.0);

      energy_file << 0 << "," << 0.0 << "," << compute_energy(u_n, v_n) << "\n";

      if (compute_error_each_step)
        {
          u_exact_ptr->set_time(0.0);
          error_file << 0 << "," << 0.0 << ","
                     << compute_L2_error(u_n, *u_exact_ptr) << ","
                     << compute_H1_error(u_n, *u_exact_ptr) << "\n";
        }

      for (unsigned int step = 0; step < n_steps; ++step)
        {
          const double time_np1 = (step + 1) * dt;

          u_pred = u_n;
          u_pred.add(dt, v_n);
          u_pred.add(dt * dt * (0.5 - beta), a_n);

          v_pred = v_n;
          v_pred.add(dt * (1.0 - gamma), a_n);

          if (options.mass_type == MassType::consistent)
            mass_matrix.vmult(rhs, u_pred);
          else
            {
              rhs = 0.0;
              for (unsigned int i = 0; i < rhs.size(); ++i)
                rhs[i] = mass_lumped[i] * u_pred[i];
            }

          for (const auto i : boundary_dofs)
            rhs[i] = 0.0;

          newmark_solver.vmult(u_np1, rhs);
          apply_dirichlet(u_np1, time_np1);

          a_n = u_np1;
          a_n.add(-1.0, u_pred);
          a_n *= (1.0 / (beta * dt * dt));
          for (const auto i : boundary_dofs)
            a_n[i] = 0.0;

          v_n = v_pred;
          v_n.add(gamma * dt, a_n);
          for (const auto i : boundary_dofs)
            v_n[i] = 0.0;

          u_n = u_np1;

          const unsigned int out_step = step + 1;
          if (output_every > 0 && (out_step % output_every == 0))
            output(u_n, out_step, time_np1);

          energy_file << out_step << "," << time_np1 << ","
                      << compute_energy(u_n, v_n) << "\n";

          if (compute_error_each_step)
            {
              u_exact_ptr->set_time(time_np1);
              error_file << out_step << "," << time_np1 << ","
                         << compute_L2_error(u_n, *u_exact_ptr) << ","
                         << compute_H1_error(u_n, *u_exact_ptr) << "\n";
            }
        }
    }

  timer_stepping.stop();
  timer_total.stop();

  u_exact_ptr->set_time(n_steps * dt);
  const double final_l2 = compute_L2_error(u_n, *u_exact_ptr);
  const double final_h1 = compute_H1_error(u_n, *u_exact_ptr);
  std::cout << "===============================================" << std::endl;
  std::cout << "  Final L2 error at T = " << n_steps * dt << ": "
            << final_l2 << std::endl;
  std::cout << "  Final H1 error at T = " << n_steps * dt << ": "
            << final_h1 << std::endl;
  if (!compute_error_each_step)
    error_file << n_steps << "," << n_steps * dt << "," << final_l2 << ","
               << final_h1 << "\n";

  {
    const double E_final = compute_energy(u_n, v_n);
    std::cout << "  Final energy: " << E_final << std::endl;
  }

  energy_file.close();
  error_file.close();

  std::cout << "===============================================" << std::endl;
  std::cout << "  Time stepping: " << timer_stepping.wall_time() << " s ("
            << n_steps << " steps, "
            << (n_steps > 0 ? timer_stepping.wall_time() / n_steps * 1e6 : 0.0)
            << " us/step)" << std::endl;
  std::cout << "  Total wall time: " << timer_total.wall_time() << " s"
            << std::endl;
  std::cout << "===============================================" << std::endl;
}
