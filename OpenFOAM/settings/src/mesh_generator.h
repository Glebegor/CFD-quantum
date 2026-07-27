#pragma once

#include <string>
#include <vector>

namespace meshgen {

struct Vec2
{
    double x;
    double y;
};

struct Config
{
    // Match the LBM's 512 nodes: OpenFOAM has one cell between each
    // adjacent node, hence 511 cells per axis.
    int resolution = 512;
    int nx = resolution - 1;
    int ny = resolution - 1;

    std::string obstacle_type = "wing";

    // Chord as a fraction of the unit domain and its physical size.
    double chord = 0.35;
    double physical_chord = 4.2;

    [[nodiscard]] double length_scale() const
    {
        return physical_chord / chord;
    }

    [[nodiscard]] double domain_width() const
    {
        return length_scale();
    }

    [[nodiscard]] double domain_height() const
    {
        return length_scale();
    }

    [[nodiscard]] double cell_size() const
    {
        return length_scale() / static_cast<double>(resolution - 1);
    }

    std::string wing_profile_file = "wing_naca_0012.txt";

    // Normalized obstacle position, matching the LBM unit-square parameters.
    // For a wing this is the leading edge; for a circle it is the centre.
    double obstacle_x = 0.25;
    double obstacle_y = 0.50;
    double circle_radius_fraction = 0.08;

    // Positive means nose-up, matching the LBM convention.
    double angle_of_attack_deg = 10.0;

    double airflow_x = 1.0;
    double airflow_y = 0.0;
    double inlet_velocity = 1.0;

    double reynolds_number = 100.0;
    double lattice_inlet_velocity = 0.05;
    double frames_per_second = 10.0;
    double duration_seconds = 4.0;
    int log_interval = 500;

    std::string output_directory = "output";
};

void generate_mesh(const Config& config);

} // namespace meshgen
