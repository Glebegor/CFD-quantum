#include "mesh_generator.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace meshgen {
namespace {
namespace fs = std::filesystem;

struct Vec3
{
    double x;
    double y;
    double z;
};

struct FaceRecord
{
    std::vector<int> points;
    int owner = -1;
    int neighbour = -1;
    std::string patch_name;
};

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

bool is_comment_or_empty(const std::string& line)
{
    const std::string trimmed = trim(line);
    if (trimmed.empty()) {
        return true;
    }
    return trimmed[0] == '#';
}

std::vector<Vec2> load_profile(const fs::path& path)
{
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Unable to open profile file: " + path.string());
    }

    std::vector<Vec2> points;
    std::string line;
    while (std::getline(input, line)) {
        if (is_comment_or_empty(line)) {
            continue;
        }

        std::istringstream stream(line);
        double x = 0.0;
        double y = 0.0;
        if (!(stream >> x >> y)) {
            continue;
        }

        // Lednicer airfoil files contain a numeric point-count line (for
        // example "66. 66.") before the actual coordinates.  Treating that
        // line as a coordinate expands the normalization range and can make
        // the obstacle thinner than a single mesh cell.
        if (points.empty() && x >= 3.0 && y >= 3.0 &&
            std::floor(x) == x && std::floor(y) == y) {
            continue;
        }

        points.push_back({x, y});
    }

    if (points.size() < 3) {
        throw std::runtime_error("Profile file does not contain enough points: " + path.string());
    }

    // Match LBM Geometry::loadSection: Lednicer files store both surfaces
    // leading-edge to trailing-edge. Reverse the lower surface to form one
    // non-self-intersecting polygon loop.
    std::size_t split = 0;
    for (std::size_t i = 1; i < points.size(); ++i) {
        if (points[i].x < points[i - 1].x - 0.5) {
            split = i;
            break;
        }
    }
    if (split != 0) {
        std::vector<Vec2> loop;
        loop.reserve(points.size());
        loop.insert(loop.end(), points.begin(), points.begin() + static_cast<std::ptrdiff_t>(split));
        for (std::size_t i = points.size(); i-- > split;) {
            loop.push_back(points[i]);
        }
        return loop;
    }
    return points;
}

std::vector<Vec2> normalize_profile(std::vector<Vec2> points)
{
    double min_x = points[0].x;
    double max_x = points[0].x;
    for (const auto& point : points) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
    }

    const double chord = max_x - min_x;
    if (std::abs(chord) < 1e-12) {
        throw std::runtime_error("Profile chord length is zero after parsing.");
    }

    for (auto& point : points) {
        point.x = (point.x - min_x) / chord;
        point.y = point.y / chord;
    }

    return points;
}

std::vector<Vec2> make_closed_polygon(const std::vector<Vec2>& profile)
{
    if (profile.size() < 3) {
        throw std::runtime_error("Profile file does not contain enough points to build a polygon.");
    }

    std::vector<Vec2> polygon = profile;

    while (polygon.size() >= 2 && polygon.front().x == polygon.back().x && polygon.front().y == polygon.back().y) {
        polygon.pop_back();
    }

    if (polygon.front().x != polygon.back().x || polygon.front().y != polygon.back().y) {
        polygon.push_back(polygon.front());
    }

    return polygon;
}

std::vector<Vec2> transform_profile(const std::vector<Vec2>& profile, const Config& config)
{
    if (config.domain_width() <= 0.0 || config.domain_height() <= 0.0) {
        throw std::runtime_error("Domain width and height must be positive.");
    }
    if (config.chord <= 0.0 || config.physical_chord <= 0.0) {
        throw std::runtime_error("Reference and physical wing chords must be positive.");
    }
    constexpr double kPi = 3.14159265358979323846;
    // LBM convention: positive AoA rotates the section clockwise about
    // the leading edge for flow in the +x direction.
    const double theta = -config.angle_of_attack_deg * kPi / 180.0;
    const double center_x = config.obstacle_x * config.domain_width();
    const double center_y = config.obstacle_y * config.domain_height();
    const double scaled_chord = config.chord * config.length_scale();

    std::vector<Vec2> transformed;
    transformed.reserve(profile.size());
    for (const auto& point : profile) {
        const double local_x = point.x * scaled_chord;
        const double local_y = point.y * scaled_chord;
        const double rotated_x = local_x * std::cos(theta) - local_y * std::sin(theta);
        const double rotated_y = local_x * std::sin(theta) + local_y * std::cos(theta);
        transformed.push_back({center_x + rotated_x, center_y + rotated_y});
    }

    return transformed;
}

bool point_in_polygon(const Vec2& point, const std::vector<Vec2>& polygon)
{
    bool inside = false;
    const std::size_t count = polygon.size();
    for (std::size_t i = 0, j = count - 1; i < count; j = i++) {
        const Vec2& a = polygon[i];
        const Vec2& b = polygon[j];
        const bool intersects = ((a.y > point.y) != (b.y > point.y)) &&
            (point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y + 1e-12) + a.x);
        if (intersects) {
            inside = !inside;
        }
    }
    return inside;
}

std::vector<std::vector<bool>> build_polygon_mask(const Config& config, const std::vector<Vec2>& polygon)
{
    std::vector<std::vector<bool>> mask(config.ny, std::vector<bool>(config.nx, false));
    const double dx = config.domain_width() / static_cast<double>(config.nx);
    const double dy = config.domain_height() / static_cast<double>(config.ny);

    for (int j = 0; j < config.ny; ++j) {
        for (int i = 0; i < config.nx; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * dx;
            const double y = (static_cast<double>(j) + 0.5) * dy;
            mask[j][i] = point_in_polygon({x, y}, polygon);
        }
    }
    return mask;
}

std::vector<std::vector<bool>> build_circle_mask(const Config& config)
{
    std::vector<std::vector<bool>> mask(config.ny, std::vector<bool>(config.nx, false));
    const double dx_cell = config.domain_width() / static_cast<double>(config.nx);
    const double dy_cell = config.domain_height() / static_cast<double>(config.ny);
    const double center_x = config.obstacle_x * config.domain_width();
    const double center_y = config.obstacle_y * config.domain_height();
    const double radius =
        config.circle_radius_fraction * std::min(config.domain_width(), config.domain_height());

    if (radius <= 0.0 ||
        center_x - radius <= 0.0 ||
        center_x + radius >= config.domain_width() ||
        center_y - radius <= 0.0 ||
        center_y + radius >= config.domain_height()) {
        throw std::runtime_error("Circle must have a positive radius and lie completely inside the domain.");
    }

    for (int j = 0; j < config.ny; ++j) {
        for (int i = 0; i < config.nx; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * dx_cell;
            const double y = (static_cast<double>(j) + 0.5) * dy_cell;
            const double dx = x - center_x;
            const double dy = y - center_y;
            mask[j][i] = (dx * dx + dy * dy) <= (radius * radius);
        }
    }
    return mask;
}

std::vector<std::vector<bool>> build_mask(const Config& config, const std::vector<Vec2>& polygon)
{
    if (config.obstacle_type == "circle") {
        return build_circle_mask(config);
    }
    if (config.obstacle_type == "wing") {
        return build_polygon_mask(config, polygon);
    }
    throw std::runtime_error("Unsupported obstacle type: " + config.obstacle_type + ". Use 'wing' or 'circle'.");
}

void fill_enclosed_fluid_cells(std::vector<std::vector<bool>>& mask)
{
    const int ny = static_cast<int>(mask.size());
    const int nx = ny == 0 ? 0 : static_cast<int>(mask.front().size());
    std::vector<std::vector<bool>> exterior(ny, std::vector<bool>(nx, false));
    std::deque<std::pair<int, int>> queue;

    const auto add_exterior = [&](const int i, const int j) {
        if (i >= 0 && i < nx && j >= 0 && j < ny &&
            !mask[j][i] && !exterior[j][i]) {
            exterior[j][i] = true;
            queue.emplace_back(i, j);
        }
    };

    for (int i = 0; i < nx; ++i) {
        add_exterior(i, 0);
        add_exterior(i, ny - 1);
    }
    for (int j = 0; j < ny; ++j) {
        add_exterior(0, j);
        add_exterior(nx - 1, j);
    }

    constexpr int offsets[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    while (!queue.empty()) {
        const auto [i, j] = queue.front();
        queue.pop_front();
        for (const auto& offset : offsets) {
            add_exterior(i + offset[0], j + offset[1]);
        }
    }

    // Rasterised airfoil outlines can enclose tiny islands of nominal fluid.
    // They are not part of the external flow domain and would create separate
    // OpenFOAM cell regions, so fold them into the obstacle.
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (!mask[j][i] && !exterior[j][i]) {
                mask[j][i] = true;
            }
        }
    }
}

fs::path resolve_profile_path(const std::string& requested)
{
    const fs::path requested_path(requested);
    if (requested_path.is_absolute()) {
        return requested_path;
    }

    const std::vector<fs::path> candidates = {
        fs::current_path() / requested_path,
        fs::current_path().parent_path() / requested_path,
        fs::current_path().parent_path().parent_path() / requested_path,
        fs::current_path().parent_path().parent_path().parent_path() / requested_path,
    };

    for (const auto& candidate : candidates) {
        if (fs::exists(candidate)) {
            return candidate;
        }
    }

    return fs::current_path() / requested_path;
}

fs::path resolve_openfoam_case_root()
{
    fs::path current = fs::current_path();
    while (!current.empty()) {
        if (fs::exists(current / "constant" / "polyMesh") && fs::exists(current / "system")) {
            return current;
        }
        const fs::path parent = current.parent_path();
        if (parent == current) {
            break;
        }
        current = parent;
    }
    return fs::current_path().parent_path();
}

void write_initial_fields(const fs::path& case_root, const Config& config)
{
    const fs::path initial_dir = case_root / "0";
    fs::create_directories(initial_dir);

    const fs::path pressure_path = initial_dir / "p";
    std::ofstream pressure(pressure_path);
    if (!pressure.is_open()) {
        throw std::runtime_error("Unable to write pressure field: " + pressure_path.string());
    }
    pressure << "FoamFile\n"
             << "{\n"
             << "    version     2.0;\n"
             << "    format      ascii;\n"
             << "    class       volScalarField;\n"
             << "    object      p;\n"
             << "}\n\n"
             << "dimensions      [0 2 -2 0 0 0 0];\n"
             << "internalField   uniform 0;\n\n"
             << "boundaryField\n"
             << "{\n"
             << "    inlet\n"
             << "    {\n"
             << "        type            zeroGradient;\n"
             << "    }\n"
             << "    outlet\n"
             << "    {\n"
             << "        type            fixedValue;\n"
             << "        value           uniform 0;\n"
             << "    }\n"
             << "    top\n"
             << "    {\n"
             << "        type            zeroGradient;\n"
             << "    }\n"
             << "    bottom\n"
             << "    {\n"
             << "        type            zeroGradient;\n"
             << "    }\n"
             << "    obstacle\n"
             << "    {\n"
             << "        type            zeroGradient;\n"
             << "    }\n"
             << "    front\n"
             << "    {\n"
             << "        type            empty;\n"
             << "    }\n"
             << "    back\n"
             << "    {\n"
             << "        type            empty;\n"
             << "    }\n"
             << "}\n";

    const fs::path velocity_path = initial_dir / "U";
    std::ofstream velocity(velocity_path);
    if (!velocity.is_open()) {
        throw std::runtime_error("Unable to write velocity field: " + velocity_path.string());
    }
    velocity << "FoamFile\n"
             << "{\n"
             << "    version     2.0;\n"
             << "    format      ascii;\n"
             << "    class       volVectorField;\n"
             << "    object      U;\n"
             << "}\n\n"
             << "dimensions      [0 1 -1 0 0 0 0];\n"
             << "internalField   uniform (" << config.airflow_x << " " << config.airflow_y << " 0);\n\n"
             << "boundaryField\n"
             << "{\n"
             << "    inlet\n"
             << "    {\n"
             << "        type            fixedValue;\n"
             << "        value           uniform (" << config.inlet_velocity << " 0 0);\n"
             << "    }\n"
             << "    outlet\n"
             << "    {\n"
             << "        type            zeroGradient;\n"
             << "    }\n"
             << "    top\n"
             << "    {\n"
             << "        type            slip;\n"
             << "    }\n"
             << "    bottom\n"
             << "    {\n"
             << "        type            slip;\n"
             << "    }\n"
             << "    obstacle\n"
             << "    {\n"
             << "        type            noSlip;\n"
             << "    }\n"
             << "    front\n"
             << "    {\n"
             << "        type            empty;\n"
             << "    }\n"
             << "    back\n"
             << "    {\n"
             << "        type            empty;\n"
             << "    }\n"
             << "}\n";
}

void verify_wing_bounds(const std::vector<Vec2>& polygon, const Config& config)
{
    double min_x = polygon.front().x;
    double max_x = polygon.front().x;
    double min_y = polygon.front().y;
    double max_y = polygon.front().y;
    for (const auto& point : polygon) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
    }

    if (min_x <= 0.0 || max_x >= config.domain_width() ||
        min_y <= 0.0 || max_y >= config.domain_height()) {
        throw std::runtime_error(
            "Wing geometry leaves the computational domain. Expected x/y inside [0, " +
            std::to_string(config.domain_width()) + "] and [0, " +
            std::to_string(config.domain_height()) + "].");
    }
}

std::string make_face_key(const std::vector<int>& points)
{
    std::vector<int> sorted = points;
    std::sort(sorted.begin(), sorted.end());
    std::ostringstream stream;
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        if (i > 0) {
            stream << ',';
        }
        stream << sorted[i];
    }
    return stream.str();
}

void write_poly_mesh(const std::vector<std::vector<bool>>& mask, const Config& config)
{
    const fs::path case_root = resolve_openfoam_case_root();
    const fs::path poly_mesh_dir = case_root / "constant" / "polyMesh";
    fs::create_directories(poly_mesh_dir);

    // OpenFOAM finite-volume cells must have non-zero volume.  Extrude the
    // two-dimensional grid by one cell in z and make its two sides empty
    // patches.  Points are allocated lazily so vertices wholly inside the
    // obstacle do not become unused mesh points.
    std::vector<Vec3> points;
    std::map<std::tuple<int, int, int>, int> point_indices;
    const double dx = config.domain_width() / static_cast<double>(config.nx);
    const double dy = config.domain_height() / static_cast<double>(config.ny);
    const double thickness = std::min(dx, dy);
    auto point_index = [&](const int i, const int j, const int k) {
        const auto key = std::make_tuple(i, j, k);
        const auto existing = point_indices.find(key);
        if (existing != point_indices.end()) {
            return existing->second;
        }
        const int index = static_cast<int>(points.size());
        points.push_back({
            static_cast<double>(i) * dx,
            static_cast<double>(j) * dy,
            static_cast<double>(k) * thickness});
        point_indices.emplace(key, index);
        return index;
    };
    auto face_vertices = [&](const int i, const int j, const int direction) {
        switch (direction) {
        case 0: // west, outward -x
            return std::vector<int>{point_index(i, j, 0), point_index(i, j, 1),
                point_index(i, j + 1, 1), point_index(i, j + 1, 0)};
        case 1: // east, outward +x
            return std::vector<int>{point_index(i + 1, j, 0), point_index(i + 1, j + 1, 0),
                point_index(i + 1, j + 1, 1), point_index(i + 1, j, 1)};
        case 2: // south, outward -y
            return std::vector<int>{point_index(i, j, 0), point_index(i + 1, j, 0),
                point_index(i + 1, j, 1), point_index(i, j, 1)};
        case 3: // north, outward +y
            return std::vector<int>{point_index(i, j + 1, 0), point_index(i, j + 1, 1),
                point_index(i + 1, j + 1, 1), point_index(i + 1, j + 1, 0)};
        case 4: // front, outward -z
            return std::vector<int>{point_index(i, j, 0), point_index(i, j + 1, 0),
                point_index(i + 1, j + 1, 0), point_index(i + 1, j, 0)};
        case 5: // back, outward +z
            return std::vector<int>{point_index(i, j, 1), point_index(i + 1, j, 1),
                point_index(i + 1, j + 1, 1), point_index(i, j + 1, 1)};
        default:
            throw std::runtime_error("Unsupported face direction.");
        }
    };

    std::vector<FaceRecord> faces;
    std::map<std::string, int> face_index_by_key;
    const std::vector<std::string> patch_order = {
        "inlet", "outlet", "bottom", "top", "obstacle", "front", "back"};

    int cell_id = 0;
    for (int j = 0; j < config.ny; ++j) {
        for (int i = 0; i < config.nx; ++i) {
            if (mask[j][i]) {
                continue;
            }

            const int cell = cell_id++;
            for (int direction = 0; direction < 6; ++direction) {
                int neighbor_i = i;
                int neighbor_j = j;
                std::string patch_name;

                switch (direction) {
                case 0: // west
                    if (i == 0) {
                        patch_name = "inlet";
                    } else {
                        neighbor_i = i - 1;
                    }
                    break;
                case 1: // east
                    if (i + 1 >= config.nx) {
                        patch_name = "outlet";
                    } else {
                        neighbor_i = i + 1;
                    }
                    break;
                case 2: // south
                    if (j == 0) {
                        patch_name = "bottom";
                    } else {
                        neighbor_j = j - 1;
                    }
                    break;
                case 3: // north
                    if (j + 1 >= config.ny) {
                        patch_name = "top";
                    } else {
                        neighbor_j = j + 1;
                    }
                    break;
                case 4:
                    patch_name = "front";
                    break;
                case 5:
                    patch_name = "back";
                    break;
                }

                if (direction < 4 && patch_name.empty()) {
                    const bool neighbor_solid = mask[neighbor_j][neighbor_i];
                    if (neighbor_solid) {
                        patch_name = "obstacle";
                    }
                }

                const std::vector<int> face_points = face_vertices(i, j, direction);
                const std::string key = make_face_key(face_points);
                const auto existing = face_index_by_key.find(key);

                if (existing == face_index_by_key.end()) {
                    FaceRecord record;
                    record.points = face_points;
                    record.owner = cell;
                    record.patch_name = patch_name;
                    faces.push_back(record);
                    face_index_by_key[key] = static_cast<int>(faces.size() - 1);
                } else {
                    faces[existing->second].neighbour = cell;
                    faces[existing->second].patch_name = "";
                }
            }
        }
    }

    std::vector<FaceRecord> ordered_faces;
    ordered_faces.reserve(faces.size());
    std::vector<int> owner_list;
    std::vector<int> neighbour_list;

    for (const auto& face : faces) {
        if (face.neighbour >= 0) {
            ordered_faces.push_back(face);
        }
    }

    // OpenFOAM requires every patch to occupy one contiguous range in the
    // face list.  Append boundary faces patch-by-patch, in the same order
    // used below when writing startFace.
    std::map<std::string, std::size_t> patch_face_counts;
    for (const auto& patch_name : patch_order) {
        for (const auto& face : faces) {
            if (face.neighbour < 0 && face.patch_name == patch_name) {
                ordered_faces.push_back(face);
                ++patch_face_counts[patch_name];
            }
        }
    }

    const auto obstacle_it = patch_face_counts.find("obstacle");
    if (obstacle_it == patch_face_counts.end() || obstacle_it->second == 0) {
        throw std::runtime_error(
            "The obstacle does not cover any mesh cells. Increase its size/resolution "
            "or check the airfoil profile format.");
    }

    owner_list.reserve(ordered_faces.size());
    neighbour_list.reserve(ordered_faces.size());
    for (const auto& face : ordered_faces) {
        owner_list.push_back(face.owner);
        if (face.neighbour >= 0) {
            neighbour_list.push_back(face.neighbour);
        }
    }

    const std::size_t n_internal_faces = neighbour_list.size();
    const std::size_t n_faces = ordered_faces.size();

    const fs::path points_path = poly_mesh_dir / "points";
    std::ofstream points_stream(points_path);
    if (!points_stream.is_open()) {
        throw std::runtime_error("Unable to write points file: " + points_path.string());
    }
    points_stream << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    points_stream << "  =========                 |\n";
    points_stream << "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    points_stream << "   \\    /   O peration     | Website: https://openfoam.org\n";
    points_stream << "    \\  /    A nd           | Version:  dev\n";
    points_stream << "     \\/     M anipulation  |\n";
    points_stream << "\\*---------------------------------------------------------------------------*/\n";
    points_stream << "FoamFile\n";
    points_stream << "{\n";
    points_stream << "    format      ascii;\n";
    points_stream << "    class       vectorField;\n";
    points_stream << "    location    \"constant/polyMesh\";\n";
    points_stream << "    object      points;\n";
    points_stream << "}\n";
    points_stream << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    points_stream << points.size() << "\n(\n";
    for (const auto& point : points) {
        points_stream << "(" << point.x << " " << point.y << " " << point.z << ")\n";
    }
    points_stream << ")\n";

    const fs::path owner_path = poly_mesh_dir / "owner";
    std::ofstream owner_stream(owner_path);
    if (!owner_stream.is_open()) {
        throw std::runtime_error("Unable to write owner file: " + owner_path.string());
    }
    owner_stream << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    owner_stream << "  =========                 |\n";
    owner_stream << "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    owner_stream << "   \\    /   O peration     | Website: https://openfoam.org\n";
    owner_stream << "    \\  /    A nd           | Version:  dev\n";
    owner_stream << "     \\/     M anipulation  |\n";
    owner_stream << "\\*---------------------------------------------------------------------------*/\n";
    owner_stream << "FoamFile\n";
    owner_stream << "{\n";
    owner_stream << "    format      ascii;\n";
    owner_stream << "    class       labelList;\n";
    owner_stream << "    location    \"constant/polyMesh\";\n";
    owner_stream << "    object      owner;\n";
    owner_stream << "}\n";
    owner_stream << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    owner_stream << n_faces << "\n(\n";
    for (const auto& label : owner_list) {
        owner_stream << label << "\n";
    }
    owner_stream << ")\n";

    const fs::path neighbour_path = poly_mesh_dir / "neighbour";
    std::ofstream neighbour_stream(neighbour_path);
    if (!neighbour_stream.is_open()) {
        throw std::runtime_error("Unable to write neighbour file: " + neighbour_path.string());
    }
    neighbour_stream << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    neighbour_stream << "  =========                 |\n";
    neighbour_stream << "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    neighbour_stream << "   \\    /   O peration     | Website: https://openfoam.org\n";
    neighbour_stream << "    \\  /    A nd           | Version:  dev\n";
    neighbour_stream << "     \\/     M anipulation  |\n";
    neighbour_stream << "\\*---------------------------------------------------------------------------*/\n";
    neighbour_stream << "FoamFile\n";
    neighbour_stream << "{\n";
    neighbour_stream << "    format      ascii;\n";
    neighbour_stream << "    class       labelList;\n";
    neighbour_stream << "    location    \"constant/polyMesh\";\n";
    neighbour_stream << "    object      neighbour;\n";
    neighbour_stream << "}\n";
    neighbour_stream << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    neighbour_stream << n_internal_faces << "\n(\n";
    for (const auto& label : neighbour_list) {
        neighbour_stream << label << "\n";
    }
    neighbour_stream << ")\n";

    const fs::path faces_path = poly_mesh_dir / "faces";
    std::ofstream faces_stream(faces_path);
    if (!faces_stream.is_open()) {
        throw std::runtime_error("Unable to write faces file: " + faces_path.string());
    }
    faces_stream << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    faces_stream << "  =========                 |\n";
    faces_stream << "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    faces_stream << "   \\    /   O peration     | Website: https://openfoam.org\n";
    faces_stream << "    \\  /    A nd           | Version:  dev\n";
    faces_stream << "     \\/     M anipulation  |\n";
    faces_stream << "\\*---------------------------------------------------------------------------*/\n";
    faces_stream << "FoamFile\n";
    faces_stream << "{\n";
    faces_stream << "    format      ascii;\n";
    faces_stream << "    class       faceList;\n";
    faces_stream << "    location    \"constant/polyMesh\";\n";
    faces_stream << "    object      faces;\n";
    faces_stream << "}\n";
    faces_stream << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    faces_stream << n_faces << "\n(\n";
    for (const auto& face : ordered_faces) {
        faces_stream << face.points.size() << "(";
        for (std::size_t i = 0; i < face.points.size(); ++i) {
            if (i > 0) {
                faces_stream << " ";
            }
            faces_stream << face.points[i];
        }
        faces_stream << ")\n";
    }
    faces_stream << ")\n";

    const fs::path boundary_path = poly_mesh_dir / "boundary";
    std::ofstream boundary_stream(boundary_path);
    if (!boundary_stream.is_open()) {
        throw std::runtime_error("Unable to write boundary file: " + boundary_path.string());
    }
    boundary_stream << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    boundary_stream << "  =========                 |\n";
    boundary_stream << "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    boundary_stream << "   \\    /   O peration     | Website: https://openfoam.org\n";
    boundary_stream << "    \\  /    A nd           | Version:  dev\n";
    boundary_stream << "     \\/     M anipulation  |\n";
    boundary_stream << "\\*---------------------------------------------------------------------------*/\n";
    boundary_stream << "FoamFile\n";
    boundary_stream << "{\n";
    boundary_stream << "    format      ascii;\n";
    boundary_stream << "    class       polyBoundaryMesh;\n";
    boundary_stream << "    location    \"constant/polyMesh\";\n";
    boundary_stream << "    object      boundary;\n";
    boundary_stream << "}\n";
    boundary_stream << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    std::size_t active_patch_count = 0;
    for (const auto& patch_name : patch_order) {
        const auto patch_it = patch_face_counts.find(patch_name);
        if (patch_it != patch_face_counts.end() && patch_it->second > 0) {
            ++active_patch_count;
        }
    }
    boundary_stream << active_patch_count << "\n(\n";

    std::size_t boundary_offset = n_internal_faces;
    for (const auto& patch_name : patch_order) {
        const auto patch_it = patch_face_counts.find(patch_name);
        if (patch_it == patch_face_counts.end() || patch_it->second == 0) {
            continue;
        }
        boundary_stream << "    " << patch_name << "\n";
        boundary_stream << "    {\n";
        boundary_stream << "        type            ";
        if (patch_name == "front" || patch_name == "back") {
            boundary_stream << "empty;\n";
        } else if (patch_name == "top" || patch_name == "bottom" || patch_name == "obstacle") {
            boundary_stream << "wall;\n";
        } else {
            boundary_stream << "patch;\n";
        }
        boundary_stream << "        nFaces          " << patch_it->second << ";\n";
        boundary_stream << "        startFace       " << boundary_offset << ";\n";
        boundary_stream << "    }\n";
        boundary_offset += patch_it->second;
    }
    boundary_stream << ")\n\n// ************************************************************************* //\n";

    write_initial_fields(case_root, config);
    std::cout << "Generated OpenFOAM polyMesh files in " << poly_mesh_dir.string() << std::endl;
}

} // namespace

void generate_mesh(const Config& config)
{
    Config effective = config;
    if (effective.resolution < 2) {
        throw std::runtime_error("resolution must be at least 2.");
    }
    effective.nx = effective.resolution - 1;
    effective.ny = effective.resolution - 1;

    const fs::path profile_path = resolve_profile_path(effective.wing_profile_file);
    std::vector<Vec2> profile = normalize_profile(load_profile(profile_path));
    std::vector<Vec2> closed_profile = make_closed_polygon(profile);

    std::vector<Vec2> transformed_profile = transform_profile(closed_profile, effective);
    verify_wing_bounds(transformed_profile, effective);

    std::vector<std::vector<bool>> mask = build_mask(effective, transformed_profile);
    fill_enclosed_fluid_cells(mask);
    write_poly_mesh(mask, effective);
}

} // namespace meshgen
