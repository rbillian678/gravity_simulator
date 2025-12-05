#include <iostream>
#include <vector>
#include <cstring> // for strncpy
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>

struct ThreeDVect
{
    double x, y, z;

    ThreeDVect() : x(0.0f), y(0.0f), z(0.0f) {}

    ThreeDVect(double x_, double y_, double z_) : x{x_}, y{y_}, z{z_} {}
};

class Planet
{
public:
    char name[10];
    ThreeDVect position;
    ThreeDVect velocity;
    ThreeDVect acceleration;
    double mass; // kg

    Planet() : name{0}, position(), velocity(), acceleration(), mass{0.0f}
    {
        name[0] = '\0';
    }

    Planet(const char *nm, ThreeDVect pos, ThreeDVect vel, ThreeDVect acc, double m)
        : position{pos}, velocity{vel}, acceleration{acc}, mass{m}
    {
        std::strncpy(this->name, nm, sizeof(this->name));
        this->name[sizeof(this->name) - 1] = '\0';
    }
};

inline void updateForce(ThreeDVect &force, const Planet &myself, const Planet &other)
{
    constexpr double G = 6.67430e-11f;
    double numerator = G * other.mass * myself.mass;
    double r_x = other.position.x - myself.position.x;
    double r_y = other.position.y - myself.position.y;
    double r_z = other.position.z - myself.position.z;
    double denom = r_x * r_x + r_y * r_y + r_z * r_z;
    force.x += numerator / (denom * sqrtf(denom)) * r_x;

    force.y += numerator / (denom * sqrtf(denom)) * r_y;

    force.z += numerator / (denom * sqrtf(denom)) * r_z;
}

inline void updatePosition(Planet &planet, double dt)
{
    planet.position.x += dt * planet.velocity.x;
    planet.position.y += dt * planet.velocity.y;
    planet.position.z += dt * planet.velocity.z;
}

inline void updateVelocity(Planet &planet, double dt)
{
    planet.velocity.x += dt * planet.acceleration.x;
    planet.velocity.y += dt * planet.acceleration.y;
    planet.velocity.z += dt * planet.acceleration.z;
}

inline void updateAcceleration(Planet &planet, const ThreeDVect &force)
{
    planet.acceleration.x = force.x / planet.mass;
    planet.acceleration.y = force.y / planet.mass;
    planet.acceleration.z = force.z / planet.mass;
}

void computePlanetData(Planet *planets, int n, double dt, int idx)
{
    if (idx < n)
    {
        ThreeDVect force(0.0f, 0.0f, 0.0f);

        Planet &me = planets[idx];

        for (int i = 0; i < n; i++)
        {
            if (i != idx)
            {
                updateForce(force, me, planets[i]);
            }
        }

        updatePosition(me, dt);
        updateVelocity(me, dt);
        updateAcceleration(me, force);
    }
}

void printUpdates(std::vector<Planet> &planets, int n, int day)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << "day: " << day
                  << ", name: " << planets[i].name
                  << ", x: " << planets[i].position.x
                  << ", y: " << planets[i].position.y
                  << ", z: " << planets[i].position.z << std::endl;
    }
}

std::vector<Planet> loadPlanetsFromCSV(const std::string &filename)
{
    std::vector<Planet> planets;

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "ERROR: Cannot open file " << filename << "\n";
        return planets;
    }

    std::string line;

    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        std::stringstream ss(line);
        std::string token;

        std::vector<std::string> cols;
        while (std::getline(ss, token, ','))
        {
            cols.push_back(token);
        }

        if (cols.size() < 11)
        {
            std::cerr << "WARNING: Bad CSV row (needs ≥11 columns): " << line << "\n";
            continue;
        }

        std::string name = cols[0];
        double x = std::stod(cols[1]);
        double y = std::stod(cols[2]);
        double z = std::stod(cols[3]);
        double vx = std::stod(cols[4]);
        double vy = std::stod(cols[5]);
        double vz = std::stod(cols[6]);
        double ax = std::stod(cols[7]);
        double ay = std::stod(cols[8]);
        double az = std::stod(cols[9]);
        double mass = std::stod(cols[10]);

        Planet p(
            name.c_str(),
            ThreeDVect(x, y, z),
            ThreeDVect(vx, vy, vz),
            ThreeDVect(ax, ay, az),
            mass);

        planets.push_back(p);
    }

    return planets;
}

int main()
{

    std::vector<Planet> planets;
    planets = loadPlanetsFromCSV("10_planets.csv");

    double earthDistanceSun = 1.496e11f; // meters
    constexpr int oneYearInSeconds = 365 * 24 * 60 * 60;
    constexpr int oneDayInSeconds = oneYearInSeconds / 365;
    constexpr double PI = 3.14159265358979323846f;
    double earthInitialSpeed = 2.0f * PI * earthDistanceSun / oneYearInSeconds; // m/s
    constexpr int dt = 100;                                                     // time step in seconds
    constexpr int N = oneYearInSeconds / dt;                                    // number of time steps

    for (int i = 0; i < N; i++)
    {
        for (int idx = 0; idx < planets.size(); idx++)
        {
            computePlanetData(planets.data(), planets.size(), dt, idx);
        }
        if (i * dt % oneDayInSeconds == 0)
        {
            printUpdates(planets, planets.size(), i * dt / oneDayInSeconds);
        }
    }

    double dx = planets[1].position.x - earthDistanceSun;
    double dy = planets[1].position.y;
    double dz = planets[1].position.z;
    double dist_err = std::sqrt(dx * dx + dy * dy + dz * dz);

    std::cout << "err: " << dist_err << std::endl;

    return 0;
}