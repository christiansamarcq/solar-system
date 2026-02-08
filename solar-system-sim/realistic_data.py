"""
Realistic solar system data from NASA
Data sources:
- NASA Planetary Fact Sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/
- JPL Solar System Dynamics: https://ssd.jpl.nasa.gov/planets/phys_par.html
- All values are in SI units unless otherwise specified
"""

import numpy as np
import math

# Scale factors for visualization
# Real solar system is too large to visualize at 1:1 scale
DISTANCE_SCALE = 1.0 / 1e9  # 1 unit = 1 billion meters (1 Gm)
RADIUS_SCALE = 1.0 / 1e9  # 1 unit = 1 billion meters (REALISTIC SCALE - planets will be tiny!)
TIME_SCALE = 1.0  # 1 simulation second = 1 real second (will be adjusted by user)

# Gravitational constant
G = 6.67430e-11  # m^3 kg^-1 s^-2

# Sun data
SUN = {
    'name': 'Sun',
    'mass': 1.989e30,  # kg
    'radius': 696000 * 1e3,  # meters (696,000 km)
    'color': (1.0, 0.9, 0.6),
    'texture': 'textures/sun.jpg',
    'is_emissive': True,
    'luminosity': 3.828e26  # watts
}

# Planet data from NASA Planetary Fact Sheet
# Format: mass (kg), radius (m), semi-major axis (m), eccentricity, inclination (deg), orbital period (days)
PLANETS = {
    'Mercury': {
        'name': 'Mercury',
        'mass': 0.33011e24,  # kg
        'radius': 2439.7 * 1e3,  # meters
        'semi_major_axis': 57.909e9,  # meters (0.387 AU)
        'eccentricity': 0.2056,
        'inclination': math.radians(7.005),  # degrees to radians
        'orbital_period': 87.97,  # days
        'rotation_period': 58.646,  # days
        'arg_periapsis': math.radians(29.124),  # argument of perihelion
        'color': (0.7, 0.7, 0.7),
        'texture': 'textures/mercury.jpg'
    },
    'Venus': {
        'name': 'Venus',
        'mass': 4.8675e24,  # kg
        'radius': 6051.8 * 1e3,  # meters
        'semi_major_axis': 108.21e9,  # meters (0.723 AU)
        'eccentricity': 0.0067,
        'inclination': math.radians(3.39),
        'orbital_period': 224.70,  # days
        'rotation_period': -243.025,  # days (retrograde, hence negative)
        'arg_periapsis': math.radians(54.884),
        'color': (0.9, 0.8, 0.6),
        'texture': 'textures/venus.jpg'
    },
    'Earth': {
        'name': 'Earth',
        'mass': 5.9724e24,  # kg
        'radius': 6371.0 * 1e3,  # meters
        'semi_major_axis': 149.60e9,  # meters (1.000 AU)
        'eccentricity': 0.0167,
        'inclination': math.radians(0.0),  # reference plane
        'orbital_period': 365.26,  # days
        'rotation_period': 0.99726968,  # days (23.9345 hours)
        'arg_periapsis': math.radians(114.20783),
        'color': (0.2, 0.4, 0.8),
        'texture': 'textures/earth.jpg'
    },
    'Mars': {
        'name': 'Mars',
        'mass': 0.64171e24,  # kg
        'radius': 3389.5 * 1e3,  # meters
        'semi_major_axis': 227.92e9,  # meters (1.524 AU)
        'eccentricity': 0.0934,
        'inclination': math.radians(1.85),
        'orbital_period': 686.98,  # days
        'rotation_period': 1.025957,  # days
        'arg_periapsis': math.radians(286.5),
        'color': (0.8, 0.4, 0.2),
        'texture': 'textures/mars.jpg'
    },
    'Jupiter': {
        'name': 'Jupiter',
        'mass': 1898.19e24,  # kg
        'radius': 69911 * 1e3,  # meters (equatorial)
        'semi_major_axis': 778.57e9,  # meters (5.204 AU)
        'eccentricity': 0.0489,
        'inclination': math.radians(1.304),
        'orbital_period': 4332.59,  # days (~11.86 years)
        'rotation_period': 0.41354,  # days (~9.925 hours)
        'arg_periapsis': math.radians(273.867),
        'color': (0.8, 0.7, 0.5),
        'texture': 'textures/jupiter.jpg'
    },
    'Saturn': {
        'name': 'Saturn',
        'mass': 568.34e24,  # kg
        'radius': 58232 * 1e3,  # meters (equatorial, excluding rings)
        'semi_major_axis': 1433.53e9,  # meters (9.582 AU)
        'eccentricity': 0.0565,
        'inclination': math.radians(2.485),
        'orbital_period': 10759.22,  # days (~29.46 years)
        'rotation_period': 0.44401,  # days (~10.656 hours)
        'arg_periapsis': math.radians(339.392),
        'color': (0.9, 0.8, 0.6),
        'texture': 'textures/saturn.jpg'
    },
    'Uranus': {
        'name': 'Uranus',
        'mass': 86.813e24,  # kg
        'radius': 25362 * 1e3,  # meters (equatorial)
        'semi_major_axis': 2872.46e9,  # meters (19.19 AU)
        'eccentricity': 0.0457,
        'inclination': math.radians(0.773),
        'orbital_period': 30688.5,  # days (~84.01 years)
        'rotation_period': -0.71833,  # days (retrograde, ~17.24 hours)
        'arg_periapsis': math.radians(96.998857),
        'color': (0.5, 0.7, 0.8),
        'texture': 'textures/uranus.jpg'
    },
    'Neptune': {
        'name': 'Neptune',
        'mass': 102.413e24,  # kg
        'radius': 24622 * 1e3,  # meters (equatorial)
        'semi_major_axis': 4495.06e9,  # meters (30.07 AU)
        'eccentricity': 0.0113,
        'inclination': math.radians(1.767),
        'orbital_period': 60182,  # days (~164.79 years)
        'rotation_period': 0.67125,  # days (~16.11 hours)
        'arg_periapsis': math.radians(273.187),
        'color': (0.3, 0.4, 0.9),
        'texture': 'textures/neptune.jpg'
    },
    'Pluto': {
        'name': 'Pluto',
        'mass': 0.01303e24,  # kg (dwarf planet)
        'radius': 1188.3 * 1e3,  # meters
        'semi_major_axis': 5906.4e9,  # meters (39.48 AU)
        'eccentricity': 0.2488,  # highly eccentric orbit
        'inclination': math.radians(17.16),  # highly inclined
        'orbital_period': 90560,  # days (~248 years)
        'rotation_period': -6.387,  # days (retrograde, ~6.4 days)
        'arg_periapsis': math.radians(113.834),
        'color': (0.8, 0.7, 0.6),  # brownish color
        'texture': 'textures/moon.jpg'  # using generic moon texture
    }
}

# Major moons with orbital data
MOONS = {
    'Moon': {  # Earth's Moon
        'name': 'Moon',
        'parent': 'Earth',
        'mass': 0.07346e24,  # kg
        'radius': 1737.4 * 1e3,  # meters
        'semi_major_axis': 384.4e6,  # meters (384,400 km from Earth)
        'eccentricity': 0.0549,
        'inclination': math.radians(5.145),  # to ecliptic
        'orbital_period': 27.3217,  # days
        'color': (0.7, 0.7, 0.7),
        'texture': 'textures/moon.jpg'
    },
    'Io': {  # Jupiter's moon
        'name': 'Io',
        'parent': 'Jupiter',
        'mass': 0.08932e24,  # kg
        'radius': 1821.6 * 1e3,  # meters
        'semi_major_axis': 421.8e6,  # meters
        'eccentricity': 0.0041,
        'inclination': math.radians(0.05),
        'orbital_period': 1.769,  # days
        'color': (0.9, 0.8, 0.3),
        'texture': 'textures/moon.jpg'  # using generic moon texture
    },
    'Europa': {  # Jupiter's moon
        'name': 'Europa',
        'parent': 'Jupiter',
        'mass': 0.04800e24,  # kg
        'radius': 1560.8 * 1e3,  # meters
        'semi_major_axis': 671.1e6,  # meters
        'eccentricity': 0.0094,
        'inclination': math.radians(0.47),
        'orbital_period': 3.551,  # days
        'color': (0.8, 0.8, 0.7),
        'texture': 'textures/moon.jpg'
    },
    'Ganymede': {  # Jupiter's moon (largest moon in solar system)
        'name': 'Ganymede',
        'parent': 'Jupiter',
        'mass': 0.1482e24,  # kg
        'radius': 2634.1 * 1e3,  # meters
        'semi_major_axis': 1070.4e6,  # meters
        'eccentricity': 0.0013,
        'inclination': math.radians(0.20),
        'orbital_period': 7.155,  # days
        'color': (0.6, 0.6, 0.5),
        'texture': 'textures/moon.jpg'
    },
    'Callisto': {  # Jupiter's moon
        'name': 'Callisto',
        'parent': 'Jupiter',
        'mass': 0.1076e24,  # kg
        'radius': 2410.3 * 1e3,  # meters
        'semi_major_axis': 1882.7e6,  # meters
        'eccentricity': 0.0074,
        'inclination': math.radians(0.51),
        'orbital_period': 16.689,  # days
        'color': (0.5, 0.5, 0.4),
        'texture': 'textures/moon.jpg'
    },
    'Titan': {  # Saturn's moon (2nd largest in solar system)
        'name': 'Titan',
        'parent': 'Saturn',
        'mass': 0.13452e24,  # kg
        'radius': 2574.7 * 1e3,  # meters
        'semi_major_axis': 1221.87e6,  # meters
        'eccentricity': 0.0288,
        'inclination': math.radians(0.34854),
        'orbital_period': 15.945,  # days
        'color': (0.8, 0.6, 0.4),
        'texture': 'textures/moon.jpg'
    }
}


def get_scaled_radius(radius_m):
    """Convert real radius in meters to scaled visualization radius"""
    return radius_m * RADIUS_SCALE


def get_scaled_distance(distance_m):
    """Convert real distance in meters to scaled visualization distance"""
    return distance_m * DISTANCE_SCALE


def get_orbital_period_seconds(period_days):
    """Convert orbital period from days to seconds"""
    return period_days * 24.0 * 3600.0


def calculate_orbital_velocity(semi_major_axis, central_mass):
    """
    Calculate approximate circular orbital velocity
    v = sqrt(G * M / r)

    Args:
        semi_major_axis: Distance from central body (meters)
        central_mass: Mass of central body (kg)

    Returns:
        Orbital velocity in m/s
    """
    return math.sqrt(G * central_mass / semi_major_axis)
