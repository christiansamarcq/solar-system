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
    'luminosity': 3.828e26,  # watts
    'rotation_period': 25.38,  # days (sidereal, equatorial)
    'axial_tilt': 7.25  # degrees (to ecliptic)
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
        'axial_tilt': 0.034,  # degrees
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
        'axial_tilt': 177.4,  # degrees (nearly upside down, retrograde)
        'arg_periapsis': math.radians(54.884),
        'color': (0.9, 0.8, 0.6),
        'texture': 'textures/venus_clouds.jpg',  # cloud tops; venus.jpg is the hidden radar surface
        'brightness': 1.05  # pale clouds: the default 2x boost would wash them out to plain white
    },
    'Earth': {
        'name': 'Earth',
        'mass': 5.9724e24,  # kg
        'radius': 6371.0 * 1e3,  # meters
        'semi_major_axis': 149.60e9,  # meters (1.000 AU)
        'eccentricity': 0.0167,
        'inclination': math.radians(0.0),  # reference plane
        'orbital_period': 365.26,  # days
        'rotation_period': 0.99726968,  # days (23.9355 hours)
        'axial_tilt': 23.44,  # degrees
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
        'axial_tilt': 25.19,  # degrees
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
        'axial_tilt': 3.13,  # degrees
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
        'axial_tilt': 26.73,  # degrees
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
        'axial_tilt': 97.77,  # degrees (rotates on its side)
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
        'axial_tilt': 28.32,  # degrees
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
        'axial_tilt': 122.53,  # degrees (tilted past sideways)
        'arg_periapsis': math.radians(113.834),
        'color': (0.8, 0.7, 0.6),  # brownish color
        'texture': 'textures/moon.jpg'  # using generic moon texture
    },
    'Ceres': {
        'name': 'Ceres',
        'mass': 9.3835e20,  # kg (dwarf planet, largest asteroid belt object)
        'radius': 473.0 * 1e3,  # meters (mean radius)
        'semi_major_axis': 413.69e9,  # meters (2.766 AU)
        'eccentricity': 0.0758,
        'inclination': math.radians(10.594),
        'orbital_period': 1681.63,  # days (~4.60 years)
        'rotation_period': 0.3781,  # days (~9.074 hours)
        'axial_tilt': 4.0,  # degrees
        'arg_periapsis': math.radians(73.597),
        'color': (0.6, 0.6, 0.55),  # grayish
        'texture': 'textures/moon.jpg'
    },
    'Eris': {
        'name': 'Eris',
        'mass': 1.6466e22,  # kg (most massive dwarf planet)
        'radius': 1163.0 * 1e3,  # meters
        'semi_major_axis': 10125.0e9,  # meters (67.67 AU)
        'eccentricity': 0.4407,  # highly eccentric
        'inclination': math.radians(44.04),  # highly inclined
        'orbital_period': 203830,  # days (~558 years)
        'rotation_period': 1.08,  # days (~25.9 hours)
        'axial_tilt': 78.0,  # degrees (estimated)
        'arg_periapsis': math.radians(151.639),
        'color': (0.85, 0.85, 0.85),  # whitish gray
        'texture': 'textures/moon.jpg'
    },
    'Haumea': {
        'name': 'Haumea',
        'mass': 4.006e21,  # kg
        'radius': 816.0 * 1e3,  # meters (mean, elongated shape)
        'semi_major_axis': 6452.0e9,  # meters (43.13 AU)
        'eccentricity': 0.1912,
        'inclination': math.radians(28.22),
        'orbital_period': 103774,  # days (~284 years)
        'rotation_period': 0.1631,  # days (~3.9 hours, fastest spinning large body)
        'axial_tilt': 126.0,  # degrees (estimated)
        'arg_periapsis': math.radians(239.041),
        'color': (0.75, 0.7, 0.65),  # slightly reddish gray
        'texture': 'textures/moon.jpg'
    },
    'Makemake': {
        'name': 'Makemake',
        'mass': 3.1e21,  # kg
        'radius': 715.0 * 1e3,  # meters (mean radius)
        'semi_major_axis': 6796.0e9,  # meters (45.43 AU)
        'eccentricity': 0.1559,
        'inclination': math.radians(28.98),
        'orbital_period': 111845,  # days (~306 years)
        'rotation_period': 0.9511,  # days (~22.8 hours)
        'axial_tilt': 0.0,  # degrees (unknown, assumed 0)
        'arg_periapsis': math.radians(294.834),
        'color': (0.7, 0.55, 0.45),  # reddish brown
        'texture': 'textures/moon.jpg'
    },
    # Dwarf planet CANDIDATE, announced Oct 27, 2026 (found in Vera C. Rubin
    # Observatory data). Orbit is measured; size, mass and spin are estimates.
    '2026 ED138': {
        'name': '2026 ED138',
        'mass': 4.0e20,  # kg (ESTIMATE: ~800 km wide, assumed icy density 1.5 g/cm^3)
        'radius': 400.0 * 1e3,  # meters (ESTIMATE: diameter likely 610-970 km)
        'semi_major_axis': 6446.4e9,  # meters (43.091 AU)
        'eccentricity': 0.1387,
        'inclination': math.radians(24.434),
        'orbital_period': 103250,  # days (~283 years)
        'rotation_period': 1.18077,  # days (~28.3 hours, tentative)
        'axial_tilt': 0.0,  # degrees (unknown, assumed 0)
        'arg_periapsis': math.radians(168.703),
        'color': (0.72, 0.58, 0.5),  # moderately red (icy, tholin-darkened)
        'texture': 'textures/moon.jpg'
    }
}

# Major moons with orbital data
MOONS = {
    'Moon': {  # Earth's Moon (tidally locked)
        'name': 'Moon',
        'parent': 'Earth',
        'mass': 0.07346e24,  # kg
        'radius': 1737.4 * 1e3,  # meters
        'semi_major_axis': 384.4e6,  # meters (384,400 km from Earth)
        'eccentricity': 0.0549,
        'inclination': math.radians(5.145),  # to ecliptic
        'orbital_period': 27.3217,  # days
        'rotation_period': 27.3217,  # days (tidally locked)
        'axial_tilt': 6.68,  # degrees (to ecliptic)
        'color': (0.7, 0.7, 0.7),
        'texture': 'textures/moon.jpg'
    },
    'Io': {  # Jupiter's moon (tidally locked)
        'name': 'Io',
        'parent': 'Jupiter',
        'mass': 0.08932e24,  # kg
        'radius': 1821.6 * 1e3,  # meters
        'semi_major_axis': 421.8e6,  # meters
        'eccentricity': 0.0041,
        'inclination': math.radians(0.05),
        'orbital_period': 1.769,  # days
        'rotation_period': 1.769,  # days (tidally locked)
        'axial_tilt': 0.0,  # degrees
        'color': (0.9, 0.8, 0.3),
        'texture': 'textures/moon.jpg'
    },
    'Europa': {  # Jupiter's moon (tidally locked)
        'name': 'Europa',
        'parent': 'Jupiter',
        'mass': 0.04800e24,  # kg
        'radius': 1560.8 * 1e3,  # meters
        'semi_major_axis': 671.1e6,  # meters
        'eccentricity': 0.0094,
        'inclination': math.radians(0.47),
        'orbital_period': 3.551,  # days
        'rotation_period': 3.551,  # days (tidally locked)
        'axial_tilt': 0.1,  # degrees
        'color': (0.8, 0.8, 0.7),
        'texture': 'textures/moon.jpg'
    },
    'Ganymede': {  # Jupiter's moon (largest moon in solar system, tidally locked)
        'name': 'Ganymede',
        'parent': 'Jupiter',
        'mass': 0.1482e24,  # kg
        'radius': 2634.1 * 1e3,  # meters
        'semi_major_axis': 1070.4e6,  # meters
        'eccentricity': 0.0013,
        'inclination': math.radians(0.20),
        'orbital_period': 7.155,  # days
        'rotation_period': 7.155,  # days (tidally locked)
        'axial_tilt': 0.165,  # degrees
        'color': (0.6, 0.6, 0.5),
        'texture': 'textures/moon.jpg'
    },
    'Callisto': {  # Jupiter's moon (tidally locked)
        'name': 'Callisto',
        'parent': 'Jupiter',
        'mass': 0.1076e24,  # kg
        'radius': 2410.3 * 1e3,  # meters
        'semi_major_axis': 1882.7e6,  # meters
        'eccentricity': 0.0074,
        'inclination': math.radians(0.51),
        'orbital_period': 16.689,  # days
        'rotation_period': 16.689,  # days (tidally locked)
        'axial_tilt': 0.0,  # degrees
        'color': (0.5, 0.5, 0.4),
        'texture': 'textures/moon.jpg'
    },
    'Titan': {  # Saturn's moon (2nd largest in solar system, tidally locked)
        'name': 'Titan',
        'parent': 'Saturn',
        'mass': 0.13452e24,  # kg
        'radius': 2574.7 * 1e3,  # meters
        'semi_major_axis': 1221.87e6,  # meters
        'eccentricity': 0.0288,
        'inclination': math.radians(0.34854),
        'orbital_period': 15.945,  # days
        'rotation_period': 15.945,  # days (tidally locked)
        'axial_tilt': 0.3,  # degrees
        'color': (0.8, 0.6, 0.4),
        'texture': 'textures/moon.jpg'
    }
}


def generate_asteroid_belt(count=100, seed=42):
    """Generate asteroid belt data between Mars and Jupiter orbits (2.1-3.3 AU)"""
    rng = np.random.RandomState(seed)
    asteroids = []

    # Asteroid belt region: 2.1 AU to 3.3 AU
    min_au = 2.1
    max_au = 3.3
    au_to_meters = 149.597870700e9  # 1 AU in meters

    # Kirkwood gaps (orbital resonances with Jupiter to avoid)
    # 2.50 AU (3:1), 2.82 AU (5:2), 2.95 AU (7:3), 3.27 AU (2:1)
    kirkwood_gaps = [2.50, 2.82, 2.95, 3.27]
    gap_width = 0.03  # AU

    for i in range(count):
        # Random semi-major axis avoiding Kirkwood gaps
        while True:
            a_au = rng.uniform(min_au, max_au)
            in_gap = any(abs(a_au - gap) < gap_width for gap in kirkwood_gaps)
            if not in_gap:
                break

        semi_major_axis = a_au * au_to_meters  # meters

        # Kepler's third law: period² = a³ (in AU and years)
        period_years = a_au ** 1.5
        period_days = period_years * 365.25

        # Random orbital elements
        eccentricity = rng.uniform(0.0, 0.3)
        inclination = math.radians(rng.uniform(0, 20))
        arg_periapsis = math.radians(rng.uniform(0, 360))
        phase = rng.uniform(0, 2 * math.pi)

        # Mass: ranges from small rocks to Ceres-sized (~9.4e20 kg)
        log_mass = rng.uniform(12, 18)  # 1e12 to 1e18 kg
        mass = 10 ** log_mass

        # Radius: 0.5 km to 500 km
        log_radius = rng.uniform(2.7, 5.7)  # log10(500) to log10(500000) meters
        radius = 10 ** log_radius

        # Gray/brown colors
        gray = rng.uniform(0.3, 0.7)
        color = (gray + rng.uniform(-0.1, 0.1),
                 gray + rng.uniform(-0.1, 0.05),
                 gray + rng.uniform(-0.15, 0.0))

        asteroids.append({
            'name': f'Asteroid-{i+1:03d}',
            'mass': mass,
            'radius': radius,
            'semi_major_axis': semi_major_axis,
            'eccentricity': eccentricity,
            'inclination': inclination,
            'orbital_period': period_days,
            'rotation_period': rng.uniform(0.1, 2.0),  # days (most spin fast)
            'axial_tilt': rng.uniform(0, 180),  # random tilt
            'arg_periapsis': arg_periapsis,
            'phase': phase,
            'color': color,
            'texture': 'textures/moon.jpg'
        })

    return asteroids


def get_scaled_radius(radius_m):
    """Convert real radius in meters to scaled visualization radius"""
    return radius_m * RADIUS_SCALE


def get_scaled_distance(distance_m):
    """Convert real distance in meters to scaled visualization distance"""
    return distance_m * DISTANCE_SCALE


def get_orbital_period_seconds(period_days):
    """Convert orbital period from days to seconds"""
    return period_days * 24.0 * 3600.0
