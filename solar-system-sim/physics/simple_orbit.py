"""
Simple orbital mechanics calculations
Handles elliptical orbits using Kepler's laws
"""

import numpy as np
import math


def solve_kepler_equation(mean_anomaly, eccentricity, tolerance=1e-6, max_iterations=100):
    """
    Solve Kepler's equation to find eccentric anomaly
    M = E - e*sin(E)

    Args:
        mean_anomaly (float): M - angle if orbit were circular (radians)
        eccentricity (float): e - orbital eccentricity (0 to 1)
        tolerance (float): Convergence tolerance
        max_iterations (int): Maximum Newton-Raphson iterations

    Returns:
        float: E - eccentric anomaly (radians)
    """
    # Initial guess
    E = mean_anomaly if eccentricity < 0.8 else math.pi

    # Newton-Raphson iteration
    for _ in range(max_iterations):
        f = E - eccentricity * math.sin(E) - mean_anomaly
        f_prime = 1 - eccentricity * math.cos(E)

        E_new = E - f / f_prime

        if abs(E_new - E) < tolerance:
            return E_new

        E = E_new

    return E


def elliptical_orbit_position(time, semi_major_axis, eccentricity, period,
                              inclination=0.0, arg_periapsis=0.0, phase=0.0):
    """
    Calculate 3D position for elliptical orbit

    Args:
        time (float): Current time
        semi_major_axis (float): a - average orbital radius
        eccentricity (float): e - orbital eccentricity
        period (float): Orbital period
        inclination (float): i - orbit tilt (radians)
        arg_periapsis (float): ω - ellipse rotation (radians)
        phase (float): Starting phase offset (radians)

    Returns:
        np.array: 3D position [x, y, z]
    """
    # Mean motion (angular velocity for circular orbit)
    n = 2 * math.pi / period

    # Mean anomaly (angle in circular orbit)
    M = n * time + phase

    # Eccentric anomaly (solve Kepler's equation)
    E = solve_kepler_equation(M, eccentricity)

    # True anomaly (actual angle in elliptical orbit)
    true_anomaly = 2 * math.atan2(
        math.sqrt(1 + eccentricity) * math.sin(E / 2),
        math.sqrt(1 - eccentricity) * math.cos(E / 2)
    )

    # Distance from focus (sun)
    r = semi_major_axis * (1 - eccentricity * math.cos(E))

    # Position in orbital plane (before rotation)
    x_orbit = r * math.cos(true_anomaly)
    y_orbit = r * math.sin(true_anomaly)

    # Apply argument of periapsis (rotate ellipse in its plane)
    cos_w = math.cos(arg_periapsis)
    sin_w = math.sin(arg_periapsis)

    x_rotated = x_orbit * cos_w - y_orbit * sin_w
    y_rotated = x_orbit * sin_w + y_orbit * cos_w

    # Apply inclination (tilt orbit out of XY plane)
    cos_i = math.cos(inclination)
    sin_i = math.sin(inclination)

    x = x_rotated
    y = y_rotated * cos_i
    z = y_rotated * sin_i

    return np.array([x, y, z])


def circular_orbit_position(time, distance, period, phase=0.0):
    """
    Calculate position for simple circular orbit (fast path)

    Args:
        time (float): Current time
        distance (float): Orbital radius
        period (float): Orbital period
        phase (float): Starting angle

    Returns:
        np.array: 3D position [x, y, z]
    """
    omega = 2 * math.pi / period
    angle = omega * time + phase

    x = distance * math.cos(angle)
    y = distance * math.sin(angle)
    z = 0.0

    return np.array([x, y, z])
