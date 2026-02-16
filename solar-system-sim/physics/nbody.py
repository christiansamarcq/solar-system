"""
N-body gravitational simulation
Calculates realistic gravitational forces between all celestial bodies

IMPORTANT: Positions are stored in scaled units (1 unit = 1 billion meters)
but physics calculations need real meters, so we convert using DISTANCE_SCALE
"""

import numpy as np
import sys
import os

# Import the real gravitational constant and scaling factors
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import realistic_data as rd

# Use the real gravitational constant from NASA data
G = rd.G  # 6.67430e-11 m^3 kg^-1 s^-2

# Distance conversion: positions are in scaled units, need to convert to meters for physics
METERS_PER_UNIT = 1.0 / rd.DISTANCE_SCALE  # 1e9 meters per unit

# Minimum distance to prevent gravitational singularity
MIN_DISTANCE = 0.1


def calculate_gravitational_force(body1_pos, body1_mass, body2_pos, body2_mass):
    """
    Calculate gravitational force between two bodies

    F = G * m1 * m2 / r²

    Args:
        body1_pos (np.array): Position of first body [x, y, z]
        body1_mass (float): Mass of first body
        body2_pos (np.array): Position of second body [x, y, z]
        body2_mass (float): Mass of second body

    Returns:
        np.array: Force vector on body1 due to body2
    """
    # Vector from body1 to body2
    r_vec = body2_pos - body1_pos

    # Distance between bodies
    r = np.linalg.norm(r_vec)

    # Avoid division by zero
    if r < MIN_DISTANCE:
        r = MIN_DISTANCE

    # Direction (unit vector)
    r_hat = r_vec / r

    # Gravitational force magnitude
    F_magnitude = G * body1_mass * body2_mass / (r * r)

    # Force vector (points toward body2)
    F_vec = F_magnitude * r_hat

    return F_vec


def calculate_all_forces(bodies):
    """
    Calculate net gravitational force on each body from all others

    Args:
        bodies (list): List of CelestialBody objects

    Returns:
        list: List of force vectors for each body
    """
    n = len(bodies)
    forces = [np.array([0.0, 0.0, 0.0]) for _ in range(n)]

    # Calculate pairwise forces (O(n²) complexity)
    for i in range(n):
        for j in range(i + 1, n):
            # Force on body i due to body j
            force_i = calculate_gravitational_force(
                bodies[i].position, bodies[i].mass,
                bodies[j].position, bodies[j].mass
            )

            # Newton's 3rd law: force on j is opposite
            force_j = -force_i

            # Accumulate forces
            forces[i] += force_i
            forces[j] += force_j

    return forces


def calculate_all_forces_vectorized(bodies):
    """
    Calculate net gravitational force on each body using vectorized NumPy operations

    This is 10-100x faster than the loop version for moderate numbers of bodies.
    Uses broadcasting to compute all pairwise forces simultaneously.

    Args:
        bodies (list): List of CelestialBody objects

    Returns:
        np.array: Array of force vectors for each body (shape: n x 3)
    """
    n = len(bodies)

    if n == 0:
        return np.array([])

    # Extract positions and masses into arrays
    positions = np.array([body.position for body in bodies])  # shape: (n, 3) in scaled units
    masses = np.array([body.mass for body in bodies])  # shape: (n,)

    # Convert positions from scaled units to meters for physics calculations
    positions_meters = positions * METERS_PER_UNIT  # shape: (n, 3) in meters

    # Compute all pairwise displacement vectors using broadcasting
    # positions_meters[:, np.newaxis, :] has shape (n, 1, 3)
    # positions_meters[np.newaxis, :, :] has shape (1, n, 3)
    # Broadcasting creates (n, n, 3) where [i, j, :] is vector from i to j
    r_vec = positions_meters[np.newaxis, :, :] - positions_meters[:, np.newaxis, :]  # (n, n, 3) in meters

    # Compute all pairwise distances (in meters)
    # Norm along axis 2 gives distance for each pair
    r = np.linalg.norm(r_vec, axis=2)  # shape: (n, n) in meters

    # Avoid division by zero - set minimum distance (0.1 billion meters = 100 million km)
    r = np.maximum(r, 0.1 * METERS_PER_UNIT)  # Minimum 0.1 scaled units = 100 million meters

    # Direction unit vectors (broadcasting r to match r_vec shape)
    r_hat = r_vec / r[:, :, np.newaxis]  # shape: (n, n, 3)

    # Gravitational force magnitudes for all pairs
    # masses[:, np.newaxis] broadcasts to (n, 1)
    # masses[np.newaxis, :] broadcasts to (1, n)
    # Result is (n, n) where [i, j] is force magnitude on i due to j
    F_magnitude = G * masses[:, np.newaxis] * masses[np.newaxis, :] / (r * r)  # (n, n)

    # Force vectors for all pairs
    # Broadcast F_magnitude to match r_hat shape
    F_vec = F_magnitude[:, :, np.newaxis] * r_hat  # shape: (n, n, 3)

    # Zero out diagonal (self-forces) to avoid NaN issues
    for dim in range(3):
        np.fill_diagonal(F_vec[:, :, dim], 0)

    # Sum forces on each body (sum over j axis, which is axis 1)
    forces = np.sum(F_vec, axis=1)  # shape: (n, 3)

    return forces


def update_nbody_physics(bodies, dt, size_scale=1.0, sun_scale=1.0):
    """
    Update positions and velocities using N-body gravity

    Uses Velocity Verlet integration for better energy conservation
    Uses vectorized force calculations for better performance

    Args:
        bodies (list): List of CelestialBody objects
        dt (float): Time step
        size_scale (float): Visual size multiplier for planet collision detection (default 1.0)
        sun_scale (float): Visual size multiplier for sun collision detection (default 1.0)

    Returns:
        list: List of collision events [(body1, body2), ...]
    """
    # Calculate forces on all bodies using vectorized version
    forces = calculate_all_forces_vectorized(bodies)

    # Extract masses for vectorized operations
    masses = np.array([body.mass for body in bodies])

    # Calculate accelerations: a = F/m
    accelerations = forces / masses[:, np.newaxis]  # Broadcasting mass to 3D

    # Velocity Verlet step 1: Update positions
    for i, body in enumerate(bodies):
        # Skip bodies with infinite mass (like the sun, if desired)
        if body.mass > 1e29:  # Treat very massive objects as fixed (sun only)
            continue

        # Velocity Verlet: v(t+dt/2) = v(t) + a(t)*dt/2 (velocity in m/s)
        half_step_velocity = body.velocity + accelerations[i] * (dt / 2)

        # Check for NaN/Inf in acceleration before updating position
        if np.any(np.isnan(accelerations[i])) or np.any(np.isinf(accelerations[i])):
            print(f"ERROR: Invalid acceleration for {body.name}: {accelerations[i]}")
            return []  # Abort physics update

        # Update position: x(t+dt) = x(t) + v(t+dt/2)*dt
        # Velocity is in m/s, but position is in scaled units
        # Convert displacement from meters to scaled units
        displacement_meters = half_step_velocity * dt
        displacement_units = displacement_meters / METERS_PER_UNIT
        new_position = body.position + displacement_units

        # Check for NaN/Inf in new position
        if np.any(np.isnan(new_position)) or np.any(np.isinf(new_position)):
            print(f"ERROR: Invalid position for {body.name}: {new_position}")
            return []  # Abort physics update

        body.position = new_position

    # Recalculate forces at new positions
    new_forces = calculate_all_forces_vectorized(bodies)
    new_accelerations = new_forces / masses[:, np.newaxis]

    # Velocity Verlet step 2: Complete velocity update
    for i, body in enumerate(bodies):
        if body.mass > 1e29:
            continue

        # Velocity Verlet: v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2
        new_velocity = body.velocity + (accelerations[i] + new_accelerations[i]) * (dt / 2)

        # Check for NaN/Inf in new velocity BEFORE updating visual
        if np.any(np.isnan(new_velocity)) or np.any(np.isinf(new_velocity)):
            print(f"ERROR: Invalid velocity for {body.name}: {new_velocity}")
            return []  # Abort physics update

        body.velocity = new_velocity

        # Only update visual if all values are valid
        body.update_visual_position()

    # Collision detection: check if any bodies are too close (using VISUAL radius)
    collisions = []
    n = len(bodies)
    for i in range(n):
        for j in range(i + 1, n):
            # Calculate distance between bodies (in scaled units)
            r_vec = bodies[j].position - bodies[i].position
            distance = np.linalg.norm(r_vec)

            # Sum of VISUAL radii (each body uses its own scale)
            uses_sun_i = getattr(bodies[i], 'uses_sun_scale', bodies[i].is_emissive)
            uses_sun_j = getattr(bodies[j], 'uses_sun_scale', bodies[j].is_emissive)
            scale_i = sun_scale if uses_sun_i else size_scale
            scale_j = sun_scale if uses_sun_j else size_scale
            collision_radius = bodies[i].radius * scale_i + bodies[j].radius * scale_j

            # Check for collision
            if distance < collision_radius:
                collisions.append((bodies[i], bodies[j]))
                print(f"COLLISION DETECTED: {bodies[i].name} and {bodies[j].name} (distance: {distance:.3f}, threshold: {collision_radius:.3f})")

    return collisions


def initialize_circular_orbit_velocity(body, central_mass, central_pos):
    """
    Set initial velocity for stable circular orbit around a central body

    Args:
        body: CelestialBody to initialize
        central_mass (float): Mass of central body
        central_pos (np.array): Position of central body
    """
    # Distance to central body (in scaled units)
    r_vec = body.position - central_pos
    r_units = np.linalg.norm(r_vec)

    if r_units < 0.1:
        return

    # Convert distance to meters for physics calculation
    r_meters = r_units * METERS_PER_UNIT

    # Orbital speed for circular orbit: v = sqrt(G*M/r) (in m/s)
    v = np.sqrt(G * central_mass / r_meters)

    # Perpendicular direction (in XY plane for now)
    # Rotate position vector by 90 degrees in XY plane
    perp = np.array([-r_vec[1], r_vec[0], 0.0])
    perp_norm = np.linalg.norm(perp)

    # Check if perpendicular vector is valid (not on Z-axis)
    if perp_norm < 1e-10:
        print(f"WARNING: {body.name} is on Z-axis, using fallback velocity direction")
        # Use a default perpendicular direction if on Z-axis
        perp = np.array([1.0, 0.0, 0.0])
    else:
        perp = perp / perp_norm

    # Set velocity
    body.velocity = v * perp
    print(f"Initialized {body.name}: v={v:.6f}, velocity={body.velocity}")


def calculate_total_energy(bodies):
    """
    Calculate total energy (kinetic + potential) of the system

    Useful for debugging and checking energy conservation

    Args:
        bodies (list): List of CelestialBody objects

    Returns:
        float: Total energy of the system
    """
    kinetic = 0.0
    potential = 0.0

    # Kinetic energy
    for body in bodies:
        v_squared = np.dot(body.velocity, body.velocity)
        kinetic += 0.5 * body.mass * v_squared

    # Potential energy (pairwise)
    n = len(bodies)
    for i in range(n):
        for j in range(i + 1, n):
            r_vec = bodies[j].position - bodies[i].position
            r = np.linalg.norm(r_vec)
            if r > 0.1:
                potential -= G * bodies[i].mass * bodies[j].mass / r

    return kinetic + potential
