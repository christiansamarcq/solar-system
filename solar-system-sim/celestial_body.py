"""
Celestial Body class for solar system simulation
Represents planets, moons, and the sun
"""

from panda3d.core import NodePath, GeomVertexFormat, GeomVertexData, GeomVertexWriter
from panda3d.core import Geom, GeomTriangles, GeomNode, Texture, TextureStage
from direct.showbase.ShowBase import ShowBase
import numpy as np
import math
from physics.simple_orbit import elliptical_orbit_position, circular_orbit_position


def create_sphere_geometry(radius=1.0, segments=32, rings=16):
    """
    Create a high-resolution sphere geometry with UV texture coordinates

    Args:
        radius (float): Radius of the sphere
        segments (int): Number of horizontal segments (longitude)
        rings (int): Number of vertical rings (latitude)

    Returns:
        GeomNode: A geometry node containing the sphere
    """
    # Create vertex data format with UV coordinates
    format = GeomVertexFormat.getV3n3t2()  # Position, Normal, Texture coords
    vdata = GeomVertexData('sphere', format, Geom.UHStatic)

    vertex = GeomVertexWriter(vdata, 'vertex')
    normal = GeomVertexWriter(vdata, 'normal')
    texcoord = GeomVertexWriter(vdata, 'texcoord')

    # Generate vertices
    for ring in range(rings + 1):
        theta = ring * math.pi / rings  # Vertical angle (0 to pi)
        sin_theta = math.sin(theta)
        cos_theta = math.cos(theta)

        for seg in range(segments + 1):
            phi = seg * 2 * math.pi / segments  # Horizontal angle (0 to 2pi)
            sin_phi = math.sin(phi)
            cos_phi = math.cos(phi)

            # Calculate vertex position
            x = radius * sin_theta * cos_phi
            y = radius * sin_theta * sin_phi
            z = radius * cos_theta

            # Vertex position
            vertex.addData3(x, y, z)

            # Normal (same as normalized position for a sphere)
            normal.addData3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta)

            # UV texture coordinates (0 to 1)
            u = seg / segments
            v = ring / rings
            texcoord.addData2(u, v)

    # Create triangles
    tris = GeomTriangles(Geom.UHStatic)

    for ring in range(rings):
        for seg in range(segments):
            # Current ring vertices
            current = ring * (segments + 1) + seg
            next_seg = ring * (segments + 1) + (seg + 1)

            # Next ring vertices
            next_ring = (ring + 1) * (segments + 1) + seg
            next_ring_next = (ring + 1) * (segments + 1) + (seg + 1)

            # First triangle
            tris.addVertices(current, next_ring, next_seg)

            # Second triangle
            tris.addVertices(next_seg, next_ring, next_ring_next)

    tris.closePrimitive()

    # Create geometry
    geom = Geom(vdata)
    geom.addPrimitive(tris)

    # Create node
    node = GeomNode('sphere')
    node.addGeom(geom)

    return node


def create_asteroid_geometry(radius=1.0, seed=0, segments=48, rings=24):
    """
    Create an irregular, potato-shaped rock with craters and lumps.

    Starts from a UV sphere and pushes each vertex in or out:
      - uneven axes (a, b, c) so it is elongated / flattened, not round
      - broad lumps and dents for a lopsided outline
      - bowl-shaped craters with slightly raised rims
      - many small bumps for a rough, rocky surface
    The displacement depends only on direction, so the texture seam and the
    poles stay closed. Normals are recomputed from the bumpy surface so the
    sun light picks out the craters.

    Args:
        radius (float): Mean radius of the rock
        seed (int): Each seed gives a different, repeatable shape
        segments (int): Horizontal detail
        rings (int): Vertical detail

    Returns:
        GeomNode: A geometry node containing the asteroid
    """
    rng = np.random.default_rng(seed)

    # Unit directions on a (rings+1) x (segments+1) grid
    theta = np.linspace(0, math.pi, rings + 1)[:, None]
    phi = np.linspace(0, 2 * math.pi, segments + 1)[None, :]
    dirs = np.stack([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta) * np.ones_like(phi),
    ], axis=-1)                                   # shape (R, S, 3)

    # 1) Uneven axes: long, medium, short
    axes = np.array([1.0, rng.uniform(0.5, 0.85), rng.uniform(0.35, 0.7)])
    axes /= np.cbrt(np.prod(axes))                # keep the same volume as the sphere

    # 2) Broad lumps and dents
    r = np.ones(dirs.shape[:2])
    for _ in range(rng.integers(8, 14)):
        center = rng.normal(size=3); center /= np.linalg.norm(center)
        width = rng.uniform(0.3, 0.9)             # angular size (radians)
        height = rng.uniform(-0.25, 0.25)
        ang = np.arccos(np.clip(dirs @ center, -1, 1))
        r += height * np.exp(-(ang / width) ** 2)

    # 3) Craters: a dip in the middle with a slightly raised rim
    for _ in range(rng.integers(6, 12)):
        center = rng.normal(size=3); center /= np.linalg.norm(center)
        size = rng.uniform(0.22, 0.5)
        depth = rng.uniform(0.08, 0.18) * (size / 0.5) ** 0.5
        x = np.arccos(np.clip(dirs @ center, -1, 1)) / size
        r += np.where(x < 1.0, -depth * (1 - x ** 2), 0.0)
        r += np.where((x >= 1.0) & (x < 1.6), depth * 0.35 * np.sin((x - 1.0) / 0.6 * math.pi), 0.0)

    # 4) Small bumps all over for a rough, rocky texture. They are kept wider
    #    than the mesh spacing (~0.13 rad); narrower ones alias into spikes.
    for _ in range(40):
        center = rng.normal(size=3); center /= np.linalg.norm(center)
        width = rng.uniform(0.16, 0.3)
        height = rng.uniform(-0.045, 0.045)
        ang = np.arccos(np.clip(dirs @ center, -1, 1))
        r += height * np.exp(-(ang / width) ** 2)

    r = np.maximum(r, 0.3)                        # never fold inside out
    pos = dirs * axes * r[..., None] * radius     # (R, S, 3)

    # Normals from the surface: cross product of neighbouring grid directions
    seg_next = np.roll(pos[:, :-1], -1, axis=1)
    seg_prev = np.roll(pos[:, :-1], 1, axis=1)
    d_phi = np.concatenate([seg_next - seg_prev, (seg_next - seg_prev)[:, :1]], axis=1)
    d_theta = np.empty_like(pos)
    d_theta[1:-1] = pos[2:] - pos[:-2]
    d_theta[0] = pos[1] - pos[0]
    d_theta[-1] = pos[-1] - pos[-2]
    normals = np.cross(d_theta, d_phi)
    lengths = np.linalg.norm(normals, axis=-1, keepdims=True)
    # Poles are degenerate on a UV sphere: fall back to the radial direction
    normals = np.where(lengths > 1e-12, normals / np.maximum(lengths, 1e-12), dirs)
    # At each pole use the average of the ring next to it, which follows the
    # bumpy surface there (a fixed up/down normal lights the pole wrongly)
    for pole, ring_next in ((0, 1), (-1, -2)):
        avg = normals[ring_next, :-1].mean(axis=0)
        normals[pole, :] = avg / max(np.linalg.norm(avg), 1e-12)
    # Make sure every normal points outward
    flip = np.sum(normals * dirs, axis=-1) < 0
    normals[flip] *= -1

    format = GeomVertexFormat.getV3n3t2()
    vdata = GeomVertexData('asteroid', format, Geom.UHStatic)
    vertex = GeomVertexWriter(vdata, 'vertex')
    normal = GeomVertexWriter(vdata, 'normal')
    texcoord = GeomVertexWriter(vdata, 'texcoord')
    for ring in range(rings + 1):
        for seg in range(segments + 1):
            p = pos[ring, seg]
            n = normals[ring, seg]
            vertex.addData3(p[0], p[1], p[2])
            normal.addData3(n[0], n[1], n[2])
            texcoord.addData2(seg / segments, ring / rings)

    tris = GeomTriangles(Geom.UHStatic)
    for ring in range(rings):
        for seg in range(segments):
            current = ring * (segments + 1) + seg
            next_seg = current + 1
            next_ring = current + segments + 1
            next_ring_next = next_ring + 1
            tris.addVertices(current, next_ring, next_seg)
            tris.addVertices(next_seg, next_ring, next_ring_next)
    tris.closePrimitive()

    geom = Geom(vdata)
    geom.addPrimitive(tris)
    node = GeomNode('asteroid')
    node.addGeom(geom)
    return node


class CelestialBody:
    """Represents a celestial body (sun, planet, moon) in the simulation"""

    # Floating origin: visuals are drawn relative to this point (usually the
    # tracked body). Panda3D stores positions as 32-bit floats, which at
    # thousands of units from the Sun can't resolve objects smaller than
    # ~0.0005 units, so far-out bodies jitter and tear. Subtracting the
    # origin here in 64-bit numpy keeps everything near the camera precise.
    render_origin = np.zeros(3)

    def __init__(self, name, mass, radius, color, position=None, velocity=None, is_emissive=False, texture_path=None):
        """
        Initialize a celestial body

        Args:
            name (str): Name of the body (e.g., "Earth", "Sun")
            mass (float): Mass in arbitrary units
            radius (float): Visual radius for rendering
            color (tuple): RGB color values (0-1 range), e.g., (1, 0, 0) for red
            position (np.array): Initial 3D position [x, y, z], defaults to origin
            velocity (np.array): Initial 3D velocity [vx, vy, vz], defaults to zero
            is_emissive (bool): If True, object glows (not affected by lighting)
            texture_path (str): Path to texture image file
        """
        self.name = name
        self.mass = mass
        self.radius = radius
        self.color = color
        self.is_emissive = is_emissive
        self.texture_path = texture_path

        # Physics properties
        self.position = position if position is not None else np.array([0.0, 0.0, 0.0])
        self.velocity = velocity if velocity is not None else np.array([0.0, 0.0, 0.0])
        self.acceleration = np.array([0.0, 0.0, 0.0])

        # Orbital parameters for simple mode
        self.orbital_distance = 0.0
        self.orbital_period = 0.0
        self.orbital_phase = 0.0  # Starting angle

        # Elliptical orbit parameters
        self.semi_major_axis = 0.0  # a (average distance)
        self.eccentricity = 0.0  # e (0=circle, 0.5=ellipse, 1=parabola)
        self.inclination = 0.0  # i (tilt angle in radians)
        self.arg_periapsis = 0.0  # ω (rotation of ellipse)

        # Rotation properties
        self.rotation_period = 0.0  # days (negative = retrograde)
        self.axial_tilt = 0.0  # degrees
        self.spin_angle = 0.0  # current rotation angle in degrees
        self.brightness = 2.0  # colour boost so bodies show up against black space

        # Energy flux properties
        self.energy_received = 0.0  # Energy flux from sun
        self.base_color = color  # Store original color

        # 3D model node (will be created when attached to scene)
        self.node = None

    def create_visual(self, app, segments=64, rings=32):
        """
        Create the 3D visual representation

        Args:
            app: Panda3D ShowBase application instance
            segments (int): Horizontal detail (more = smoother)
            rings (int): Vertical detail (more = smoother)
        """
        # Irregular rock for asteroids, smooth sphere for everything else
        shape_seed = getattr(self, 'shape_seed', None)
        if shape_seed is not None:
            sphere_geom = create_asteroid_geometry(radius=self.radius, seed=shape_seed)
        else:
            sphere_geom = create_sphere_geometry(radius=self.radius, segments=segments, rings=rings)

        # Attach to scene
        self.node = app.render.attachNewNode(sphere_geom)

        # Apply texture if available
        if self.texture_path:
            try:
                tex = app.loader.loadTexture(self.texture_path)
                self.node.setTexture(tex)
                # Set white color to not tint the texture. Remember it as the
                # base so the heat-map reset doesn't tint it with self.color.
                self.node.setColor(1, 1, 1, 1)
                self.base_color = (1.0, 1.0, 1.0)
            except Exception as e:
                print(f"Warning: Could not load texture {self.texture_path}: {e}")
                # Fallback to solid color
                self.node.setColor(self.color[0], self.color[1], self.color[2], 1)
        else:
            # Set color if no texture
            self.node.setColor(self.color[0], self.color[1], self.color[2], 1)

        # If emissive (like the sun), make it bright and unaffected by lighting
        if self.is_emissive:
            # Disable lighting so it always appears bright
            self.node.setLightOff()
            # Boost brightness
            self.node.setColorScale(1.5, 1.5, 1.2, 1)  # Extra bright yellow-white
        else:
            # For planets/moons, increase base brightness to make them visible against black space
            b = self.brightness
            self.node.setColorScale(b, b, b, 1)

        # Set initial position
        self.update_visual_position()

    def update_visual_position(self):
        """Update the 3D model position and rotation to match physics"""
        if self.node:
            rel = self.position - CelestialBody.render_origin
            self.node.setPos(rel[0], rel[1], rel[2])
            # Apply axial tilt (pitch) and spin (heading)
            # H = spin around axis, P = axial tilt, R = 0
            self.node.setHpr(self.spin_angle, self.axial_tilt, 0)

    def update_spin(self, dt_seconds):
        """Update spin angle based on rotation period"""
        if self.rotation_period != 0:
            # rotation_period is in days, dt is in seconds
            period_seconds = abs(self.rotation_period) * 86400.0
            degrees_per_second = 360.0 / period_seconds
            # Negative rotation_period = retrograde = spin in opposite direction
            if self.rotation_period < 0:
                degrees_per_second = -degrees_per_second
            self.spin_angle += degrees_per_second * dt_seconds
            self.spin_angle %= 360.0

    def set_simple_orbit(self, distance, period, phase=0.0):
        """
        Configure simple circular orbit parameters

        Args:
            distance (float): Distance from center (orbital radius)
            period (float): Time for one complete orbit
            phase (float): Starting angle in radians
        """
        self.orbital_distance = distance
        self.orbital_period = period
        self.orbital_phase = phase

        # Set elliptical parameters for circular orbit
        self.semi_major_axis = distance
        self.eccentricity = 0.0

    def set_elliptical_orbit(self, semi_major_axis, eccentricity, period,
                            inclination=0.0, arg_periapsis=0.0, phase=0.0):
        """
        Configure elliptical orbit parameters

        Args:
            semi_major_axis (float): a - average orbital radius
            eccentricity (float): e - 0=circle, 0.5=ellipse (0 to <1)
            period (float): Orbital period
            inclination (float): Tilt angle in radians (0=flat)
            arg_periapsis (float): Rotation of ellipse in radians
            phase (float): Starting position along orbit
        """
        self.semi_major_axis = semi_major_axis
        self.eccentricity = eccentricity
        self.orbital_period = period
        self.inclination = inclination
        self.arg_periapsis = arg_periapsis
        self.orbital_phase = phase

        # For compatibility
        self.orbital_distance = semi_major_axis

    def update_simple_orbit(self, time):
        """
        Update position based on orbital parameters (circular or elliptical)

        Args:
            time (float): Current simulation time
        """
        if self.orbital_period > 0:
            # Use elliptical orbit calculation if eccentricity > 0 or inclination > 0
            if self.eccentricity > 0.001 or abs(self.inclination) > 0.001:
                self.position = elliptical_orbit_position(
                    time,
                    self.semi_major_axis,
                    self.eccentricity,
                    self.orbital_period,
                    self.inclination,
                    self.arg_periapsis,
                    self.orbital_phase
                )
            else:
                # Fast path for circular orbits
                self.position = circular_orbit_position(
                    time,
                    self.orbital_distance,
                    self.orbital_period,
                    self.orbital_phase
                )

            self.update_visual_position()

    def apply_force(self, force):
        """
        Apply a force to this body (for N-body simulation)

        Args:
            force (np.array): Force vector [fx, fy, fz]
        """
        # F = ma, so a = F/m
        self.acceleration += force / self.mass

    def update_physics(self, dt):
        """
        Update position and velocity based on current acceleration (N-body mode)

        Args:
            dt (float): Time step
        """
        # Velocity Verlet integration
        self.velocity += self.acceleration * dt
        self.position += self.velocity * dt

        # Reset acceleration for next frame
        self.acceleration = np.array([0.0, 0.0, 0.0])

        self.update_visual_position()

    def calculate_energy_flux(self, sun_position, sun_luminosity=1000.0):
        """
        Calculate energy flux received from the sun (inverse-square law)

        E = L / (4π * r²)

        Args:
            sun_position (np.array): Position of the sun
            sun_luminosity (float): Luminosity of the sun (arbitrary units)

        Returns:
            float: Energy flux received
        """
        # Distance to sun
        r_vec = self.position - sun_position
        distance = np.linalg.norm(r_vec)

        if distance < 0.1:  # Minimum distance to prevent singularity
            distance = 0.1

        # Inverse square law
        self.energy_received = sun_luminosity / (4 * np.pi * distance * distance)

        return self.energy_received

    def update_heat_visualization(self, energy_vis_enabled=True):
        """
        Update planet color based on energy received
        Creates a heat map effect - closer to sun = warmer colors

        Args:
            energy_vis_enabled (bool): Whether to apply heat visualization
        """
        if not self.node or self.is_emissive:
            return

        if not energy_vis_enabled:
            # Reset to original color
            self.node.setColor(self.base_color[0], self.base_color[1], self.base_color[2], 1)
            return

        # Normalize energy (0 = cold, 1 = hot)
        # Typical range: 0.1 to 10.0 arbitrary units
        energy_normalized = min(1.0, max(0.0, (self.energy_received - 0.1) / 10.0))

        # Create heat gradient
        # Cold (blue tint) -> Warm (red/orange tint)
        if energy_normalized < 0.5:
            # Cold side: add blue, reduce red
            red_boost = 0.7
            green_boost = 0.8
            blue_boost = 1.0 + (1.0 - energy_normalized * 2) * 0.5
        else:
            # Hot side: add red/orange
            red_boost = 1.0 + (energy_normalized - 0.5) * 1.0
            green_boost = 1.0 + (energy_normalized - 0.5) * 0.3
            blue_boost = 1.0 - (energy_normalized - 0.5) * 0.5

        # Apply heat tint to base color
        heated_color = (
            min(1.0, self.base_color[0] * red_boost),
            min(1.0, self.base_color[1] * green_boost),
            min(1.0, self.base_color[2] * blue_boost)
        )

        self.node.setColor(heated_color[0], heated_color[1], heated_color[2], 1)

    def __repr__(self):
        return f"CelestialBody({self.name}, mass={self.mass}, pos={self.position})"
