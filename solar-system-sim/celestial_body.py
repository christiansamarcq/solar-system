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


class CelestialBody:
    """Represents a celestial body (sun, planet, moon) in the simulation"""

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
        # Create high-resolution sphere geometry
        sphere_geom = create_sphere_geometry(radius=self.radius, segments=segments, rings=rings)

        # Attach to scene
        self.node = app.render.attachNewNode(sphere_geom)

        # Apply texture if available
        if self.texture_path:
            try:
                tex = app.loader.loadTexture(self.texture_path)
                self.node.setTexture(tex)
                # Set white color to not tint the texture
                self.node.setColor(1, 1, 1, 1)
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

        # Set initial position
        self.update_visual_position()

    def update_visual_position(self):
        """Update the 3D model position to match physics position"""
        if self.node:
            # Panda3D uses different coordinate convention
            # Our physics: x=right, y=forward, z=up
            # Panda3D: x=right, y=forward, z=up (same, actually)
            self.node.setPos(self.position[0], self.position[1], self.position[2])

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

        if distance < 0.1:
            distance = 0.1  # Avoid division by zero

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
