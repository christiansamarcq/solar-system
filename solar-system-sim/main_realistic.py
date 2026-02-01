"""
Realistic 3D Solar System Simulator with NASA data
Uses real planetary masses, radii, orbital parameters, and textures
"""

from direct.showbase.ShowBase import ShowBase
from panda3d.core import AmbientLight, DirectionalLight, PointLight, Texture
from panda3d.core import LVector3, GeomNode, CardMaker
from direct.gui.DirectGui import DirectOptionMenu
import numpy as np
import math
from celestial_body import CelestialBody, create_sphere_geometry
from physics import nbody
from rendering.trail_renderer import TrailManager
from ui.controls import ControlPanel
import realistic_data as rd


class RealisticSolarSystemApp(ShowBase):
    """Main application class for realistic solar system simulator"""

    def __init__(self):
        ShowBase.__init__(self)

        # Simulation parameters
        self.time = 0.0
        self.time_scale = 1.0  # Speed multiplier (realtime = 1.0 second per second)
        self.paused = False
        self.physics_mode = "SIMPLE"  # "SIMPLE" or "NBODY"
        self.energy_vis_enabled = False  # Energy flux visualization
        self.size_scale = 1.0  # Visual size multiplier (1.0 = realistic, higher = exaggerated)

        # Celestial bodies
        self.bodies = []
        self.sun = None
        self.body_dict = {}  # Dictionary for quick lookup by name

        # Camera tracking
        self.tracked_body = None  # Body to center camera on

        # Trail manager
        self.trail_manager = None

        # UI Control panel
        self.control_panel = None
        self.planet_menu = None

        # Setup
        self.setup_camera()
        self.setup_lights()
        self.create_starfield()
        self.create_realistic_solar_system()
        self.setup_trails()
        self.setup_controls()

        # Set default tracked body (Sun)
        self.tracked_body = self.sun

        # Create UI control panel
        self.control_panel = ControlPanel(self)

        # Create planet selection menu
        self.create_planet_menu()

        # Update task
        self.taskMgr.add(self.update, "update_task")

        print("=" * 70)
        print("REALISTIC 3D SOLAR SYSTEM SIMULATOR - NASA Data")
        print("=" * 70)
        print("DATA SOURCES:")
        print("  - NASA Planetary Fact Sheet")
        print("  - JPL Solar System Dynamics")
        print("  - Textures from Solar System Scope")
        print()
        print("SCALE:")
        print(f"  Distance: 1 unit = {1.0/rd.DISTANCE_SCALE/1e9:.1f} billion meters")
        print(f"  Radius: 1 unit = {1.0/rd.RADIUS_SCALE/1e9:.1f} billion meters (REALISTIC!)")
        print()
        print("CAMERA CONTROLS:")
        print("  Mouse drag          - Rotate camera around tracked body")
        print("  Mouse wheel         - Zoom in/out")
        print("  Planet menu (top left) - Select body to track (Sun default)")
        print()
        print("SIMULATION CONTROLS:")
        print("  Space               - Pause/Unpause")
        print("  + / -               - Increase/Decrease time scale")
        print("  M                   - Toggle physics mode (SIMPLE/N-BODY)")
        print("  T                   - Toggle orbital trails")
        print("  C                   - Clear trails")
        print("  E                   - Toggle energy flux visualization")
        print("  ESC                 - Quit")
        print()
        print("BODIES:")
        print(f"  - Sun + {len(rd.PLANETS)} planets + {len(rd.MOONS)} moons")
        print("=" * 70)

    def setup_camera(self):
        """Configure camera position and controls"""
        self.disableMouse()

        # Start far out to see whole solar system
        self.camera_distance = 500
        self.camera_angle_h = 0
        self.camera_angle_v = 30

        self.update_camera_position()

        self.last_mouse_x = 0
        self.last_mouse_y = 0
        self.mouse_dragging = False

    def update_camera_position(self):
        """Update camera position based on angles and distance"""
        # Get target position (tracked body or origin)
        if self.tracked_body is not None:
            target_pos = self.tracked_body.position
        else:
            target_pos = np.array([0.0, 0.0, 0.0])

        # Convert to radians
        h_rad = math.radians(self.camera_angle_h)
        v_rad = math.radians(self.camera_angle_v)

        # Spherical to Cartesian offset from target
        x = self.camera_distance * math.cos(v_rad) * math.sin(h_rad)
        y = -self.camera_distance * math.cos(v_rad) * math.cos(h_rad)
        z = self.camera_distance * math.sin(v_rad)

        # Set camera position relative to tracked body
        self.camera.setPos(
            target_pos[0] + x,
            target_pos[1] + y,
            target_pos[2] + z
        )
        self.camera.lookAt(target_pos[0], target_pos[1], target_pos[2])

    def setup_lights(self):
        """Setup scene lighting"""
        # Very minimal ambient light for dark space
        alight = AmbientLight('alight')
        alight.setColor((0.05, 0.05, 0.05, 1))
        alnp = self.render.attachNewNode(alight)
        self.render.setLight(alnp)

        # Point light at sun
        plight = PointLight('sunlight')
        plight.setColor((1, 1, 0.9, 1))
        plnp = self.render.attachNewNode(plight)
        plnp.setPos(0, 0, 0)
        self.render.setLight(plnp)

    def create_starfield(self):
        """Create starfield background sphere"""
        # Set background color to black
        self.setBackgroundColor(0, 0, 0)

        # Create a large inverted sphere for the starfield
        # Use a radius large enough to encompass the entire solar system view
        starfield_radius = 50000

        # Create sphere geometry (inverted normals by reversing winding order)
        from panda3d.core import GeomVertexFormat, GeomVertexData, GeomVertexWriter
        from panda3d.core import Geom, GeomTriangles, GeomNode

        format = GeomVertexFormat.getV3n3t2()
        vdata = GeomVertexData('starfield', format, Geom.UHStatic)

        vertex = GeomVertexWriter(vdata, 'vertex')
        normal = GeomVertexWriter(vdata, 'normal')
        texcoord = GeomVertexWriter(vdata, 'texcoord')

        # Generate sphere vertices
        segments = 32
        rings = 16

        for ring in range(rings + 1):
            theta = ring * math.pi / rings
            sin_theta = math.sin(theta)
            cos_theta = math.cos(theta)

            for seg in range(segments + 1):
                phi = seg * 2 * math.pi / segments
                sin_phi = math.sin(phi)
                cos_phi = math.cos(phi)

                x = starfield_radius * sin_theta * cos_phi
                y = starfield_radius * sin_theta * sin_phi
                z = starfield_radius * cos_theta

                vertex.addData3(x, y, z)
                # Inverted normals (point inward)
                normal.addData3(-sin_theta * cos_phi, -sin_theta * sin_phi, -cos_theta)

                u = seg / segments
                v = ring / rings
                texcoord.addData2(u, v)

        # Create triangles (reversed winding for inverted sphere)
        tris = GeomTriangles(Geom.UHStatic)

        for ring in range(rings):
            for seg in range(segments):
                current = ring * (segments + 1) + seg
                next_seg = ring * (segments + 1) + (seg + 1)
                next_ring = (ring + 1) * (segments + 1) + seg
                next_ring_next = (ring + 1) * (segments + 1) + (seg + 1)

                # Reversed winding order
                tris.addVertices(current, next_seg, next_ring)
                tris.addVertices(next_seg, next_ring_next, next_ring)

        tris.closePrimitive()

        geom = Geom(vdata)
        geom.addPrimitive(tris)

        node = GeomNode('starfield')
        node.addGeom(geom)

        # Attach to render
        self.starfield = self.render.attachNewNode(node)

        # Load and apply starfield texture
        try:
            tex = self.loader.loadTexture('textures/stars_milky_way.jpg')
            self.starfield.setTexture(tex)
        except:
            print("Warning: Could not load starfield texture")
            # Fallback: set to dark blue/black
            self.starfield.setColor(0.01, 0.01, 0.02, 1)

        # Disable lighting on starfield
        self.starfield.setLightOff()
        # Make sure it renders behind everything
        self.starfield.setBin('background', 0)
        self.starfield.setDepthWrite(False)

    def create_realistic_solar_system(self):
        """Create solar system using NASA data"""
        # Create the Sun
        sun_data = rd.SUN
        self.sun = CelestialBody(
            name=sun_data['name'],
            mass=sun_data['mass'],
            radius=rd.get_scaled_radius(sun_data['radius']),
            color=sun_data['color'],
            position=np.array([0.0, 0.0, 0.0]),
            velocity=np.array([0.0, 0.0, 0.0]),
            is_emissive=True,
            texture_path=sun_data['texture']
        )
        self.sun.create_visual(self, segments=64, rings=32)
        self.bodies.append(self.sun)
        self.body_dict[self.sun.name] = self.sun

        # Create planets
        for planet_name, planet_data in rd.PLANETS.items():
            # Scale parameters
            mass = planet_data['mass']
            radius = rd.get_scaled_radius(planet_data['radius'])
            semi_major_axis = rd.get_scaled_distance(planet_data['semi_major_axis'])
            period = rd.get_orbital_period_seconds(planet_data['orbital_period'])

            # Create planet
            planet = CelestialBody(
                name=planet_data['name'],
                mass=mass,
                radius=radius,
                color=planet_data['color'],
                texture_path=planet_data['texture']
            )

            # Set elliptical orbit parameters
            planet.set_elliptical_orbit(
                semi_major_axis=semi_major_axis,
                eccentricity=planet_data['eccentricity'],
                period=period,
                inclination=planet_data['inclination'],
                arg_periapsis=planet_data['arg_periapsis'],
                phase=0.0
            )

            # Create visual
            planet.create_visual(self, segments=64, rings=32)
            self.bodies.append(planet)
            self.body_dict[planet.name] = planet

            print(f"Created {planet_name}: {semi_major_axis:.1f} units, period {period/(86400):.1f} days")

        # Create moons
        for moon_name, moon_data in rd.MOONS.items():
            # Find parent planet
            parent_name = moon_data['parent']
            parent = None
            for body in self.bodies:
                if body.name == parent_name:
                    parent = body
                    break

            if parent is None:
                print(f"Warning: Could not find parent {parent_name} for moon {moon_name}")
                continue

            # Scale parameters
            mass = moon_data['mass']
            radius = rd.get_scaled_radius(moon_data['radius'])
            semi_major_axis = rd.get_scaled_distance(moon_data['semi_major_axis'])
            period = rd.get_orbital_period_seconds(moon_data['orbital_period'])

            # Create moon
            moon = CelestialBody(
                name=moon_data['name'],
                mass=mass,
                radius=radius,
                color=moon_data['color'],
                texture_path=moon_data['texture']
            )

            # Set orbit around parent (will need custom update logic)
            moon.set_elliptical_orbit(
                semi_major_axis=semi_major_axis,
                eccentricity=moon_data['eccentricity'],
                period=period,
                inclination=moon_data['inclination'],
                phase=0.0
            )

            # Store parent reference
            moon.parent = parent

            # Create visual (smaller detail for moons)
            moon.create_visual(self, segments=32, rings=16)
            self.bodies.append(moon)
            self.body_dict[moon.name] = moon

            print(f"Created {moon_name} orbiting {parent_name}")

    def setup_trails(self):
        """Initialize trail rendering for all bodies"""
        self.trail_manager = TrailManager(self.render)

        # Add trails for planets (not sun) - all white for visibility on black background
        for body in self.bodies:
            if body != self.sun:
                # Longer trails for outer planets
                max_points = 2000 if body.name in ['Jupiter', 'Saturn', 'Uranus', 'Neptune'] else 1000
                # White trails for all bodies
                self.trail_manager.add_trail(body, max_points=max_points, color=(1, 1, 1), thickness=2.0)

    def setup_controls(self):
        """Setup keyboard and mouse controls"""
        # Keyboard
        self.accept('space', self.toggle_pause)
        self.accept('+', lambda: setattr(self, 'time_scale', min(1000.0, self.time_scale * 2.0)))
        self.accept('=', lambda: setattr(self, 'time_scale', min(1000.0, self.time_scale * 2.0)))
        self.accept('-', lambda: setattr(self, 'time_scale', max(0.1, self.time_scale / 2.0)))
        self.accept('m', self.toggle_physics_mode)
        self.accept('t', self.toggle_trails)
        self.accept('c', self.clear_trails)
        self.accept('e', self.toggle_energy_vis)
        self.accept('escape', exit)

        # Mouse camera controls
        self.accept('mouse1', self.start_camera_drag)
        self.accept('mouse1-up', self.stop_camera_drag)
        self.accept('wheel_up', self.zoom_in)
        self.accept('wheel_down', self.zoom_out)

    def start_camera_drag(self):
        """Start camera rotation"""
        if self.mouseWatcherNode.hasMouse():
            self.mouse_dragging = True
            self.last_mouse_x = self.mouseWatcherNode.getMouseX()
            self.last_mouse_y = self.mouseWatcherNode.getMouseY()

    def stop_camera_drag(self):
        """Stop camera rotation"""
        self.mouse_dragging = False

    def zoom_in(self):
        """Zoom camera in"""
        self.camera_distance = max(10, self.camera_distance * 0.9)
        self.update_camera_position()

    def zoom_out(self):
        """Zoom camera out"""
        self.camera_distance = min(5000, self.camera_distance * 1.1)
        self.update_camera_position()

    def create_planet_menu(self):
        """Create dropdown menu for planet selection"""
        # Create list of body names
        body_names = [body.name for body in self.bodies]

        # Create dropdown menu
        self.planet_menu = DirectOptionMenu(
            text="Track:",
            scale=0.07,
            items=body_names,
            initialitem=0,  # Start with Sun (first body)
            highlightColor=(0.65, 0.65, 0.65, 1),
            command=self.on_planet_selected,
            pos=(-1.5, 0, 0.9),
            textMayChange=1
        )

    def on_planet_selected(self, selected_name):
        """Handle planet selection from menu"""
        if selected_name in self.body_dict:
            self.tracked_body = self.body_dict[selected_name]

            # Auto-zoom based on body size
            body_radius = self.tracked_body.radius * self.size_scale

            if self.tracked_body == self.sun:
                # Zoom out to see sun comfortably
                self.camera_distance = max(50, body_radius * 5)
            else:
                # Zoom in to see planets/moons (planets are tiny at realistic scale)
                # At realistic scale, need to zoom very close
                self.camera_distance = max(1, body_radius * 10)

            self.update_camera_position()
            print(f"Now tracking: {selected_name} (distance: {self.camera_distance:.2f} units)")

    def toggle_pause(self):
        """Toggle simulation pause"""
        self.paused = not self.paused
        print(f"Simulation {'PAUSED' if self.paused else 'RUNNING'}")

    def toggle_physics_mode(self):
        """Switch between simple and N-body physics"""
        if self.physics_mode == "SIMPLE":
            self.physics_mode = "NBODY"
            # Initialize velocities for N-body
            print("Switched to N-BODY physics mode")
            print("Initializing circular orbit velocities...")
            for body in self.bodies:
                if body != self.sun and not hasattr(body, 'parent'):
                    # Planet orbiting sun
                    nbody.initialize_circular_orbit_velocity(body, self.sun.mass, self.sun.position)
                elif hasattr(body, 'parent'):
                    # Moon orbiting planet
                    nbody.initialize_circular_orbit_velocity(body, body.parent.mass, body.parent.position)
        else:
            self.physics_mode = "SIMPLE"
            print("Switched to SIMPLE physics mode")

        self.clear_trails()

    def toggle_trails(self):
        """Toggle trail rendering"""
        if self.trail_manager:
            self.trail_manager.toggle()

    def clear_trails(self):
        """Clear all trails"""
        if self.trail_manager:
            self.trail_manager.clear_all()
            print("Trails cleared")

    def toggle_energy_vis(self):
        """Toggle energy flux visualization"""
        self.energy_vis_enabled = not self.energy_vis_enabled
        print(f"Energy flux visualization: {'ON' if self.energy_vis_enabled else 'OFF'}")

    def update_body_sizes(self, new_scale):
        """Update visual size of all bodies"""
        self.size_scale = new_scale

        for body in self.bodies:
            if body.node:
                if body == self.sun:
                    # Don't scale the sun, or scale it minimally
                    body.node.setScale(1.0)
                else:
                    # Scale planets and moons
                    body.node.setScale(new_scale)

                    # Make planets brighter as they get bigger (easier to see)
                    if not body.is_emissive:
                        # Increase color intensity for visibility
                        brightness = min(2.0, 1.0 + (new_scale / 5000.0))
                        body.node.setColorScale(brightness, brightness, brightness, 1)

        print(f"Planet size scale: {new_scale:.1f}x ({'realistic' if new_scale == 1.0 else 'exaggerated'})")

    def update(self, task):
        """Main update loop"""
        dt = globalClock.getDt()

        if not self.paused:
            # Mouse camera control
            if self.mouse_dragging and self.mouseWatcherNode.hasMouse():
                mouse_x = self.mouseWatcherNode.getMouseX()
                mouse_y = self.mouseWatcherNode.getMouseY()

                dx = mouse_x - self.last_mouse_x
                dy = mouse_y - self.last_mouse_y

                self.camera_angle_h += dx * 100
                self.camera_angle_v = max(-89, min(89, self.camera_angle_v - dy * 100))

                self.update_camera_position()

                self.last_mouse_x = mouse_x
                self.last_mouse_y = mouse_y

            # Physics update
            if self.physics_mode == "SIMPLE":
                # Simple orbital mechanics
                for body in self.bodies:
                    if body != self.sun:
                        if hasattr(body, 'parent'):
                            # Moon: orbit relative to parent
                            body.update_simple_orbit(self.time)
                            # Offset by parent position
                            body.position += body.parent.position
                            body.update_visual_position()
                        else:
                            # Planet: orbit around sun
                            body.update_simple_orbit(self.time)

                self.time += dt * self.time_scale * 86400.0  # Scale to days

            elif self.physics_mode == "NBODY":
                # N-body gravitational simulation with substeps
                total_dt = dt * self.time_scale
                max_substep = 1000.0  # seconds

                if total_dt <= max_substep:
                    nbody.update_nbody_physics(self.bodies, total_dt)
                else:
                    num_substeps = int(np.ceil(total_dt / max_substep))
                    substep_dt = total_dt / num_substeps
                    for _ in range(num_substeps):
                        nbody.update_nbody_physics(self.bodies, substep_dt)

            # Update energy flux visualization
            if self.energy_vis_enabled:
                for body in self.bodies:
                    if body != self.sun:
                        body.calculate_energy_flux(self.sun.position, self.sun.mass)
                        body.update_heat_visualization(True)
            else:
                for body in self.bodies:
                    body.update_heat_visualization(False)

            # Update trails
            if self.trail_manager:
                self.trail_manager.update()

        # Update camera to follow tracked body (even when paused, for manual rotation)
        self.update_camera_position()

        return task.cont


if __name__ == "__main__":
    app = RealisticSolarSystemApp()
    app.run()
