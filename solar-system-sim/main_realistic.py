"""
Realistic 3D Solar System Simulator with NASA data
Uses real planetary masses, radii, orbital parameters, and textures
"""

from direct.showbase.ShowBase import ShowBase
from panda3d.core import AmbientLight, PointLight, Texture
from panda3d.core import LVector3, GeomNode, CardMaker
from direct.gui.DirectGui import DirectOptionMenu
import numpy as np
import math
from celestial_body import CelestialBody, create_sphere_geometry
from physics import nbody
from rendering.trail_renderer import TrailManager
from ui.controls import ControlPanel, format_time_rate
import realistic_data as rd


class RealisticSolarSystemApp(ShowBase):
    """Main application class for realistic solar system simulator"""

    def __init__(self):
        ShowBase.__init__(self)

        # Simulation parameters
        self.time = 0.0
        self.time_scale = 1.0  # Simulated days per real second (1.0 = 1 day/s, Earth year ~6 min)
        self.paused = False
        self.physics_mode = "SIMPLE"  # "SIMPLE" or "NBODY"
        self.energy_vis_enabled = False  # Energy flux visualization
        self.collisions_enabled = False  # Enable/disable collision detection (default: off)
        self.size_scale = 1.0  # Visual size multiplier (1.0 = realistic, higher = exaggerated)
        self.sun_scale_rate = 0.005  # Sun scaling rate (0.5% by default)

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

        # Floating-origin root: the sun's light and explosions hang off
        # this node, which is shifted by -origin each frame (see _apply_render_origin)
        self.world_root = self.render.attachNewNode('world_root')
        self.starfield = None

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
        print("  Mouse wheel         - Zoom in/out (hold Shift for fine steps)")
        print("  Touchpad            - Two-finger scroll or pinch to zoom")
        print("  Zoom buttons        - Hold - / + in the control panel")
        print("  Up / Down arrows    - Hold to zoom in/out smoothly")
        print("  F                   - Frame tracked body (zoom to fit)")
        print("  Planet menu (top left) - Select body to track (Sun default)")
        print()
        print("SIMULATION CONTROLS:")
        print("  Space               - Pause/Unpause")
        print("  + / -               - Speed up / slow down time (starts at 1 day per second)")
        print("  M                   - Toggle physics mode (SIMPLE/N-BODY)")
        print("  T                   - Toggle orbital trails")
        print("  C                   - Clear trails")
        print("  E                   - Toggle energy flux visualization")
        print("  K                   - Toggle collision detection")
        print("  R                   - Restart simulation (reset to initial state)")
        print("  ESC                 - Quit")
        print()
        print("BODIES:")
        print(f"  - Sun + {len(rd.PLANETS)} planets + {len(rd.MOONS)} moons")
        print("=" * 70)

    def setup_camera(self):
        """Configure camera position and controls"""
        self.disableMouse()

        # Adjust camera clipping planes to allow very close and very far viewing
        # Near plane: how close objects can be before being clipped (default ~1.0)
        # Far plane: how far objects can be before being clipped (default ~1000)
        # Balance between close zoom and depth buffer precision (avoid z-fighting)
        # Using more conservative ratio to prevent flickering on small objects
        self.camLens.setNear(0.1)     # Prevent z-fighting (keep near:far ratio under 1:10M)
        self.camLens.setFar(1000000)  # Allow seeing whole solar system when zoomed far out

        # Start far out to see whole solar system
        self.camera_distance = 500
        # Zoom eases toward this distance each frame (see _update_zoom)
        self.target_camera_distance = 500
        self.zoom_keys = {'in': False, 'out': False}
        self.camera_angle_h = 0
        self.camera_angle_v = 30

        self.update_camera_position()

        self.last_mouse_x = 0
        self.last_mouse_y = 0
        self.mouse_dragging = False

    def update_camera_position(self):
        """Update camera position based on angles and distance"""
        self._apply_render_origin()
        # Target is the render origin, so it sits at (0, 0, 0) in scene space
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

        # Dynamically adjust near/far clip planes based on zoom distance
        # Keeps ratio under ~100,000:1 to prevent z-fighting
        # Near plane follows the tracked body's surface so even tiny bodies can
        # be viewed up close. The depth buffer breaks down (bodies vanish) once
        # far/near exceeds ~1e7, so the far plane is capped to that ratio.
        surface_gap = max(self.camera_distance - self._tracked_visual_radius(), 0.0)
        near = max(1e-7, 0.05 * surface_gap)
        far = min(max(1000, self.camera_distance * 100), near * 1e7)
        self.camLens.setNear(near)
        self.camLens.setFar(far)

        # Keep the starfield centred on the camera and inside the far plane,
        # so it never gets clipped away when zoomed in
        if self.starfield is not None:
            self.starfield.setPos(self.camera.getPos())
            self.starfield.setScale(far * 0.5 / 50000.0)

    def _apply_render_origin(self):
        """Re-centre the scene on the tracked body (floating origin)"""
        if self.tracked_body is not None:
            origin = np.array(self.tracked_body.position, dtype=np.float64)
        else:
            origin = np.zeros(3)
        CelestialBody.render_origin = origin
        for body in self.bodies:
            body.update_visual_position()
        self.world_root.setPos(-origin[0], -origin[1], -origin[2])
        if self.trail_manager:
            self.trail_manager.apply_origin(origin)

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
        plnp = self.world_root.attachNewNode(plight)
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
        self.sun.rotation_period = sun_data.get('rotation_period', 0)
        self.sun.axial_tilt = sun_data.get('axial_tilt', 0)
        self.sun.spin_angle = np.random.uniform(0, 360)
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

            # Set rotation properties
            planet.rotation_period = planet_data.get('rotation_period', 0)
            planet.axial_tilt = planet_data.get('axial_tilt', 0)
            planet.brightness = planet_data.get('brightness', 2.0)
            planet.spin_angle = np.random.uniform(0, 360)

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

            # Set rotation properties
            moon.rotation_period = moon_data.get('rotation_period', 0)
            moon.axial_tilt = moon_data.get('axial_tilt', 0)
            moon.spin_angle = np.random.uniform(0, 360)

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

        # Create asteroid belt
        asteroids = rd.generate_asteroid_belt(count=100)
        for ast_index, ast_data in enumerate(asteroids):
            # Use a visible minimum radius (real asteroids are too tiny to see)
            radius = max(rd.get_scaled_radius(ast_data['radius']), 0.002)
            semi_major_axis = rd.get_scaled_distance(ast_data['semi_major_axis'])
            period = rd.get_orbital_period_seconds(ast_data['orbital_period'])

            asteroid = CelestialBody(
                name=ast_data['name'],
                mass=ast_data['mass'],
                radius=radius,
                color=ast_data['color'],
                texture_path=ast_data['texture']
            )

            asteroid.rotation_period = ast_data['rotation_period']
            asteroid.axial_tilt = ast_data['axial_tilt']
            asteroid.spin_angle = np.random.uniform(0, 360)

            asteroid.set_elliptical_orbit(
                semi_major_axis=semi_major_axis,
                eccentricity=ast_data['eccentricity'],
                period=period,
                inclination=ast_data['inclination'],
                arg_periapsis=ast_data['arg_periapsis'],
                phase=ast_data['phase']
            )

            # Irregular rocky shape, different for each asteroid
            asteroid.shape_seed = ast_index
            asteroid.create_visual(self)
            self.bodies.append(asteroid)
            self.body_dict[asteroid.name] = asteroid

        print(f"Created {len(asteroids)} asteroids in the asteroid belt")

    def setup_trails(self):
        """Initialize trail rendering for all bodies"""
        self.trail_manager = TrailManager(self.render)
        self._add_body_trails()

    def _add_body_trails(self):
        """Add trail tracking for all non-sun, non-asteroid bodies"""
        for body in self.bodies:
            if body != self.sun and not body.name.startswith('Asteroid'):
                max_points = 2000 if body.name in ['Jupiter', 'Saturn', 'Uranus', 'Neptune'] else 1000
                self.trail_manager.add_trail(body, max_points=max_points, color=(1, 1, 1), thickness=4.0)

        # Add trails for a few representative asteroids to show the belt region
        asteroid_bodies = [b for b in self.bodies if b.name.startswith('Asteroid')]
        for ast in asteroid_bodies[::20]:  # every 20th asteroid = 5 trails
            self.trail_manager.add_trail(ast, max_points=1500, color=(0.7, 0.6, 0.4), thickness=2.0)

    def setup_controls(self):
        """Setup keyboard and mouse controls"""
        # Keyboard
        self.accept('space', self.toggle_pause)
        self.accept('+', self.change_time_scale, [2.0])
        self.accept('=', self.change_time_scale, [2.0])
        self.accept('-', self.change_time_scale, [0.5])
        self.accept('m', self.toggle_physics_mode)
        self.accept('t', self.toggle_trails)
        self.accept('c', self.clear_trails)
        self.accept('e', self.toggle_energy_vis)
        self.accept('k', self.toggle_collisions)
        self.accept('r', self.reset_simulation)
        self.accept('escape', exit)

        # Mouse camera controls
        self.accept('mouse1', self.start_camera_drag)
        self.accept('mouse1-up', self.stop_camera_drag)
        # Wheel / touchpad. Windows sends touchpad two-finger scroll as a fast
        # stream of small wheel events, and pinch as Ctrl + wheel, which
        # Panda3D reports as 'control-wheel_up' / 'control-wheel_down'.
        for prefix in ('', 'control-', 'alt-'):
            self.accept(prefix + 'wheel_up', self.on_wheel, [1, False])
            self.accept(prefix + 'wheel_down', self.on_wheel, [-1, False])
        self.accept('shift-wheel_up', self.on_wheel, [1, True])
        self.accept('shift-wheel_down', self.on_wheel, [-1, True])
        self._last_wheel_time = -1.0

        # Hold arrow keys for continuous zoom; F frames the tracked body
        self.accept('arrow_up', self.zoom_keys.__setitem__, ['in', True])
        self.accept('arrow_up-up', self.zoom_keys.__setitem__, ['in', False])
        self.accept('arrow_down', self.zoom_keys.__setitem__, ['out', True])
        self.accept('arrow_down-up', self.zoom_keys.__setitem__, ['out', False])
        self.accept('f', self.frame_tracked_body)

    def start_camera_drag(self):
        """Start camera rotation"""
        # Ignore clicks on the GUI (buttons, sliders) so they don't spin the camera
        if self.mouseWatcherNode.hasMouse() and not self.mouseWatcherNode.isOverRegion():
            self.mouse_dragging = True
            self.last_mouse_x = self.mouseWatcherNode.getMouseX()
            self.last_mouse_y = self.mouseWatcherNode.getMouseY()

    def stop_camera_drag(self):
        """Stop camera rotation"""
        self.mouse_dragging = False

    MAX_CAMERA_DISTANCE = 1000000

    def _tracked_visual_radius(self):
        """Radius of the tracked body as drawn (includes the size sliders)"""
        body = self.tracked_body
        if body is None:
            return 0.0
        uses_sun = getattr(body, 'uses_sun_scale', body.is_emissive)
        return body.radius * (self._get_sun_scale() if uses_sun else self.size_scale)

    def _min_camera_distance(self):
        """Closest zoom: just outside the tracked body's surface"""
        return max(1e-5, self._tracked_visual_radius() * 1.2)

    def _set_zoom_target(self, distance):
        self.target_camera_distance = min(self.MAX_CAMERA_DISTANCE,
                                          max(self._min_camera_distance(), distance))

    def on_wheel(self, direction, fine):
        """Wheel or touchpad zoom step (direction +1 = in, -1 = out).

        A mouse wheel sends one event per notch, so each moves 20%. A touchpad
        sends many events close together; those are treated as a continuous
        gesture with small 4% steps so it zooms smoothly instead of racing.
        """
        now = globalClock.getRealTime()
        streaming = (now - self._last_wheel_time) < 0.06
        self._last_wheel_time = now
        if fine:
            factor = 0.95
        elif streaming:
            factor = 0.96
        else:
            factor = 0.8
        if direction > 0:
            self.zoom_in(factor)
        else:
            self.zoom_out(factor)

    def zoom_in(self, factor=0.8):
        """Zoom camera in (each wheel step moves 20% closer)"""
        self._set_zoom_target(self.target_camera_distance * factor)

    def zoom_out(self, factor=0.8):
        """Zoom camera out"""
        self._set_zoom_target(self.target_camera_distance / factor)

    def frame_tracked_body(self):
        """Zoom so the tracked body fills a comfortable part of the screen"""
        if self.tracked_body == self.sun:
            self._set_zoom_target(max(50, self._tracked_visual_radius() * 5))
        else:
            self._set_zoom_target(self._tracked_visual_radius() * 4)

    def _update_zoom(self, dt):
        """Hold-to-zoom keys, plus smooth easing toward the target distance"""
        dt = min(dt, 0.1)
        # If the body grew (size sliders), back the camera out of it
        if self.target_camera_distance < self._min_camera_distance():
            self._set_zoom_target(self.target_camera_distance)
        if self.zoom_keys['in']:
            self.zoom_in(0.33 ** dt)   # ~3x closer per second
        if self.zoom_keys['out']:
            self.zoom_out(0.33 ** dt)
        # Ease in log space so zooming feels the same at every scale
        current = math.log(self.camera_distance)
        target = math.log(self.target_camera_distance)
        blend = 1.0 - math.exp(-12.0 * dt)
        self.camera_distance = math.exp(current + (target - current) * blend)

    def create_planet_menu(self):
        """Create dropdown menu for planet selection"""
        # Create list of body names (exclude asteroids except first one as representative)
        body_names = [body.name for body in self.bodies if not body.name.startswith('Asteroid')]
        asteroid_bodies = [body for body in self.bodies if body.name.startswith('Asteroid')]
        if asteroid_bodies:
            body_names.append(asteroid_bodies[0].name)

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

            # Auto-zoom so the body is nicely framed, however small it is
            self.frame_tracked_body()
            self.camera_distance = self.target_camera_distance
            self.update_camera_position()
            print(f"Now tracking: {selected_name} (distance: {self.camera_distance:.3g} units)")

    def change_time_scale(self, factor):
        """Speed time up or down (keyboard + / -), keeping the slider and label in sync"""
        self.time_scale = min(1000.0, max(0.01, self.time_scale * factor))
        if self.control_panel:
            # Moving the slider runs on_time_scale_changed, which updates the label
            self.control_panel.time_slider['value'] = self.time_scale
        print(format_time_rate(self.time_scale))

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
            # First initialize planets (they orbit the sun)
            for body in self.bodies:
                if body != self.sun and not hasattr(body, 'parent'):
                    nbody.initialize_circular_orbit_velocity(body, self.sun.mass, self.sun.position)
            # Then initialize moons - add parent's velocity so they move WITH their planet
            for body in self.bodies:
                if hasattr(body, 'parent'):
                    nbody.initialize_circular_orbit_velocity(body, body.parent.mass, body.parent.position)
                    body.velocity += body.parent.velocity  # Add parent's orbital velocity
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

    def reset_simulation(self):
        """Reset simulation to initial state"""
        print("\n" + "=" * 70)
        print("RESETTING SIMULATION...")
        print("=" * 70)

        # Remember trail state before reset
        trails_were_enabled = self.trail_manager and self.trail_manager.enabled

        # Remove all existing bodies
        for body in self.bodies[:]:  # Use slice to avoid modifying list during iteration
            if body.node:
                body.node.removeNode()

        # Clear lists and dictionaries
        self.bodies.clear()
        self.body_dict.clear()

        # Reset simulation parameters
        self.time = 0.0
        self.time_scale = 1.0
        self.paused = False
        self.physics_mode = "SIMPLE"
        self.size_scale = 1.0
        self.collisions_enabled = False
        self.sun_scale_rate = 0.005

        # Clear trail history but preserve enabled state
        self.clear_trails()

        # Recreate solar system
        self.create_realistic_solar_system()

        # Re-add trails for new bodies if trails were enabled
        if trails_were_enabled and self.trail_manager:
            self.trail_manager.enabled = True
            self._add_body_trails()

        # Reset camera to track sun
        self.tracked_body = self.sun
        self.camera_distance = 50.0
        self.target_camera_distance = 50.0
        self.update_camera_position()

        # Update control panel sliders to default values
        if self.control_panel:
            self.control_panel.time_slider['value'] = 1.0  # also updates the speed label
            self.control_panel.size_slider['value'] = 1.0
            self.control_panel.size_label['text'] = "Planet Size: 1.0x"
            self.control_panel.sun_rate_slider['value'] = 1.0
            self.control_panel.sun_rate_label['text'] = "Sun Size: 1.0x"
            self.control_panel.physics_button['text'] = "Mode: SIMPLE"
            self.control_panel.pause_button['text'] = "Pause"
            self.control_panel.collision_button['text'] = "Collisions: OFF"
            self.control_panel.collision_button['frameColor'] = (0.8, 0.4, 0.2, 1)
            # Update trails button to reflect current state
            if trails_were_enabled:
                self.control_panel.trails_button['text'] = "Trails: ON"
            else:
                self.control_panel.trails_button['text'] = "Trails: OFF"

        # Recreate planet menu with new bodies
        if self.planet_menu:
            self.planet_menu.destroy()
        self.create_planet_menu()

        print("Simulation reset complete!")
        print("=" * 70 + "\n")

    def handle_collisions(self, collisions):
        """Handle collision events between celestial bodies"""
        bodies_to_remove = set()
        bodies_to_add = []

        for body1, body2 in collisions:
            if body1 in bodies_to_remove or body2 in bodies_to_remove:
                continue

            print(f"\n{'='*60}")
            print(f"COLLISION EVENT: {body1.name} + {body2.name}")
            print(f"{'='*60}")

            # Create merged body from collision
            merged_body = self._create_merged_body(body1, body2)
            merged_body.create_visual(self, segments=64, rings=32)

            # Apply the winner's scale type
            if merged_body.node:
                if merged_body.uses_sun_scale:
                    merged_body.node.setScale(self.control_panel.sun_rate_slider['value'] if self.control_panel else 1.0)
                else:
                    merged_body.node.setScale(self.size_scale)

            bodies_to_remove.add(body1)
            bodies_to_remove.add(body2)
            bodies_to_add.append(merged_body)

        # Remove collided bodies and add merged ones
        for body in bodies_to_remove:
            if body.node:
                body.node.removeNode()
            if body in self.bodies:
                self.bodies.remove(body)
            if body.name in self.body_dict:
                del self.body_dict[body.name]
            # Remove trail
            if self.trail_manager:
                self.trail_manager.remove_trail(body)
            print(f"  Removed: {body.name}")

        for body in bodies_to_add:
            self.bodies.append(body)
            self.body_dict[body.name] = body
            # If merged body is emissive (absorbed the sun), update self.sun
            if body.is_emissive:
                self.sun = body
            else:
                # Only add trail for non-sun bodies
                if self.trail_manager and self.trail_manager.enabled:
                    self.trail_manager.add_trail(body, max_points=1000, color=(1, 1, 1), thickness=4.0)
            print(f"  Added: {body.name}")

        # Refresh planet menu if any collisions happened
        if bodies_to_remove:
            if self.planet_menu:
                self.planet_menu.destroy()
            self.create_planet_menu()

    def _create_merged_body(self, body1, body2):
        """Create a new body from two colliding bodies using conservation laws"""
        total_mass = body1.mass + body2.mass
        merged_velocity = (body1.velocity * body1.mass + body2.velocity * body2.mass) / total_mass
        merged_position = (body1.position * body1.mass + body2.position * body2.mass) / total_mass
        # Additive radius: larger body + half of smaller body for visible growth
        merged_radius = max(body1.radius, body2.radius) + min(body1.radius, body2.radius) * 0.5

        # Larger body (by visual size) determines name and texture
        visual1 = body1.radius * (self._get_sun_scale() if body1.is_emissive else self.size_scale)
        visual2 = body2.radius * (self._get_sun_scale() if body2.is_emissive else self.size_scale)
        larger, smaller = (body1, body2) if visual1 >= visual2 else (body2, body1)
        mass_ratio = larger.mass / total_mass
        merged_color = tuple(
            larger.color[i] * mass_ratio + smaller.color[i] * (1 - mass_ratio)
            for i in range(3)
        )

        print(f"  Merged body: {larger.name}+{smaller.name}")
        print(f"  Total mass: {total_mass:.3e} kg")
        print(f"  Merged radius: {merged_radius:.3f} units")

        merged = CelestialBody(
            name=f"{larger.name}+{smaller.name}",
            mass=total_mass,
            radius=merged_radius,
            color=merged_color,
            position=merged_position.copy(),
            velocity=merged_velocity.copy(),
            is_emissive=larger.is_emissive or smaller.is_emissive,
            texture_path=larger.texture_path
        )
        # Use the winner's scale type (planet ate sun = use planet scale)
        merged.uses_sun_scale = larger.uses_sun_scale if hasattr(larger, 'uses_sun_scale') else larger.is_emissive
        return merged

    def create_explosion(self, position, radius):
        """
        Create visual explosion effect at given position

        Args:
            position (np.array): Position of explosion in world coordinates
            radius (float): Size of explosion effect
        """
        from panda3d.core import Point3, Vec4
        from direct.interval.IntervalGlobal import Sequence, Parallel, LerpScaleInterval, LerpColorScaleInterval, Func

        # Create explosion sphere
        explosion_node = self.world_root.attachNewNode("explosion")
        explosion_node.setPos(Point3(position[0], position[1], position[2]))

        # Create sphere geometry for explosion
        from panda3d.core import CardMaker
        card = CardMaker("explosion_sprite")
        card.setFrame(-radius, radius, -radius, radius)
        explosion_sprite = explosion_node.attachNewNode(card.generate())
        explosion_sprite.setBillboardPointEye()

        # Set initial appearance (bright yellow/orange)
        explosion_sprite.setColor(Vec4(1.0, 0.8, 0.2, 1.0))
        explosion_sprite.setScale(0.1)

        # Animate explosion: expand and fade out
        explosion_sequence = Sequence(
            Parallel(
                LerpScaleInterval(explosion_sprite, 1.5, radius * 3, startScale=0.1),
                LerpColorScaleInterval(explosion_sprite, 1.5, Vec4(1.0, 0.3, 0.0, 0.0), startColorScale=Vec4(1.0, 0.8, 0.2, 1.0))
            ),
            Func(explosion_node.removeNode)
        )
        explosion_sequence.start()

        print(f"  EXPLOSION created at position {position} with radius {radius}")

    def toggle_energy_vis(self):
        """Toggle energy flux visualization"""
        self.energy_vis_enabled = not self.energy_vis_enabled
        print(f"Energy flux visualization: {'ON' if self.energy_vis_enabled else 'OFF'}")

    def toggle_collisions(self):
        """Toggle collision detection"""
        self.collisions_enabled = not self.collisions_enabled
        print(f"Collision detection: {'ON' if self.collisions_enabled else 'OFF'}")

    def _get_sun_scale(self):
        """Get the current sun visual scale from the slider"""
        if self.control_panel:
            return self.control_panel.sun_rate_slider['value']
        return 1.0

    def update_body_sizes(self, new_scale):
        """Update visual size of all bodies"""
        self.size_scale = new_scale

        for body in self.bodies:
            uses_sun = getattr(body, 'uses_sun_scale', body.is_emissive)
            if uses_sun or not body.node:
                continue
            # Scale planets and moons
            body.node.setScale(new_scale)

            # Make planets brighter as they get bigger (easier to see)
            if not body.is_emissive:
                brightness = min(4.0, body.brightness + (new_scale / 2500.0))
                body.node.setColorScale(brightness, brightness, brightness, 1)

        print(f"Planet size scale: {new_scale:.1f}x ({'realistic' if new_scale == 1.0 else 'exaggerated'})")

    def update(self, task):
        """Main update loop"""
        dt = globalClock.getDt()

        if not self.paused:
            self._update_mouse_camera()
            self._update_physics(dt)
            self._update_energy_visualization()
            self._update_body_rotations(dt)

            if self.trail_manager:
                self.trail_manager.update()

        # Update camera to follow tracked body (even when paused, for manual rotation)
        self._update_zoom(dt)
        self.update_camera_position()

        return task.cont

    def _update_mouse_camera(self):
        """Handle mouse drag camera rotation"""
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

    def _update_physics(self, dt):
        """Step physics simulation forward"""
        if self.physics_mode == "SIMPLE":
            for body in self.bodies:
                if body != self.sun:
                    if hasattr(body, 'parent'):
                        body.update_simple_orbit(self.time)
                        body.position += body.parent.position
                        body.update_visual_position()
                    else:
                        body.update_simple_orbit(self.time)

            self.time += dt * self.time_scale * 86400.0

        elif self.physics_mode == "NBODY":
            total_dt = dt * self.time_scale * 86400.0
            max_substep = 50.0
            max_substeps = 50

            all_collisions = []
            if total_dt <= max_substep:
                all_collisions = nbody.update_nbody_physics(self.bodies, total_dt, self.size_scale, self._get_sun_scale())
            else:
                num_substeps = min(int(np.ceil(total_dt / max_substep)), max_substeps)
                substep_dt = total_dt / num_substeps
                for _ in range(num_substeps):
                    collisions = nbody.update_nbody_physics(self.bodies, substep_dt, self.size_scale, self._get_sun_scale())
                    all_collisions.extend(collisions)

            if self.collisions_enabled and all_collisions:
                self.handle_collisions(all_collisions)

            self._check_stability()

    def _check_stability(self):
        """Check for NaN values indicating numerical instability"""
        for body in self.bodies:
            if np.any(np.isnan(body.position)) or np.any(np.isnan(body.velocity)):
                print("WARNING: N-body simulation became unstable (NaN detected)")
                print("Automatically switching back to SIMPLE mode")
                self.physics_mode = "SIMPLE"
                self.time = 0.0
                for reset_body in self.bodies:
                    if hasattr(reset_body, 'orbital_params'):
                        reset_body.update_simple_orbit(0.0)
                if self.control_panel:
                    self.control_panel.physics_button['text'] = "Mode: SIMPLE"
                self.clear_trails()
                break

    def _update_body_rotations(self, dt):
        """Update axial spin for all bodies"""
        sim_dt = dt * self.time_scale * 86400.0  # convert to simulation seconds
        for body in self.bodies:
            body.update_spin(sim_dt)

    def _update_energy_visualization(self):
        """Update energy flux visualization for all bodies"""
        if self.energy_vis_enabled:
            for body in self.bodies:
                if body != self.sun:
                    body.calculate_energy_flux(self.sun.position, self.sun.mass)
                    body.update_heat_visualization(True)
        else:
            for body in self.bodies:
                body.update_heat_visualization(False)


if __name__ == "__main__":
    app = RealisticSolarSystemApp()
    app.run()
