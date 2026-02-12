"""
Orbital trail rendering
Visualizes the path taken by celestial bodies
"""

from panda3d.core import LineSegs, NodePath
from collections import deque
import numpy as np


class TrailRenderer:
    """Renders orbital trails for a celestial body"""

    def __init__(self, body, max_points=1000, color=None, thickness=2.0):
        """
        Initialize trail renderer

        Args:
            body: CelestialBody to track
            max_points (int): Maximum number of trail points to store
            color (tuple): RGB color (0-1), defaults to body's color
            thickness (float): Line thickness
        """
        self.body = body
        self.max_points = max_points
        self.color = color if color is not None else body.color
        self.thickness = thickness

        # Ring buffer of positions
        self.positions = deque(maxlen=max_points)

        # Panda3D line segments
        self.line_segs = None
        self.line_node = None

        # Control
        self.enabled = True
        self.update_counter = 0
        self.update_interval = 3  # Add point every 3 frames to reduce geometry rebuilds
        self._dirty = False

    def update(self, render):
        """
        Update trail with current body position

        Args:
            render: Panda3D render node
        """
        if not self.enabled:
            return

        self.update_counter += 1

        # Only add a new point every N frames
        if self.update_counter % self.update_interval != 0:
            return

        # Add current position to trail
        self.positions.append(self.body.position.copy())
        self._dirty = True

        # Only rebuild geometry when new points were added
        if self._dirty and len(self.positions) > 1:
            self._rebuild_trail(render)
            self._dirty = False

    def _rebuild_trail(self, render):
        """
        Rebuild the trail geometry

        Args:
            render: Panda3D render node
        """
        # Remove old trail
        if self.line_node is not None:
            self.line_node.removeNode()

        # Create new line segments
        self.line_segs = LineSegs()
        self.line_segs.setThickness(self.thickness)

        # Set color with fade effect (older points are more transparent)
        num_points = len(self.positions)

        # Draw line segments
        for i in range(num_points - 1):
            # Fade factor (0 at oldest, 1 at newest)
            fade = (i + 1) / num_points

            # Set color with alpha fade
            self.line_segs.setColor(
                self.color[0],
                self.color[1],
                self.color[2],
                fade  # Full opacity for maximum visibility
            )

            # Draw line segment
            p1 = self.positions[i]
            p2 = self.positions[i + 1]

            self.line_segs.moveTo(p1[0], p1[1], p1[2])
            self.line_segs.drawTo(p2[0], p2[1], p2[2])

        # Create node and attach to scene
        self.line_node = render.attachNewNode(self.line_segs.create())
        # Make trails unaffected by lighting for maximum brightness
        self.line_node.setLightOff()
        # Increase brightness with color scale for better visibility
        self.line_node.setColorScale(1.5, 1.5, 1.5, 1)

    def clear(self):
        """Clear all trail points"""
        self.positions.clear()
        if self.line_node is not None:
            self.line_node.removeNode()
            self.line_node = None

    def set_enabled(self, enabled):
        """Enable or disable trail rendering"""
        self.enabled = enabled
        if not enabled and self.line_node is not None:
            self.line_node.removeNode()
            self.line_node = None

    def set_max_points(self, max_points):
        """Change maximum number of trail points"""
        self.max_points = max_points
        # Create new deque with new max length
        old_positions = list(self.positions)
        self.positions = deque(old_positions, maxlen=max_points)


class TrailManager:
    """Manages trails for all celestial bodies"""

    def __init__(self, render):
        """
        Initialize trail manager

        Args:
            render: Panda3D render node
        """
        self.render = render
        self.trails = {}
        self.enabled = True

    def add_trail(self, body, max_points=1000, color=None, thickness=2.0):
        """
        Add a trail for a celestial body

        Args:
            body: CelestialBody to track
            max_points (int): Maximum trail points
            color (tuple): Trail color (RGB)
            thickness (float): Line thickness

        Returns:
            TrailRenderer: The created trail renderer
        """
        trail = TrailRenderer(body, max_points, color, thickness)
        self.trails[body.name] = trail
        return trail

    def remove_trail(self, body):
        """Remove trail for a celestial body"""
        if body.name in self.trails:
            self.trails[body.name].clear()
            del self.trails[body.name]

    def update(self):
        """Update all trails"""
        if not self.enabled:
            return

        for trail in self.trails.values():
            trail.update(self.render)

    def clear_all(self):
        """Clear all trails"""
        for trail in self.trails.values():
            trail.clear()

    def set_enabled(self, enabled):
        """Enable or disable all trails"""
        self.enabled = enabled
        for trail in self.trails.values():
            trail.set_enabled(enabled)

    def toggle(self):
        """Toggle trails on/off"""
        self.set_enabled(not self.enabled)
