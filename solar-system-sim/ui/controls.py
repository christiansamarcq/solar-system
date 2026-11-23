"""
UI controls using Panda3D's DirectGUI
Provides sliders and buttons for parameter adjustment
"""

from direct.gui.DirectGui import DirectSlider, DirectButton, DirectLabel, DirectFrame, DGG


def format_time_rate(time_scale):
    """Describe a time scale as simulated time per real second.

    The simulation advances time_scale days for every real second
    (1x = 1 day/s, so Earth orbits the Sun in about 6 minutes).
    """
    def amount(value, unit):
        text = f"{value:.3g}"
        return f"Speed: {text} {unit if text == '1' else unit + 's'}/s"

    days = time_scale
    if days < 1.0:
        return amount(days * 24, "hour")
    if days < 365.25:
        return amount(days, "day")
    return amount(days / 365.25, "year")
from panda3d.core import TextNode


class ControlPanel:
    """Main control panel for simulation parameters"""

    def __init__(self, app):
        """
        Initialize control panel

        Args:
            app: SolarSystemApp instance
        """
        self.app = app
        self.visible = True

        # Create background panel (expanded to fit size slider and new buttons)
        self.panel = DirectFrame(
            frameColor=(0.1, 0.1, 0.1, 0.8),
            frameSize=(-0.3, 0.3, -0.8, 0.5),
            pos=(0.85, 0, 0)
        )

        # Title
        self.title = DirectLabel(
            text="CONTROLS",
            scale=0.06,
            pos=(0, 0, 0.45),
            parent=self.panel,
            text_fg=(1, 1, 1, 1),
            frameColor=(0, 0, 0, 0)
        )

        # Time scale slider
        self.time_label = DirectLabel(
            text=format_time_rate(1.0),
            scale=0.045,
            pos=(0, 0, 0.35),
            parent=self.panel,
            text_fg=(1, 1, 1, 1),
            frameColor=(0, 0, 0, 0)
        )

        self.time_slider = DirectSlider(
            range=(0.1, 1000.0),
            value=1.0,
            pageSize=1.0,
            command=self.on_time_scale_changed,
            pos=(0, 0, 0.28),
            parent=self.panel,
            scale=0.5
        )

        # Size scale slider
        self.size_label = DirectLabel(
            text="Planet Size: 1.0x",
            scale=0.045,
            pos=(0, 0, 0.18),
            parent=self.panel,
            text_fg=(1, 1, 1, 1),
            frameColor=(0, 0, 0, 0)
        )

        self.size_slider = DirectSlider(
            range=(1.0, 10000.0),
            value=1.0,
            pageSize=100.0,
            command=self.on_size_scale_changed,
            pos=(0, 0, 0.11),
            parent=self.panel,
            scale=0.5
        )

        # Sun size slider
        self.sun_rate_label = DirectLabel(
            text="Sun Size: 1.0x",
            scale=0.045,
            pos=(0, 0, 0.01),
            parent=self.panel,
            text_fg=(1, 1, 1, 1),
            frameColor=(0, 0, 0, 0)
        )

        self.sun_rate_slider = DirectSlider(
            range=(0.1, 100.0),
            value=1.0,
            pageSize=0.5,
            command=self.on_sun_rate_changed,
            pos=(0, 0, -0.06),
            parent=self.panel,
            scale=0.5
        )

        # Physics mode toggle button
        self.physics_button = DirectButton(
            text="Mode: SIMPLE",
            scale=0.055,
            pos=(0, 0, -0.13),
            parent=self.panel,
            command=self.on_physics_toggle,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.3, 0.5, 0.7, 1)
        )

        # Trails toggle button
        self.trails_button = DirectButton(
            text="Trails: ON",
            scale=0.055,
            pos=(0, 0, -0.23),
            parent=self.panel,
            command=self.on_trails_toggle,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.3, 0.7, 0.5, 1)
        )

        # Clear trails button
        self.clear_button = DirectButton(
            text="Clear Trails",
            scale=0.055,
            pos=(0, 0, -0.33),
            parent=self.panel,
            command=self.on_clear_trails,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.7, 0.5, 0.3, 1)
        )

        # Pause button
        self.pause_button = DirectButton(
            text="Pause",
            scale=0.055,
            pos=(0, 0, -0.43),
            parent=self.panel,
            command=self.on_pause_toggle,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.7, 0.3, 0.3, 1)
        )

        # Collision toggle button
        self.collision_button = DirectButton(
            text="Collisions: OFF",
            scale=0.055,
            pos=(0, 0, -0.51),
            parent=self.panel,
            command=self.on_collision_toggle,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.8, 0.4, 0.2, 1)
        )

        # Reset button
        self.reset_button = DirectButton(
            text="Reset Simulation",
            scale=0.055,
            pos=(0, 0, -0.59),
            parent=self.panel,
            command=self.on_reset,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.5, 0.3, 0.7, 1)
        )

        # Zoom buttons: hold - / + to zoom (handy on a touchpad), Fit frames the body
        self.zoom_out_button = self._make_hold_button("-", -0.17, 'out')
        self.fit_button = DirectButton(
            text="Fit",
            scale=0.055,
            pos=(0, 0, -0.68),
            parent=self.panel,
            command=self.app.frame_tracked_body,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.3, 0.5, 0.6, 1)
        )
        self.zoom_in_button = self._make_hold_button("+", 0.17, 'in')

        # Hide panel button (small, moved to bottom)
        self.hide_button = DirectButton(
            text="Hide UI",
            scale=0.04,
            pos=(0, 0, -0.76),
            parent=self.panel,
            command=self.toggle_visibility,
            text_fg=(0.7, 0.7, 0.7, 1),
            frameColor=(0.2, 0.2, 0.2, 0.8)
        )

    def on_time_scale_changed(self):
        """Handle time scale slider change"""
        value = self.time_slider['value']
        self.app.time_scale = value
        self.time_label['text'] = format_time_rate(value)

    def on_size_scale_changed(self):
        """Handle size scale slider change"""
        value = self.size_slider['value']
        self.app.update_body_sizes(value)
        self.size_label['text'] = f"Planet Size: {value:.0f}x"

    def on_sun_rate_changed(self):
        """Handle sun growth slider - scales the sun independently"""
        value = self.sun_rate_slider['value']
        self.sun_rate_label['text'] = f"Sun Size: {value:.1f}x"
        if self.app.sun and self.app.sun.node:
            self.app.sun.node.setScale(max(0.1, value))

    def on_physics_toggle(self):
        """Toggle physics mode"""
        self.app.toggle_physics_mode()

        # Update button text
        mode_text = "Mode: " + self.app.physics_mode
        self.physics_button['text'] = mode_text

    def on_trails_toggle(self):
        """Toggle trails"""
        self.app.toggle_trails()

        # Update button text
        if self.app.trail_manager and self.app.trail_manager.enabled:
            self.trails_button['text'] = "Trails: ON"
        else:
            self.trails_button['text'] = "Trails: OFF"

    def on_clear_trails(self):
        """Clear all trails"""
        self.app.clear_trails()

    def on_pause_toggle(self):
        """Toggle pause"""
        self.app.toggle_pause()

        # Update button text
        if self.app.paused:
            self.pause_button['text'] = "Resume"
        else:
            self.pause_button['text'] = "Pause"

    def on_collision_toggle(self):
        """Toggle collision detection"""
        self.app.toggle_collisions()

        # Update button text
        if self.app.collisions_enabled:
            self.collision_button['text'] = "Collisions: ON"
            self.collision_button['frameColor'] = (0.8, 0.2, 0.2, 1)  # Red when ON
        else:
            self.collision_button['text'] = "Collisions: OFF"
            self.collision_button['frameColor'] = (0.8, 0.4, 0.2, 1)  # Orange when OFF

    def on_reset(self):
        """Reset the simulation"""
        self.app.reset_simulation()

    def _make_hold_button(self, label, x, direction):
        """Button that zooms continuously while the mouse button is held on it"""
        button = DirectButton(
            text=label,
            scale=0.065,
            pad=(0.5, 0.15),
            pos=(x, 0, -0.68),
            parent=self.panel,
            text_fg=(1, 1, 1, 1),
            frameColor=(0.3, 0.5, 0.6, 1)
        )
        button.bind(DGG.B1PRESS, lambda _evt: self.app.zoom_keys.__setitem__(direction, True))
        button.bind(DGG.B1RELEASE, lambda _evt: self.app.zoom_keys.__setitem__(direction, False))
        # Also stop if the pointer slides off the button while held
        button.bind(DGG.EXIT, lambda _evt: self.app.zoom_keys.__setitem__(direction, False))
        return button

    def toggle_visibility(self):
        """Show/hide the control panel"""
        if self.visible:
            self.panel.hide()
            self.visible = False
        else:
            self.panel.show()
            self.visible = True

    def destroy(self):
        """Clean up UI elements"""
        self.panel.destroy()