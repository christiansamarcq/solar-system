# 3D Solar System Simulator

A real-time 3D visualization of a solar system built with Panda3D, featuring switchable physics modes, interactive controls, and educational visualizations.

## Features

### Current (Phase 1)
- Real-time 3D rendering of sun and planets
- Interactive camera controls (orbit, zoom)
- Simple orbital mechanics (circular orbits)
- Time scale control (speed up/slow down)
- Pause functionality

### Planned
- [ ] N-body gravitational physics simulation
- [ ] Orbital trail visualization
- [ ] Parameter adjustment UI
- [ ] Energy flux visualization from the sun
- [ ] Elliptical orbit support

## Installation

1. Install dependencies:
```bash
pip install -r requirements.txt
```

2. Run the simulator:
```bash
python main.py
```

## Controls

### Camera
- **Mouse Drag**: Rotate camera around the solar system
- **Mouse Wheel Up**: Zoom in
- **Mouse Wheel Down**: Zoom out

### Simulation
- **Space**: Pause/Unpause simulation
- **+/=**: Increase time scale (speed up)
- **-**: Decrease time scale (slow down)
- **ESC**: Quit

## Project Structure

```
solar-system-sim/
├── main.py                 # Entry point, main simulation loop
├── celestial_body.py      # CelestialBody class definition
├── requirements.txt        # Python dependencies
├── physics/               # Physics simulation modules
│   ├── simple_orbit.py   # Simple orbital mechanics (future)
│   └── nbody.py          # N-body simulation (future)
├── rendering/             # Visualization modules
│   ├── trail_renderer.py # Orbital trails (future)
│   └── flux_shader.py    # Energy visualization (future)
└── ui/                    # User interface
    └── controls.py        # Parameter controls (future)
```

## Learning Goals

This project demonstrates:
- 3D graphics programming with Panda3D
- Physics simulation (orbital mechanics, gravity)
- Real-time interactive applications
- Scientific visualization
- Object-oriented programming in Python

## Physics Modes

### Simple Mode (Current)
- Planets follow predetermined circular orbits
- Orbital period and distance are fixed
- Stable and predictable
- Good for understanding basic orbital concepts

### N-Body Mode (Planned)
- Realistic gravitational interactions between all bodies
- Planets affect each other's orbits
- More complex but scientifically accurate
- Demonstrates chaotic systems

## Requirements

- Python 3.7+
- Panda3D 1.10.13+
- NumPy 1.24.0+

## Development Phases

1. **Phase 1: Foundation** ✓ (Current)
   - Basic 3D scene and camera
   - Simple orbital motion
   - User controls

2. **Phase 2: Simple Physics** (Next)
   - Elliptical orbits
   - Kepler's laws implementation

3. **Phase 3: N-Body Physics**
   - Gravitational force calculations
   - Numerical integration
   - Mode switching

4. **Phase 4: Visualization**
   - Orbital trails
   - Trail fade effects

5. **Phase 5: UI Controls**
   - Parameter adjustment sliders
   - Real-time physics tuning

6. **Phase 6: Energy Visualization**
   - Solar flux calculations
   - Heat map on planets
   - Sun glow effects

## License

Educational project - Free to use and modify

## Author

Created as a learning project to explore physics simulation and 3D graphics
