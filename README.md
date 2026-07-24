# Quantum CFD

## Core
### Project goal
- Determine whether CFD simulations can be performed using quantum methods.
- Compare the classical CFD approach with quantum methods, including QLB, HHL, and related techniques.
- Identify the technologies, algorithms, and resources required for this approach.
- Evaluate the advantages and limitations of quantum-based CFD solutions.
- Compare the quantum approach with LBM as well.

### Shape optimization
- Optimization of wing, rocket, diffuser, and similar geometries using machine learning.
- Comparison of quantum and classical models.

### Measure accuracy
- Classical method
- Classical samples
- Classical median error (bps)
- QPU median point error (bps)
- QPU median point + CL (bps)
- QPU certified success

### Attributes
- 2D space
- 512 × 512 px area size
- Wing and upper section of the rocket
- Wing profile
- LBM
- Ansys/OpenFOAM
- QLBM
- NACA 0012 profile
- Length: 128 cells
- Angle: 0°, 5°, 10°
- Left to right
- Laminar, incompressible
- 512 × 512
- Determine how to add the wing to the quantum state

![alt text](<Images/Screenshot 2026-07-24 at 13.16.16.png>)



### Area boundary attributes
- 512 × 512 cells for LBM and QLBM; Ansys/OpenFOAM: 256 × 256 meters
- 10 shots = 1 second, 1 shot = 0.1 second, velocity = 1 m/s in ANSYS/OpenFOAM, LBM, and QLBM
- Default temperature: 23°C

### Constants
- Air density: 1.225 kg/m³
- Dynamic viscosity of air: 1.81 × 10⁻⁵ Pa·s
- Kinematic viscosity of air: 1.48 × 10⁻⁵ m²/s
- Speed of sound in air: 343 m/s
- Specific heat ratio: 1.4
- Gas constant for air: 287 J/(kg·K)
- Reference temperature: 23°C