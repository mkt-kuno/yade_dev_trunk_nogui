"""
Check test for ThermalEngine with periodic conditions, considering the following:
- Particles with different initial temperatures.
- Thermal conduction between particles.
- Thermal expansion of particles.
"""
from yade import pack, ymport, plot, utils, export, timing
import numpy as np

if ('THERMAL' in features):
# PARAMETERS ==============================================================================
    # Material properties
    mat_young           = 1e6
    mat_poisson         = 0.3
    mat_friction        = 0.5
    mat_density         = 1.0
    mat_heat_capacity   = 10.0
    mat_conductivity    = 100.0
    mat_expansion_coeff = 1e-4

    # Sizes
    radius_avg = 0.08
    radius_var = 0.30
    cell_size  = 1.0

    # Analysis settings
    time_step_ratio = 0.1
    limit_steps     = 10000
    temp_ref        = 364.1620967973047
    tolerance       = 1e-2

# EXECUTION ===============================================================================
    # Create materials
    mat_particles = FrictMat(young=mat_young, poisson=mat_poisson, density=mat_density, frictionAngle=mat_friction)
    O.materials.append(mat_particles)

    # Create particles
    O.periodic = True
    sp = pack.randomPeriPack(initSize=Vector3(cell_size,cell_size,cell_size), radius=radius_avg, rRelFuzz=radius_var, seed=1)
    sp.toSimulation(material=mat_particles)

    # Set analysis parameters
    timing.reset()
    O.dt = time_step_ratio * PWaveTimeStep()
    O.dynDt = False
    O.trackEnergy = True

    # Set engines
    O.engines = [
        ForceResetter(),
        InsertionSortCollider(
            [Bo1_Sphere_Aabb()]
        ),
        InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom()],
            [Ip2_FrictMat_FrictMat_FrictPhys()],
            [Law2_ScGeom_FrictPhys_CundallStrack()]
        ),
        FlowEngine(label='engineFlow', dead=True, updateTriangulation=False),
        ThermalEngine(label = 'engineThermal', dead = False,
            particleT0          = float('nan'), # If NaN, set individually for each particle
            particleCp          = float('nan'), # If NaN, set individually for each particle
            particleK           = float('nan'), # If NaN, set individually for each particle
            particleAlpha       = float('nan'), # If NaN, set individually for each particle
            tsSafetyFactor      = 0.0,
            conduction          = True,
            useBoBMethod        = True,
            thermoMech          = True,
            advection           = False,
            fluidBeta           = 0.0,
            boundarySet         = True,
            flowTempBoundarySet = True,
        ),
        NewtonIntegrator(label='engineNI', dead=False,
            gravity = (0,0,0),
            damping = 0.0
        )
    ]

    # Set particle thermal properties individually
    for b in O.bodies:
        if isinstance(b.shape, Sphere):
            b.state.temp  = 400.0 if 0.0 < b.state.pos[0] < 0.5 * cell_size else 300.0
            b.state.Cp    = mat_heat_capacity
            b.state.k     = mat_conductivity
            b.state.alpha = mat_expansion_coeff

    # Run simulation
    O.run(limit_steps, wait=True)

    # Get average temperature
    temp_avg = np.average([b.state.temp for b in O.bodies if isinstance(b.shape, Sphere)])
    print('Average temperature = ', temp_avg)

    # Compare current and reference results
    if abs(temp_avg-temp_ref)/temp_ref > tolerance:
        raise YadeCheckError('ThermalEngine checktest: average temperature incorrect')

else:
    print("This checkThermalEngine.py cannot be executed because ENABLE_THERMAL is disabled")
