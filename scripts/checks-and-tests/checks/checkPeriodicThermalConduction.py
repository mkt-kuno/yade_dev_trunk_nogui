"""
Check test for ThermalEngine with periodic conditions,
The test involves the compression of a periodic assembly of particles considering the following:
- Particles with different initial temperatures.
- Thermal conduction between particles.
- Thermal expansion of particles.
- Heat generation by energy dissipation during particle interactions.
"""
from yade import pack, ymport, plot, utils, export, timing
import numpy as np

if ('THERMAL' in features):
# PARAMETERS ==============================================================================
    # Material properties
    mat_young           = 1e6
    mat_poisson         = 0.3
    mat_friction        = 0.1
    mat_density         = 1.0
    mat_heat_capacity   = 10.0
    mat_conductivity    = 1000.0
    mat_expansion_coeff = 1e-4

    # Sizes
    radius_avg = 0.08
    radius_var = 0.30
    cell_size  = 1.0

    # Analysis settings
    time_step_ratio       = 0.1
    thermal_frequency     = 2
    temperature_avg       = 300.0
    temperature_dev       = 10.0
    heat_generation_ratio = 10.0
    strain_rate           = -0.1
    limit_steps           = 10000
    temp_ref              = 312.2663320345086
    tolerance             = 1e-2

# EXECUTION ===============================================================================
    # Create materials
    mat_particles = FrictMat(young=mat_young, poisson=mat_poisson, density=mat_density, frictionAngle=mat_friction)
    O.materials.append(mat_particles)

    # Set boundary properties
    O.periodic = True
    O.cell.velGrad = Matrix3(strain_rate,strain_rate,strain_rate,strain_rate,strain_rate,strain_rate,strain_rate,strain_rate,strain_rate)
    
    # Create particles
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
            thermalFreq         = thermal_frequency,
            heatGenerationRatio = heat_generation_ratio,
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
    spheres = [b for b in O.bodies if isinstance(b.shape, Sphere)]
    np.random.seed(1)
    temps = np.random.normal(loc=temperature_avg, scale=temperature_dev, size=len(spheres))
    for b, temp in zip(spheres, temps):
        b.state.temp  = temp
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
