include("../src/systems.jl")
include("../src/visuals.jl")

# Create 3 pendulum objects with some different properties
p1 = Pendulum(1; L=[1.5], θ=[π/4], ω=[0.0]) # A single pendulum
p2 = Pendulum(2; L=[1.0, 1.0], θ=[π/6, π/3], ω=[0.0, 0.0]) # A double pendulum
p3 = Pendulum(3; L=[0.8, 0.8, 0.8], θ=[π/8, π/4, π/2], ω=[0.0, 0.0, 0.0]) # A triple pendulum

# Create a `Simulation` object containing said pendulums
sim = Simulation([p1, p2, p3])

# Simulate those pendulums!
animated_multi_pendulum(sim)


# Creating some similar double/triple with slightly perturbed initial conditions can be quite interesting!
p1 = Pendulum(2; L=[1.0, 1.0], θ=[π/6, π/3], ω=[0.0, 0.0])
p2 = Pendulum(2; L=[1.0, 1.0], θ=[π/6, π/3] .+ 0.01, ω=[0.0, 0.0])
p3 = Pendulum(2; L=[1.0, 1.0], θ=[π/6, π/3] .+ 0.01, ω=[0.0, 0.0])
sim = Simulation([p1, p2, p3])
animated_multi_pendulum(sim)