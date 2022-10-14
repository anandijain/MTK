# https://archimede.uniba.it/~testset/problems/transamp.php
using OrdinaryDiffEq, ModelingToolkit, Symbolics
using LinearAlgebra
using ModelingToolkitStandardLibrary
using Plots

const MTK = ModelingToolkit
const MSL = ModelingToolkitStandardLibrary

using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Constant

using DifferentialEquations, ModelingToolkit

Uₑ(t) = 0.1sin(200π * t)
g(x) = β * (exp(x / UF) - 1)
@parameters t
D = Differential(t)
ModelingToolkit.@parameters begin
    Ub = 6.0
    UF = 0.026
    α = 0.99
    β = 1e-6
    R0 = 1e3
    Rk[1:9] = [9e3 for _ in 1:9]
    Ck[1:5] = [i * 1e-6 for i in 1:5]
end
Capacitor

C = []
for i in 1:5
    c = Capacitor(; name=Symbol(:C, i), C=Ck[i])
    push!(C, c)
end

rc_eqs = [
    connect(constant.output, source.V)
    connect(source.p, resistor.p)
    connect(resistor.n, capacitor.p)
    connect(capacitor.n, source.n, ground.g)
]

@named o = Constant(; k=1.0)
@named o = Sine(; frequency=Inf)
@named r = Resistor(; R=100)
@named v = Voltage()
@named c = Capacitor(; C=100)
@named g = Ground()
@named l = Inductor(L=100)

rlc_eqs = [
    connect(o.output, v.V)
    connect(v.p, r.p)
    connect(r.n, l.p)
    connect(l.n, c.p)
    connect(c.n, v.n, g.g)
]

eqs = [
    connect(o.output, v.V)
    # connect(v.p, r.p, c.p)

    connect(v.p, r.p)
    connect(v.p, c.p)


    # connect(c.n, r.n, v.n, g.g)

    connect(c.n, r.n)
    connect(r.n, v.n)
    connect(v.n, g.g)
]


@named sys = ODESystem(eqs, t; systems=[r, c, o, v, g])

@named sys = ODESystem(rlc_eqs, t; systems=[o, v, r, l, c, g])
ssys = structural_simplify(sys)
equations(ssys)
prob = ODAEProblem(ssys, [], (0, 10.0))
sol = solve(prob, Tsit5(); saveat=0:0.01:10)
plot(sol[o.output.u])

