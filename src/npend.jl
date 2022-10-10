using ModelingToolkit, LinearAlgebra, DifferentialEquations, Plots, Test 

function Pendulum(;name, g=9.81, r=1.)
    gval, rval = g, r
    @parameters g r
    @variables x(t) y(t) T(t)
    D2 = Differential(t)^2
    eqs = [ D2(x) ~ T * x
            D2(y) ~ T * y - g ]
    ODESystem(eqs, t, [x, y, T], [g, r]; name=name, defaults=(g => gval, r => rval))
end

function connect_pendulums(pends)
    eqs = [0 ~ pends[1].x^2 + pends[1].y^2 - pends[1].r^2]
    for i in 2:length(pends)
        push!(eqs, connect_eq(pends[i - 1:i]...))
    end
    eqs
end

connect_eq(a, b) = 0 ~ (a.x - b.x)^2 + (a.y - b.y)^2 - b.r^2 

@parameters t
D = Differential(t)

N = 3 # change this to add or remove pendulumns
pends = []
for i in 1:N
    push!(pends, Pendulum(;name=Symbol("pendulum", i)))
end

npendeqs = connect_pendulums(pends)

@named npend = compose(ODESystem(npendeqs, t; name=:connections), pends...)
@named npend2 = ODESystem(npendeqs, t; systems=pends)

u0 = Pair[]
for (i, p) in enumerate(pends)
    pairs = [
        p.x => i
        p.y => 0
        D(p.x) => 0
        D(p.y) => 0
        p.T => 0
    ]
    append!(u0, pairs)
end

pendulum_sys = structural_simplify(ode_order_lowering(flatten(npend)))
prob = ODAEProblem(pendulum_sys, u0, (0, 100.0))
sol = solve(prob, Tsit5())
plot(sol)
sol = solve(prob)