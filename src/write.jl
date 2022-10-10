using ModelingToolkit, OrdinaryDiffEq

@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(D(x)) ~ σ*(y-x),
        D(y) ~ x*(ρ-z)-y,
        D(z) ~ x*y - β*z]

@named sys = ODESystem(eqs)
tspan=(0,10.)

sys = ode_order_lowering(sys)

open("sys.jl", "w") do io
    write(io, sys)
end

io = IOBuffer()
s = String(take!(io))
