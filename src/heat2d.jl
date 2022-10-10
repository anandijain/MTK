using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets
# Method of Manufactured Solutions: exact solution
# u_exact = (x,t) -> exp.(-t) * cos.(x)

# Parameters, variables, and derivatives
N = 2
@parameters t x y v[1:N]
@variables u(..)

Dt = Differential(t)
Ds = Differential.(v)
ds = collect(Ds)
Dx = Differential(x)
Dy = Differential(y)
Δ =  sum(map(x->(x^2)(u(t, v...)), ds))
# function Δ(u, t, v)

eq = Dt(u(t, v...)) ~ Δ

# Dxx = Dx^2
# Dyy = Dy^2

# ∇ = (Dx, Dy)
# ∇² = (Dxx, Dyy)
# d2 = ∇²
# d2 .* u(t, x, y)

# eq  = Dt(u(t, x, y)) ~ Dxx(u(t, x, y) + 

bcs = [u(0, v...) ~ cos(x),

        u(t, 0, v[2]) ~ exp(-t),
        u(t, v[1], 0) ~ exp(-t),

        u(t, 1, 1) ~ exp(-t) * cos(1)
        ]

# Space and time domains
domains = [t ∈ Interval(0.0, 1.0),
           v[1] ∈ Interval(0.0, 1.0),
           v[2] ∈ Interval(0.0, 1.0)
           ]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)])

# Method of lines discretization
dx = dy = 0.1
order = 2
discretization = MOLFiniteDifference([v[1] => dx, v[2]=> dy], t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob, Tsit5(), saveat=0.2)

# Plot results and compare with exact solution
discrete_x = sol[x]
discrete_t = sol[t]
solu = sol[u(t, x)]