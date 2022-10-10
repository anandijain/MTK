using ModelingToolkit

@parameters x
@variables t u(..)
Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dt = Differential(t)

#2D PDE
C=1
eq  = Dtt(u(t,x)) ~ C^2*Dxx(u(t,x))

# Initial and boundary conditions
bcs = [u(t,0) ~ 0.,# for all t > 0
       u(t,1) ~ 0.,# for all t > 0
       u(0,x) ~ x*(1. - x), #for all 0 < x < 1
       Dt(u(0,x)) ~ 0. ] #for all  0 < x < 1]

# Space and time domains
domains = [t ∈ (0.0,1.0),
           x ∈ (0.0,1.0)]

@named pde_system = PDESystem(eq,bcs,domains,[t,x],[u])


# heat 
using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, DomainSets

# Parameters, variables, and derivatives
@parameters t x y
@variables u(..)
Dt = Differential(t)
Dxx = Differential(x)^2
# Dyy = Differential(y)^2

# 1D PDE and boundary conditions
eq  = Dt(u(t,x)) ~ Dxx(u(t,x))
bcs = [u(0,x) ~ cos(x),
       u(t,0) ~ exp(-t),
       u(t,Float64(pi)) ~ -exp(-t)]

# Space and time domains
domains = [t ∈ Interval(0.0,1.0),
           x ∈ Interval(0.0,Float64(pi))]

# PDE system
@named pdesys = PDESystem(eq,bcs,domains,[t,x],[u(t,x)])

# Method of lines discretization
dx = 0.1
order = 2
discretization = MOLFiniteDifference([x=>dx],t;centered_order=order)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
sol = solve(prob,Tsit5(),saveat=0.1)


x = (0:0.1:1)[2:end-1]
t = sol.t



using Plots
anim = @animate for i in 1:length(t)
       p1 = plot(x,sol.u[i][1:9],label="u, t=$(t[i])";legend=false,xlabel="x",ylabel="u",ylim=[0,1])
       p2 = plot(x,sol.u[i][10:end],label="v, t=$(t[i])";legend=false,xlabel="x",ylabel="v",ylim=[0,1])
       plot(p1,p2)
end
gif(anim, "plot.gif",fps=30)