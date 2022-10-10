using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, DomainSets
@parameters t x
@parameters Dn, Dp
@variables u(..) v(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

eqs  = [Dt(u(t,x)) ~ Dn * Dxx(u(t,x)) + u(t,x)*v(t,x), 
        Dt(v(t,x)) ~ Dp * Dxx(v(t,x)) - u(t,x)*v(t,x)]
bcs = [u(0,x) ~ sin(pi*x/2),
       v(0,x) ~ sin(pi*x/2),
       u(t,0) ~ 0.0, Dx(u(t,1)) ~ 0.0,
       v(t,0) ~ 0.0, Dx(v(t,1)) ~ 0.0]

domains = [t ∈ Interval(0.0,1.0),
           x ∈ Interval(0.0,1.0)]

@named pdesys = PDESystem(eqs,bcs,domains,[t,x],[u(t,x),v(t,x)],[Dn=>0.5, Dp=>2])
discretization = MOLFiniteDifference([x=>0.1],t)
prob = discretize(pdesys,discretization)
sol = solve(prob,Tsit5())
