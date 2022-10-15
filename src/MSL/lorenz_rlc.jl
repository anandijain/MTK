using DifferentialEquations, ModelingToolkit, Symbolics, ModelingToolkitStandardLibrary, LinearAlgebra
using ModelingToolkitStandardLibrary: Blocks, Blocks.Constant
using ModelingToolkitStandardLibrary.Electrical

using Plots

const MTK = ModelingToolkit
const MSL = ModelingToolkitStandardLibrary



@parameters t
D = Differential(t)
@variables Lt(t) Ct(t) = 100 Vt(t) = 100

function TimeDependentInductor(; name, i_start=0.0)
    @named oneport = OnePort(; i_start=i_start)
    @unpack v, i = oneport
    sts = @variables Lt(t)
    eqs = [
        D(i) ~ 1 / Lt * v,
    ]
    extend(ODESystem(eqs, t, sts, []; name=name), oneport)
end

function TimeDependentCapacitor(; name, v_start=0.0)
    @named oneport = OnePort(; v_start=v_start)
    @unpack v, i = oneport
    # pars = @parameters C = C
    sts = @variables Ct(t)
    eqs = [
        D(v) ~ i / Ct,
    ]
    extend(ODESystem(eqs, t, sts, []; name=name), oneport)
end

function TimeDependentResistor(; name)
    @named oneport = OnePort()
    @unpack v, i = oneport
    # pars = @parameters R = R
    sts = @variables Rt(t)
    eqs = [
        v ~ i * Rt,
    ]
    extend(ODESystem(eqs, t, sts, []; name=name), oneport)
end



@named o = Constant(; k=1.0)
# @named o = Sine(; frequency=1)
@named r = TimeDependentResistor()
@named v = Voltage()
@named c = TimeDependentCapacitor()
@named g = Ground()
@named l = TimeDependentInductor()

using ModelingToolkit, DifferentialEquations, Plots
@parameters t sig = 10 rho = 28.0 beta = 8 / 3
@variables x(t) = 100 y(t) = 1.0 z(t) = 1
D = Differential(t)
leqs = [D(x) ~ sig * (y - x),
    D(y) ~ x * (rho - z) - y,
    D(z) ~ x * y - beta * z]

@named lorenz = ODESystem(leqs)

# prob = ODAEProblem(lorenz, [0], (0, 10.0))
# prob = ODEProblem(lorenz, [], (0, 100.0))
# sol = solve(prob) # err
# plot(sol)
# prob = ODEProblem(lorenz, [0], (0, 10.0)) ArgumentError: Equations (3), states (3), and initial conditions (1) are of different lengths
# sol = solve(prob)

Lcircut = [
    connect(o.output, v.V)
    connect(v.p, l.p)
    connect(l.n, v.n, g.g)]

@named lcirc = ODESystem(Lcircut, t; systems=[o, v, l, g])

lleqs = [
    lcirc.l.Lt ~ (abs(lorenz.x) + 1) / 1000
]

@named llsys = ODESystem(lleqs, t)
@named con = compose(llsys, [lcirc, lorenz])

scon = ModelingToolkit.structural_simplify(con)
equations(scon)
# prob = ODAEProblem(scon, [0], (0, 10.0))
prob = ODAEProblem(scon, [], (0, 1000.0))
sol = solve(prob)
plot(sol)
# plot(sol, vars=[x, y, z])

rlc_eqs = [
    connect(o.output, v.V)
    connect(v.p, r.p)
    connect(r.n, l.p)
    connect(l.n, c.p)
    connect(c.n, v.n, g.g)]

@named rlc = ODESystem(rlc_eqs, t; systems=[o, v, r, l, c, g])

lleqs = [
    rlc.l.Lt ~ (abs(lorenz.x) + 1) / 1000,
    # D(rlc.r.Rt) ~ 0,
    rlc.c.Ct ~ (abs(lorenz.z) + 1) / 1000,
    # rlc.v.V.u ~ (abs(lorenz.x) + 1) / 1000,
    rlc.r.Rt ~ (abs(lorenz.x) + 1) / 1000,
]


@named llsys = ODESystem(lleqs, t)
@named con = compose(llsys, [rlc, lorenz])

scon = ModelingToolkit.structural_simplify(con)
equations(scon)
prob = ODAEProblem(scon, [rlc.r.Rt => 100], (0, 10.0))
probode = ODEProblem(scon, [rlc.r.Rt => 100], (0, 10.))
sol = solve(prob)
sol = solve(probode)
plot(sol, vars=[rlc.v.i, rlc.c.v])

sys = scon
# function add_zero_diffeqs(sys)
sts = states(sys)
varmap = ModelingToolkit.defaults(sys)
missingvars = setdiff(varlist, keys(varmap))

# end