using ModelingToolkit, OrdinaryDiffEq, Symbolics
using Symbolics: solve_for, solve_single_eq
using Plots
# using DifferentialEquations
using Unitful
@info "usings"
Lunit = u"m"
# 1d elastic collision with 2 masses with radius
const to = TimerOutput()
# for My in [1, 10, 100]
My = 1
Mx = 1
@parameters t rx = 0.5 ry = 0.5 mx = Mx my = My wall1 = 0 wall2 = 10
sts = @variables x(t) = 2.0 y(t) = 5.0 vx(t) = 0.0 vy(t) = -1.0 #n_collisions(t) = 0
D = Differential(t)
eqs = [
    D(y) ~ vy,
    D(vy) ~ 0,
    D(x) ~ vx,
    D(vx) ~ 0,
    # D(n_collisions) ~ 0
]

# global N_COLLISIONS = 0

function affect!(integ, u, p, ctx)
    @variables vxx vyy
    vy = integ.u[u.vy]
    vx = integ.u[u.vx]
    mx, my = integ.p[[p.mx, p.my]]
    eq1 = mx * vx + my * vy ~ mx * vxx + my * vyy
    eq2 = 1 // 2 * mx * (vx^2) + 1 // 2 * my * (vy^2) ~ 1 // 2 * mx * (vxx^2) + 1 // 2 * my * (vyy^2)

    eq2 = substitute(eq2, Dict(vxx => solve_for(eq1, vxx)))
    vyysol = simplify.(solve_single_eq(eq2, vyy.val))
    @info u p vyysol
    filt = filter(x -> x.rhs != vy, vyysol)
    new_vy = filt[1].rhs
    new_vx = solve_for(substitute(eq1, Dict(vyy => new_vy)), vxx)

    integ.u[u.vx] = Symbolics.value(new_vx)
    integ.u[u.vy] = Symbolics.value(new_vy)
    # integ.u[u.n_collisions] += 1
    # N_COLLISIONS += 1
    nothing
end

continuous_events = [
    [(x - rx) ~ 0] => [vx ~ -vx],
    # , n_collisions ~ n_collisions + 1
    # ],
    # [(x + rx) ~ (y - ry)] => (affect!, [vx, vy, n_collisions], [mx, my], nothing)
    # [(y - ry) ~ (x + rx)] => (affect!, [vx, vy, n_collisions], [mx, my], nothing)
    [(x + rx) ~ (y - ry)] => (affect!, [vx, vy], [mx, my], nothing),
    [(y - ry) ~ (x + rx)] => (affect!, [vx, vy], [mx, my], nothing)
    # [(y + ry) ~ 10] => [vy ~ -vy, n_collisions ~ n_collisions + 1]
]
@named elastic = ODESystem(eqs, t, sts, [rx, ry, mx, my]; continuous_events)
@info "sys $My"
old_sts = states(elastic)
sys = structural_simplify(elastic) # simplify deletes `n_collisions` state
tspan = (0, 50)
prob = ODEProblem(sys, [], tspan)
solve(prob,Rosenbrock23())

try 
    solve(prob,Rosenbrock23())
catch 
end
    # @timeit to "$My" solve(prob, Rosenbrock23())
# end

# why don't these work 
prob = ODEProblem(elastic, [], tspan, [my => 1]; saveat=0.1)
prob = ODEProblem(elastic, [], tspan, [My => 1]; saveat=0.1)
solve(prob, Rosenbrock23()) # MethodError: no method matching resize!(::Vector{Equation}, ::SymbolicUtils.Add{Int64, Int64, Dict{Any, Number}, Nothing})

probs = [ODEProblem(elastic, [], tspan, [my => n]; saveat=0.1) for n in [1, 10, 100]]
@time solve(probs[1], Rosenbrock23())
# f = m->solve(, Rosenbrock23())
function time_scaling(f, xs) #; save=false)
    to = TimerOutput()
    # ys = []
    for (i, x) in enumerate(xs)
        y = @timeit to "$i" f(x)
    end
    timerout_to_df(to)
end
@time f(1)
# time_scaling(f, [1, 10, 100])

prob = ODEProblem(elastic, [], tspan, [My=>100]; saveat=0.1)
prob = ODEProblem(elastic, [], tspan, [My=>100]; saveat=0.1)
prob = ODEProblem(elastic, [], tspan, [My=>100]; saveat=0.1)
prob2 = ODEProblem(elastic, [], tspan, [My=>1000]; saveat=0.1)
# prob = ODEProblem(elastic, [], tspan; saveat=0.1)
@info "prob"
# prob = ODEProblem(sys, [], tspan; saveat=0.1)
sol = solve(prob, Rosenbrock23())
@info "sol"
# plot(sol)
nsteps = 500
saves = range(tspan..., length=nsteps)
solitp = sol(saves)

@test_throws Any solm = Matrix(sol)
solm = Array(sol)
@test solm isa Matrix

# anim_mat = reduce(hcat, sol(saves))'

anim = @animate for i in 1:length(solitp)
    cur_vx = solitp[vx][i]
    cur_vy = solitp[vy][i]
    last_ten = maximum([1, i - 10]):i
    nlast = length(last_ten)
    ncs = round(Int, solitp[n_collisions][i])
    curt = solitp.t[i]
    plot_t = round(Int, curt)
    myplot = plot(solitp[x][last_ten], zeros(nlast); title="collisions: $(ncs) | t = $(plot_t)", xlims=(0, 10), ylims=(-5, 5))#, $([solitp.u[1, i], solitp.u[2, i]])")
    scatter!(myplot, [solitp[x][i]], [0.0]; ms=log(Mx) + 1)#, texts=["$cur_vx"])
    plot!(myplot, solitp[y][last_ten], zeros(nlast))
    scatter!(myplot, [solitp[y][i]], [0.0]; ms=log(My) + 1)#, texts=["$cur_vy"])
end
@info "anim"

mp4(anim, "elastic_with_radius_$(My).mp4", fps=60)

## indexing interface 
# :ODECompositeSolution
# indexing with time(s) or 
@test SciMLBase.AbstractTimeseriesSolution <: AbstractDiffEqArray

sol1 = sol[1] # first timestep
solitp = sol(0) # values of u for t=0
@test sol1 == solitp
@which sol[1]
@which sol(0)

# we should allow
sol[1][x] # get symbolic x at first timestep
sol(0)[x] # same, but with interpolation
sol[idxs, :][x] # get symbolic x at all timesteps in idxs
sol[[x, vx], :]
sol[[x, vx], [1, 2, 3]]
sol([x, vx])(0:0.1:tspan[2])

# 1) i think this just means removing the special sol(0) and make it sol([0])
# i dont think that would be breaking 

# 2) make getindex methods return DiffEqArrays 

ts = [1, 2, 3]
@which sol(ts)

solitp2 = sol(ts)
solitp2[x] # interesting that when passing ts as vector, it stays a diffeqarray, allowing symbolic indexing
@test_throws Any sol(ts)[[x, vx]]

idxs = [2, 3, 4]

# this is why DataFrames disallows df[1], df[[1,2,3]]
# if sol[j] is allowed, its strange that sol[[j]] is not allowed
sol[idxs] # solution_slice methoderror. definitely need better error message here
# what you want here is 
sol[idxs, :]


ss = sol[idxs, :] # remember time is last component, unlike DataFrame(sol)
@test_throws ArgumentError ss[x]
sol[:, [1, 2, 3]] # get all states, for timesteps 1, 2, 3

@test_throws Any sol[:, [1.5, 2, 3]] # doesn't implicitly interpolate when using getindex

sol(5.5; idxs=idxs)

arr = sol(saves)
arr[vx]
sol[[x, vx], 1:10]


@test_throws Any sol[[x, vx], 1:10]

sol(ts)[x]
solt

xvx = sol[[x, vx]]

# so you can interpolate only a subset of states (with indices or symbolic vars)
sol(1, idxs=[x, vx])
sol([1, 2], idxs=[x, vx])

# is it faster
@btime sol(0:0.1:tspan[2], idxs=[x])
@btime sol(0:0.1:tspan[2])

# we should allow people to go the other way
sol([x, y])(0:0.1:tspan[2])



# this prints out a hilarious amount of stuff
sol.interp
