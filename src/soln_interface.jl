
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

Matrix(sol)
Array(sol) isa Matrix
# sol[[1,2], :]
sol[[x, vx], [1, 2, 3]]

sol[[x, vx]](0:0.1:tspan[2])


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


using Test, ModelingToolkit, OrdinaryDiffEq

@parameters t σ=1 ρ=2 β=3
@variables x(t)=1 y(t)=1 z(t)=1
D = Differential(t)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]
@named sys = ODESystem(eqs)
prob = ODEProblem(sys, [], (0., 10.))
sol = solve(prob, Tsit5())
@test_throws Any Matrix(sol)
@test Array(sol) isa Matrix