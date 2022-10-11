using ModelingToolkit, OrdinaryDiffEq, Symbolics
using Symbolics: solve_for, solve_single_eq
using Plots
using DifferentialEquations

# ball bouncing in circle
Mx = 1
My = 100
R1 = 1 # start with no width
R2 = 5
G = 10

# My = 101 # with 101, the simulation never finishes
@parameters t rx = 0.5 ry = 0.5 mx = Mx my = My r = R1 g = G
sts = @variables x(t) = 1 y(t) = -1 vx(t) = 0.0 vy(t) = 0
D = Differential(t)
eqs = [
    D(y) ~ vy,
    D(x) ~ vx,
    D(vy) ~ -g,
    D(vx) ~ 0
]

function affect_bounce!(integ, u, p, ctx)
    x_t, vx_t, y_t, vy_t = integ.u[[u.x, u.vx, u.y, u.vy]]
    theta = atan(y_t, x_t)
    rm = [cos(theta) sin(theta); -sin(theta) cos(theta)]
    rmi = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    vxn, vyn = rm * [vx_t; vy_t]
    vyn -= vyn
    vxn, vyn = rmi * [vxn; vyn]
    integ.u[u.vx] = vxn
    integ.u[u.vy] = vyn
    # @info "triggered" x_t, vx_t, y_t, vy_t, vxn, vyn

    global N_COLLISIONS += 1
    nothing
end

continuous_events = [
    [(x^2 + (y)^2)^1 / 2 ~ r] => (affect_bounce!, [vx, vy, x, y], [], nothing)
    # [y ~ -5] => [vy ~ -vy] # this affect works 
]

global N_COLLISIONS = 0
@named ball = ODESystem(eqs, t, sts, [mx, g, r]; continuous_events)
# @named ball = ODESystem(eqs; continuous_events)
old_sts = states(ball)
sys = structural_simplify(ball)
tspan = (0, 5)
prob = ODEProblem(sys, [], tspan; saveat=0.1)
sol = solve(prob)
N_COLLISIONS
# plot(sol)
nsteps = 100
saves = range(tspan..., length=nsteps)
solitp = sol(saves)

anim = @animate for i in 1:length(solitp)

    cur_vx = solitp[vx][i]
    cur_vy = solitp[vy][i]
    cur_x = solitp[x][i]
    cur_y = solitp[y][i]

    last_ten = maximum([1, i - 10]):i
    nlast = length(last_ten)
    # ncs = round(Int, solitp[n_collisions][i])
    curt = solitp.t[i]
    plot_t = round(Int, curt)
    # myplot = plot(solitp[x][last_ten], , zeros(nlast); title="t = $(plot_t)", xlims=(0, 10), ylims=(-5, 5))#, $([solitp.u[1, i], solitp.u[2, i]])")
    hyp = sqrt(cur_x ^2 + cur_y^2)
    myplot = scatter([cur_x], [cur_y], title="t = $(plot_t), $hyp"; xlims=(-6, 6), ylims=(-6, 6))
    plot!(myplot, 5 * cos.(0:0.01:2pi), 5 * sin.(0:0.01:2pi); xlims=(-6, 6), ylims=(-6, 6))
    # plot(5 * cos.(0:0.01:2pi), 5* sin.(0:0.01:2pi); xlims=(-6, 6), ylims=(-6, 6))
    # scatter!(myplot, [solitp[x][i]], [0.0]; ms=log(Mx) + 1)#, texts=["$cur_vx"])
    # plot!(myplot, solitp[y][last_ten], zeros(nlast))
    # scatter!(myplot, [solitp[y][i]], [0.0]; ms=log(My) + 1)#, texts=["$cur_vy"])
end
mp4(anim, "ball_in_circle.mp4", fps=30)
@info "anim"
