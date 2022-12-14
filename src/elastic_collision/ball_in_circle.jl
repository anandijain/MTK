using ModelingToolkit, OrdinaryDiffEq, Symbolics
using Symbolics: solve_for, solve_single_eq
using Plots
using DifferentialEquations

# ball bouncing in circle
Mx = 1
R1 = 5
G = 10
G = 10
x0 = R1 / sqrt(2)
# My = 101 # with 101, the simulation never finishes
@parameters t rx = 0.5 ry = 0.5 mx = Mx r = R1 g = G
sts = @variables x(t) = x0 y(t) = 0 vx(t) = 0.0 vy(t) = 0
D = Differential(t)
eqs = [
    D(y) ~ vy,
    D(x) ~ vx,
    D(vy) ~ -g/mx,
    # D(vy) ~ 0,
    D(vx) ~ 0
]

function condition(u, t, integrator)
    x = u[3]
    y = u[1]
    r = integrator.p[3]
    x^2 + y^2 - r^2
end

function affect_bounce!(integ, u, p, ctx)
    x_t, vx_t, y_t, vy_t = integ.u[[u.x, u.vx, u.y, u.vy]]
    theta = atan(y_t, x_t)
    rm = [cos(theta) sin(theta); -sin(theta) cos(theta)]
    rmi = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    vxn, vyn = rm * [vx_t, vy_t]
    # vyn -= vyn
    vyn *= -1
    # vxn -= vxn
    vxn, vyn = rmi * [vxn, vyn]
    integ.u[u.vx] = -vxn
    integ.u[u.vy] = -vyn
    # integ.u[u.vx] = vx_t*cos(theta) + vy_t*sin(theta)
    # integ.u[u.vy] = -vy_t*sin(theta) + vx_t*cos(theta)
    @info "triggered" theta (x_t, y_t), (vx_t, vy_t), vxn, vyn

    global N_COLLISIONS += 1
    nothing
end

continuous_events = [
    [(x^2 + y^2) ~ r^2] => (affect_bounce!, [vx, vy, x, y], [], nothing)
    # [y ~ -5] => [vy ~ -vy] # this affect works 
]
# cbs = CallbackSet(ContinuousCallback(condition, affect_bounce!))
global N_COLLISIONS = 0
@named ball = ODESystem(eqs, t, sts, [mx, g, r]; continuous_events)
# @named ball = ODESystem(eqs, t, sts, [mx, g, r])# ; continuous_events)
old_sts = states(ball)
sys = structural_simplify(ball)
tspan = (0, 10)
prob = ODEProblem(sys, [], tspan; saveat=0.01)#, callback=cbs)
sol = solve(prob)
N_COLLISIONS
# plot(sol)
nsteps = 1000
saves = range(tspan..., length=nsteps)
solitp = sol(saves)

anim = @animate for i in 1:length(solitp)

    cur_vx = solitp[vx][i]
    cur_vy = solitp[vy][i]
    cur_x = solitp[x][i]
    cur_y = solitp[y][i]

    last_ten = maximum([1, i - 10]):i
    nlast = length(last_ten)
    curt = solitp.t[i]
    plot_t = round(Int, curt)
    hyp = sqrt(cur_x ^2 + cur_y^2)
    myplot = scatter([cur_x], [cur_y], title="t = $(plot_t), $hyp"; xlims=(-6, 6), ylims=(-6, 6), aspect_ratio=:equal)
    plot!(myplot, solitp[x][last_ten], solitp[y][last_ten], title="t = $(plot_t), $hyp"; xlims=(-6, 6), ylims=(-6, 6))
    plot!(myplot, 5 * cos.(0:0.01:2pi), 5 * sin.(0:0.01:2pi); xlims=(-6, 6), ylims=(-6, 6))
end
mp4(anim, "ball_in_circle.mp4", fps=60)
@info "anim"
