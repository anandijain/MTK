using ModelingToolkit, OrdinaryDiffEq, Symbolics
using Symbolics: solve_for, solve_single_eq
using Plots
using DifferentialEquations

# 1d elastic collision with 2 masses
Mx = 1
My = 1000
@parameters t rx = 0.5 ry = 0.5 mx = Mx my = My wall1 = 0 wall2 = 10
sts = @variables x(t) = 2.0 y(t) = 5.0 vx(t) = 0.0 vy(t) = -1.0
D = Differential(t)
eqs = [
    D(y) ~ vy,
    D(vy) ~ 0,
    D(x) ~ vx,
    D(vx) ~ 0
]

function affect!(integ, u, p, ctx)
    @variables vxx vyy
    vy = integ.u[u.vy]
    vx = integ.u[u.vx]
    mx, my = integ.p[[p.mx, p.my]]
    eq1 = mx * vx + my * vy ~ mx * vxx + my * vyy
    eq2 = 1 // 2 * mx * (vx^2) + 1 // 2 * my * (vy^2) ~ 1 // 2 * mx * (vxx^2) + 1 // 2 * my * (vyy^2)

    eq2 = substitute(eq2, Dict(vxx => solve_for(eq1, vxx)))
    vyysol = simplify.(solve_single_eq(eq2, vyy.val))
    new_vy = filter(x -> x.rhs != vy, vyysol)[1].rhs
    new_vx = solve_for(substitute(eq1, Dict(vyy => new_vy)), vxx)

    integ.u[u.vx] = Symbolics.value(new_vx)
    integ.u[u.vy] = Symbolics.value(new_vy)
    global N_COLLISIONS += 1
    @info (integ.u[u.x], vx), (integ.u[u.y], vy), N_COLLISIONS

    term_affect!(integ, u, p, ctx)
    nothing
end

# condition(u, t, integrator) = u[2] > 0
function term_affect!(integrator, u, p, ctx)
    vx_t = integrator.u[u.vx]
    foo = sqrt(p.mx) * vx_t
    if foo > 0 && foo < vx_t
        @info "terminating at $(integrator.t)"
        terminate!(integrator)
    end
    nothing
end

function reflect_affect!(integrator, u, p, ctx)
    integrator.u[u.vx] = -integrator.u[u.vx]
    global N_COLLISIONS += 1
    term_affect!(integrator, u, p, ctx)
    nothing
end

continuous_events = [
    [x ~ 0] => (reflect_affect!, [vx, vy], [mx, my], nothing),
    [x ~ y, y ~ x] => (affect!, [vx, vy, x, y], [mx, my], nothing)
    # [y ~ 10] => [vy ~ -vy]
]

global N_COLLISIONS = 0
@named elastic = ODESystem(eqs, t, sts, [mx, my]; continuous_events)
old_sts = states(elastic)
sys = structural_simplify(elastic)
tspan = (0, 1e3)
prob = ODEProblem(sys, [], tspan; saveat=0.1)
sol = solve(prob)
# plot(sol)
nsteps = 500
saves = range(tspan..., length=nsteps)
solitp = sol(saves)
N_COLLISIONS

anim = @animate for i in 1:length(solitp)
    cur_vx = solitp[vx][i]
    cur_vy = solitp[vy][i]
    last_ten = maximum([1, i - 10]):i
    nlast = length(last_ten)
    # ncs = round(Int, solitp[n_collisions][i])
    curt = solitp.t[i]
    plot_t = round(Int, curt)
    myplot = plot(solitp[x][last_ten], zeros(nlast); title="t = $(plot_t)", xlims=(0, 10), ylims=(-5, 5))#, $([solitp.u[1, i], solitp.u[2, i]])")
    scatter!(myplot, [solitp[x][i]], [0.0]; ms=log(Mx) + 1)#, texts=["$cur_vx"])
    plot!(myplot, solitp[y][last_ten], zeros(nlast))
    scatter!(myplot, [solitp[y][i]], [0.0]; ms=log(My) + 1)#, texts=["$cur_vy"])
end
mp4(anim, "elastic_with_radius_$(My).mp4", fps=60)
@info "anim"

