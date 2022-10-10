using ModelingToolkit, OrdinaryDiffEq, Symbolics
using Symbolics: solve_for, solve_single_eq
using Plots
using DifferentialEquations

# 1d elastic collision with 2 masses
Mx = 1
My = 10
@parameters t rx = 0.5 ry = 0.5 mx = Mx my = My wall1 = 0 wall2 = 10
sts = @variables x(t) = 2.0 y(t) = 5.0 vx(t) = 0.0 vy(t) = -1.0 n_collisions(t)=0
D = Differential(t)
eqs = [
    D(y) ~ vy,
    D(vy) ~ 0,
    D(x) ~ vx,
    D(vx) ~ 0,
    D(n_collisions) ~ 0
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
    integ.u[u.n_collisions] += 1
    nothing
end

continuous_events = [
    [(x - rx) ~ 0] => [vx ~ -vx, n_collisions ~ n_collisions + 1],
    [(x + rx) ~ (y - ry), (y - ry) ~ (x + rx)] => (affect!, [vx, vy, n_collisions], [mx, my], nothing),
    [(y + ry) ~ 10] => [vy ~ -vy, n_collisions ~ n_collisions + 1]
]

@named elastic = ODESystem(eqs, t, sts, [mx, my, rx, ry]; continuous_events)
old_sts = states(elastic)
sys = structural_simplify(elastic)
tspan = (0, 50)
prob = ODEProblem(sys, [], tspan; saveat=0.1)
sol = solve(prob)
# sol = solve(prob, Rosenbrock23())
# plot(sol)
nsteps = 500
saves = range(tspan..., length=nsteps)
solitp = sol(saves)

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

mp4(anim, "elastic_with_radius_and_counter_$(My).mp4", fps=60)
