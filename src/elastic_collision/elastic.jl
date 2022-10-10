using ModelingToolkit, OrdinaryDiffEq, Symbolics
using Symbolics: solve_for, solve_single_eq
using Plots
using DifferentialEquations

# 1d elastic collision with 2 masses
Mx = 1
My = 100
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
    nothing
end

continuous_events = [
    [x ~ 0] => [vx ~ -vx]
    [x ~ y, y ~ x] => (affect!, [vx, vy], [mx, my], nothing)
    [y ~ 10] => [vy ~ -vy]
]

@named elastic = ODESystem(eqs, t, sts, [mx, my]; continuous_events)
old_sts = states(elastic)
sys = structural_simplify(elastic)
tspan = (0, 50)
prob = ODEProblem(sys, [], tspan; saveat=0.1)
sol = solve(prob)
# plot(sol)

anim = @animate for i in 1:length(sol)
    cur_vx = sol[vx][i]
    cur_vy = sol[vy][i]
    myplot = plot(sol[x][1:i], zeros(i); title="$(sol.t[i])", xlims=(0, 10), ylims=(-5, 5))#, $([sol.u[1, i], sol.u[2, i]])")
    scatter!(myplot, [sol[x][i]], [0.0]; ms=Mx)#, texts=["$cur_vx"])
    plot!(myplot, sol[y][1:i], zeros(i))
    scatter!(myplot, [sol[y][i]], [0.0]; ms=My)#, texts=["$cur_vy"])
end
mp4(anim, "elastic.mp4", fps=60)

