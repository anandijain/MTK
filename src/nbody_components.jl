using ModelingToolkit, Distributions, OrdinaryDiffEq
@parameters t G=10
D = Differential(t)
const GAUSSIAN = Normal(0, 1)

function random_point_on_n_sphere(n)
    x = rand(GAUSSIAN, n)
    x .* 1 / sqrt(sum((x .^ 2)))
end

function Mass(; name, m=1.0, spatial_dims=3, r0=random_point_on_n_sphere(spatial_dims), v0=random_point_on_n_sphere(spatial_dims))
    ps = @parameters m = m
    sts = @variables r(t)[1:spatial_dims] = r0 v(t)[1:spatial_dims] = v0
    eqs = Equation[]
    ODESystem(eqs, t, [r..., v...], ps; name)
end
N = 16
@named ms[1:N] = Mass()

function connect_masses()
    
end