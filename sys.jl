begin
    var"##iv#273" = (@variables(t))[1]
    var"##sts#274" = (collect)(@variables(xˍt(t), y(t), z(t), x(t)))
    var"##ps#275" = (collect)(@parameters(σ, ρ, β))
    var"##eqs#276" = [
        (Differential(t))(xˍt) ~ σ * (-1x + y)
        (Differential(t))(y) ~ (ρ + -1z) * x + -1y
        (Differential(t))(z) ~ x * y + -1 * β * z
        (Differential(t))(x) ~ xˍt
    ]
    var"##defs#277" = (Dict)()
    var"##iv#278" = (@variables(t))[1]
    (ODESystem)(
        var"##eqs#276",
        var"##iv#278",
        var"##sts#274",
        var"##ps#275";
        defaults = var"##defs#277",
        name = :sys,
        checks = false,
    )
end