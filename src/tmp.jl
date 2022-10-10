using Symbolics, DataInterpolations

@register Base.isless(x::AbstractFloat, y)
@register isless(x::AbstractFloat, y)
# @register isless(x::AbstractFloat, y::Num)

function foo(a, b, c)
    if isless(3.0, a) 
        b = c
    end
end

function foo2(a, b, c)
    if isless(3.0, a) 
        b = c
    end
end

function bar(a, b)
    u = rand(5)
    t = 0:4
    itp = LinearInterpolation(u,t)
    itp(a)
end

function bar2(a, b, c)
    itp = LinearInterpolation(b, c)
    itp(a)
end

@register foo(a, b, c)
@register bar(a, b::AbstractFloat)

@variables x
u = rand(5)
t = collect(0.:4)

@register bar2(x, u::Vector{AbstractFloat}, t::Vector{AbstractFloat})
@register bar2(x, u::AbstractVector, t::AbstractVector)

bar2(x, u, t)



foo(x, 4, 5)
foo2(x, 4, 5)
bar(x, 5.)

h(x, y) = x^2 + y
f(x, y) = h(x, y) + y
@register h(x, y)

