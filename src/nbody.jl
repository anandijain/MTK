using ModelingToolkit, OrdinaryDiffEq, Symbolics
using LinearAlgebra, Statistics
using Unitful, Measurements
using Plots
using SparseArrays
using DataFrames
using DifferentialEquations
using Distributions
using WAV, Colors, ImageInTerminal
using TimerOutputs
using Images
using ImageTransformations

const GAUSSIAN = Normal(0, 1)

function random_point_on_sphere()
    x = rand(GAUSSIAN, 3)
    x .* 1 / sqrt(sum((x .^ 2)))
end
function mat_to_wav(arr, fn; sr=44100)
    for i in 1:size(arr, 2)
        arr[:, i] .= arr[:, i] ./ maximum(arr[:, i])
    end
    arr = 2 * arr .- 1

    for i in 1:2:size(arr, 2)-1
        fp = "$(fn)_$(i).wav"
        @info fp
        wavwrite(arr[:, i:i+1], sr, fp)
    end
    arr
end
# assumes parameters are constant
function my_getindex(sol, num)
    f = sol.prob.f
    hasproperty(f, :sys) || error()
    sol[substitute(Symbolics.scalarize(num), Dict(parameters(f.sys) .=> sol.prob.p))]
end

const to = TimerOutput()
N = 16
# for N in 10:10:50 
radius = 200
tspan = (0, 10000.0)
Npts = 10000

ps = collect((random_point_on_sphere() for i in 1:N))
initial_positions = radius * stack(ps)
masses = ones(N)

# scatter3d(PM[1, :], PM[2, :], PM[3, :], seriestype=:scatter, legend=false)

@parameters t M[1:N] = masses G = 10
@variables R(t)[1:3, 1:N] = initial_positions


D = Differential(t)
D2 = D^2
R = collect(R)
eqs = Equation[]
for i in 1:N
    Fi = zeros(3)
    for j in 1:N
        i == j && continue
        # the force on mass i due to mass j
        Fij = (G * M[i] * M[j] * (R[:, j] - R[:, i])) / norm(R[:, j] - R[:, i])^3
        Fi += Fij
    end

    feqs = D2.(R[:, i]) .~ Fi / M[i]
    append!(eqs, feqs)
end

# @timeit to "nbody $N cse true" nbody = ODESystem(eqs; cse=true, name=Symbol(:nbody, N))
@timeit to "nbody $N" nbody = ODESystem(eqs; name=Symbol(:nbody, N))
@info "sys"
old_sts = states(nbody)
# @info N 
# @timeit to "ode_lower $N" sys = ode_order_lowering(nbody) 
@timeit to "ss $N" ssys = structural_simplify(nbody)
@info "ss"
sys = ssys
# end

sts = states(sys)
new_sts = setdiff(sts, old_sts)
velocity_u0 = new_sts .=> 0.0

saveat = abs(-(tspan...) / Npts)
prob = ODEProblem(sys, velocity_u0, tspan; saveat)
@info "prob"

using OrdinaryDiffEq
const cfun = ODEFunction{true}(nbody);
prob = SecondOrderODEProblem((out, du, u, p, t) -> cfun.f(out, u, p, t), zeros(length(states(nbody))), map(v -> ModelingToolkit.defaults(nbody)[v], states(nbody)), (0, 5000.0), map(v -> ModelingToolkit.defaults(nbody)[v], parameters(nbody)));
sol = solve(prob, VelocityVerlet(), dt=1e-2);
sol = solve(prob, DPRKN12());

# sol = solve(prob)
@info "sol $(sol.retcode)"

# sol = solve(prob, Tsit5())
df = DataFrame(sol(range(0, 5000, length=30 * 5)))

# plot(sol)
# plot(sol[old_sts])
# plot(sol[R[:, 1]], vars=R[:, 1])

ts = df.timestamp
vel = df[:, 3:2:end]
pos = df[:, 2:2:end]

anim = @animate for i in 1:size(df, 1)
    myplot = Plots.plot(; title=df.timestamp[i])
    for j in 1:3:(3*N-2)
        plot3d!(myplot, pos[1:i, j], pos[1:i, j+1], pos[1:i, j+2], label="$j")
        scatter3d!(myplot, [pos[i, j]], [pos[i, j+1]], [pos[i, j+2]]; markersize=masses[ceil(Int, j // 3)])
    end
end

mp4(anim, "hi.mp4", fps=30)
@info "gif"

# wavwrite(arr[1:2, :]', 44100, "foo.wav")
pos_mat = Matrix(pos)
vel_mat = Matrix(vel)
# arr = Array(sol)
arr = pos_mat .+ abs(minimum(pos_mat))
arr2 = vel_mat .+ abs(minimum(vel_mat))

mat_to_wav(arr, joinpath("/Users/anand/Music/gen/", "foo"))
mat_to_wav(arr2, joinpath("/Users/anand/Music/gen/", "vel"))
@info "wav"
function df_to_pixarr(df)
    map(x -> RGB(x...), eachrow(df))
end

pos_dfs = [pos[:, i:i+2] for i in 1:3:ncol(pos)]
vel_dfs = [vel[:, i:i+2] for i in 1:3:ncol(vel)]
pic = map(df_to_pixarr, pos_dfs)
pics = map(df_to_pixarr, vel_dfs)


function build_frames(pics)
    @assert allequal(length.(pics))
    # frames = Array{RGB{Float64}, 3}(undef, size(pics[1], 1), size(pics[1], 2), length(pics))
    frames = []
    for i in eachindex(pics[1])
        f = reshape(getindex.(pics, i), 4, 4)
        push!(frames, f)
    end
    frames
end

fs = build_frames(pics)
fs = build_frames(pic)

save_pix_t = RGB{N0f8}
function convert_frame(f)
    f = map(clamp01nan, f)
    convert.(save_pix_t, f)
end
function convert_frames(fs; save_pix_t=RGB{N0f8})
    frames = Matrix{save_pix_t}[]
    for i in 1:length(fs)
        # frame = reshape(getindex.(pixs, i), (N, N))
        # frame = reshape(getindex.(pixs, i), img_size)
        push!(frames, convert_frame(fs[i]))
        # save("data/out3/$i.png", frame)
    end
    frames
end

using VideoIO
# encoder_options = (crf=23, preset="medium")
# VideoIO.save("video.mp4", frames, framerate=60, encoder_options=encoder_options)
# vid = VideoIO.load("video.mp4")
# vid2 = imresize.(vid, ratio=100)
# VideoIO.save("video.mp4", vid2, framerate=60)
# vid_name = "video"
# cmd = Cmd(["ffmpeg -framerate 60 -start_number 1 -i 'data/out4/%d.png' -r 60", "-y", "data/$(vid_name).mp4"])
# run(cmd)

function my_resize(img, amt)
    new_size = size(img) .* amt
    # img = imresize(img, ratio=ratio)
    img = convert_frame(img)
    new_img = zeros(RGB{N0f8}, new_size...)
    for i in 1:size(img, 1)
        for j in 1:size(img, 2)
            Binds = CartesianIndices(((i-1)*amt+1:i*amt, (j-1)*amt+1:j*amt))
            # @info Binds
            new_img[Binds] .= img[i, j]
        end
    end
    new_img
end

newfs = map(f -> my_resize(f, 100), fs)
VideoIO.save("wacky.mp4", newfs, framerate=60)

xwacky = [Array{RGB{N0f8},2}(undef, new_size...) for _ in 1:60*10]
wacky[3]
wackys = convert_frames(wacky)


pic = load("vel_pic.png")
pic[1, :]



# adding elastic collision

