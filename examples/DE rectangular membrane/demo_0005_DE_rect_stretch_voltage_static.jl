using BifurcationKit
using Comodo.GLMakie
#=
Result for 
Method to Control Dynamic Snap-Through Instability of Dielectric Elastomers
DOI: 10.1103/PhysRevApplied.6.064012
Figure 2
=#
GLMakie.closeall()
# function to record information from a solution
recordFromSolution(x, p; k...) = (u1 = x[1])
function DE!(dz, z, p, t=0)
    ## p1 == Vdc
    ## p2 == Jm
    ## p3 == P
    (; p1, p2, p3) = p
    u1, u2 = z
    dz[1] = u2
    dz[2] = -p2 * (u1 - u1^-5) / (p2 - 2 * u1^2 - u1^-4 + 3) + p1^2 * u1^3 + p3
    dz
end

par_DE = (p1=0., p2=100.0, p3=1.0)
z0 = [1.0, 0.0]
prob = BifurcationProblem(DE!, z0, par_DE,
    # specify the continuation parameter
    (@optic _.p1), record_from_solution=recordFromSolution)

opts_br = ContinuationPar(p_min=0., p_max=1.0, dsmax=0.001, dsmin=1e-6, ds=0.001, max_steps=5000,
    newton_options=NewtonPar(tol=1e-9, max_iterations=100))

br = continuation(prob, PALC(), opts_br; bothside=true)

ps = br.param
xs = br.x
stable = br.stable
# Split branch into contiguous stable/unstable segments
function contiguous_segments(ps, xs, stable)
    segs = []
    start = 1
    for i in 2:length(ps)
        if stable[i] != stable[i-1]
            push!(segs, (ps[start:i-1], xs[start:i-1], stable[i-1]))
            start = i
        end
    end
    push!(segs, (ps[start:end], xs[start:end], stable[end]))
    return segs
end

segs = contiguous_segments(ps, xs, stable)

# GLMakie plot
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1]; xlabel="λ", ylabel="Vdc") # the plot style that is common in DE (Dielectric Extractlastomers)
ax2 = Axis(fig[1, 2]; xlabel="Vdc", ylabel="λ") # the plot style that is common in MEMS (micro electromechanical systems)

function plot_first(segs)
    first_stable = true
    first_unstable = true

    for (param_segment, sol_segment, is_stable) in segs
        if is_stable
            line_color = :black
            line_style = :solid
            linewidth = 2
            label = first_stable ? "Stable" : nothing
            first_stable = false
        else
            line_color = :blue
            line_style = :dash
            linewidth = 2
            label = first_unstable ? "Unstable" : nothing
            first_unstable = false
        end

        lines!(ax, sol_segment, param_segment;
            color=line_color,
            linestyle=line_style,
            linewidth=linewidth,
            label=label)
    end
end

plot_(segs)
# Extract branch points (bp) only
# How to find the type:
# # Iterate over special points
# for sp in br.specialpoint
#     println("Index: ", sp.idx, ", Type: ", sp.type, ", Status: ", sp.status)
# end
bps = filter(sp -> sp.type == :bp, br.specialpoint)

if !isempty(bps)
    x_bp = [ps[sp.idx] for sp in bps]
    y_bp = [xs[sp.idx] for sp in bps]
    scatter!(ax, y_bp, x_bp; color=:red, markersize=12, label="Branch Point (bp)")
end

function plot_second(segs)
    first_stable = true
    first_unstable = true
    for (param_segment, sol_segment, is_stable) in segs
        if is_stable
            line_color = :black
            line_style = :solid
            linewidth = 2
            label = first_stable ? "Stable" : nothing
            first_stable = false
        else
            line_color = :blue
            line_style = :dash
            linewidth = 2
            label = first_unstable ? "Unstable" : nothing
            first_unstable = false
        end

        lines!(ax2, param_segment, sol_segment;
            color=line_color,
            linestyle=line_style,
            linewidth=linewidth,
            label=label)

    end
end
plot_second(segs)
# Extract branch points (bp) only
bps = filter(sp -> sp.type == :bp, br.specialpoint)

if !isempty(bps)
    x_bp = [ps[sp.idx] for sp in bps]
    y_bp = [xs[sp.idx] for sp in bps]
    scatter!(ax2, x_bp, y_bp; color=:red, markersize=12, label="Branch Point (bp)")
end


axislegend(ax; position=:rb)
display(GLMakie.Screen(), fig)