using DynamicalSystems, Comodo.GLMakie
import DifferentialEquations as DE
# -----------------------------
# Parameters
# -----------------------------
Jm = 100.0
Ω = 1.58
Vdc = 0.1
Vac = 0.5
P = 0.5

x0 = [1.4, 0.0]        # initial condition
T = 2π / Ω             # forcing period
tmax = 1000.0          # total integration time
Ttr = 200.0            # transient to discard
nplot = 500            # number of Poincaré points to plot after transient

# Range of damping parameter c
c_values = 0.0:0.001:0.2

# Store bifurcation points
bif_points = Float64[]

c_plot = Float64[]
diffeq = (alg = DE.Tsit5(), abstol = 1e-12, reltol = 1e-12)
# -----------------------------
# Loop over c
# -----------------------------
for c in c_values
    # define ODE system with current c
    function DE!(du, u, p, t)
        du[1] = u[2]
        stiffness = Jm * (u[1] - u[1]^(-5)) / (Jm - 2 * u[1]^2 - u[1]^(-4) + 3)
        du[2] = -c * u[2] + P + Vdc * (1 + Vac * sin(Ω * t))^2 * u[1]^3 - stiffness
    end

    # continuous system
    ds = ContinuousDynamicalSystem(DE!, x0; diffeq)

    # trajectory
    Y, t = trajectory(ds, tmax; Δt=T, Ttr=Ttr)

    npoints = min(nplot, size(Y, 1))
    start_idx = size(Y, 1) - npoints + 1  # compute start index
    for i = start_idx:size(Y, 1)
        push!(bif_points, Y[i, 1])   # displacement x
        push!(c_plot, c)
    end

end

@info "Done."
# -----------------------------
# Plot bifurcation diagram
# -----------------------------
GLMakie.closeall()
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], xlabel="c", ylabel="x (Poincaré map)")
scatter!(ax, c_plot, bif_points, markersize=2, color=:black)
display(fig)
