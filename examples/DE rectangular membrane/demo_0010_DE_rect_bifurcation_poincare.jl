using DifferentialEquations
using GLMakie
GLMakie.closeall()

# -----------------------------
# ODE system (unchanged)
# -----------------------------
function DE!(du, u, p, t)
    (Jm, Ω, Vdc, Vac, P, c) = p
    du[1] = u[2]
    stiffness = Jm * (u[1] - u[1]^(-5)) / (Jm - 2 * u[1]^2 - u[1]^(-4) + 3)
    du[2] = -c * u[2] + P + Vdc * (1 + Vac * sin(Ω * t))^2 * u[1]^3 - stiffness
end

# -----------------------------
# Parameters
# -----------------------------
Jm  = 100.0
Ω   = 1.0
Vdc = 0.3
P   = 0.5
c = 0.1

u0 = [1.88, 0.0]

Vac_values = 0.0:0.001:0.5

# -----------------------------
# Time settings
# -----------------------------
T     = 2π / Ω
tspan = (0.0, 1000.0)
ts    = 0.0:T:tspan[2]

Ttr  = 800.0     # transient time
Npts = 100       # points per Vac after transient

# -----------------------------
# Storage
# -----------------------------
bif_points = Float64[]
Vac_plot   = Float64[]

# -----------------------------
# Main loop
# -----------------------------
for Vac in Vac_values

    p = (Jm, Ω, Vdc, Vac, P, c)
    prob = ODEProblem(DE!, u0, tspan, p)

    sol = solve( prob, Tsit5(), reltol = 1e-8, abstol = 1e-8,  maxiters = 10^8, saveat = ts)

    idx = findall(t -> t > Ttr, sol.t)
    n = min(Npts, length(idx))
    idx = idx[end-n+1:end]

    for i in idx
        push!(bif_points, sol[1, i])
        push!(Vac_plot, Vac)
    end
    u0 = sol.u[end] |> collect
end

# -----------------------------
# Plot
# -----------------------------

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel = "Vac", ylabel = L"$\lambda$ (Poincaré map)")
scatter!(ax, Vac_plot, bif_points; markersize = 2, color = :black)
DataInspector(fig) 
display(fig)
