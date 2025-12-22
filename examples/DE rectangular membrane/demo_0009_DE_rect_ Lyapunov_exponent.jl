using DynamicalSystems
using GLMakie
import DifferentialEquations as DE
GLMakie.closeall()
# -----------------------------
# ODE system
# -----------------------------
function DE!(du, u, p, t)
    (Jm, Ω, Vdc, Vac, P, c) = p
    du[1] = u[2]
    stiffness = Jm * (u[1] - u[1]^(-5)) / (Jm - 2*u[1]^2 - u[1]^(-4) + 3)
    du[2] = -c*u[2] + P + Vdc*(1 + Vac*sin(Ω*t))^2 * u[1]^3 - stiffness
end

# -----------------------------
# Initial conditions & constants
# -----------------------------
u0 = [1.88, 0.0]
Jm = 100.0
Ω = 1.0
Vdc = 0.3
c = 0.1
P = 0.5

T = 0.1
# -----------------------------
# Parameter scan over c
# -----------------------------
Vac_values = 0.0:0.001:0.5
Lya = zeros(length(Vac_values))

diffeq = (alg = DE.Tsit5(), abstol = 1e-10, reltol = 1e-10)
@info "Computing Lyapunov exponents..."

for (i, Vac) in enumerate(Vac_values)
    p = (Jm, Ω, Vdc, Vac, P, c)
    ds = ContinuousDynamicalSystem(DE!, u0, p; diffeq)
    Lya[i] = lyapunov(ds, 5000.0; Δt=T, Ttr = 4800.0)
end

@info "Done."


# plot
fig = Figure(size = (800,600))
ax1 = Axis(fig[1, 1], xlabel="c", ylabel="Lyapunov exponent")
lines!(ax1, Vac_values, Lya, linewidth=2)
#ylims!(ax1, -.002, 0.02)
DataInspector(fig) 
display(fig)
