using DynamicalSystems
using Comodo.GLMakie
import DifferentialEquations as DE
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
u0 = [1.4, 0.0]
Jm = 100.0
Ω = 1.58
Vdc = 0.1
Vac = 0.5
P = 0.5

T = 2π / Ω  # Forcing period 
# -----------------------------
# Parameter scan over c
# -----------------------------
c_values = 0.0:0.001:0.2
Lya = zeros(length(c_values))

diffeq = (alg = DE.Tsit5(), abstol = 1e-8, reltol = 1e-8)
@info "Computing Lyapunov exponents..."

for (i, c) in enumerate(c_values)
    p = (Jm, Ω, Vdc, Vac, P, c)
    ds = ContinuousDynamicalSystem(DE!, u0, p; diffeq)

    Lya[i] = lyapunov(ds, 1000.0; Δt=T, Ttr = 600.0)

end

@info "Done."

GLMakie.closeall()
# plot
fig = Figure(size = (800,600))
ax1 = Axis(fig[1, 1], xlabel="c", ylabel="Lyapunov exponent")
lines!(ax1, c_values, Lya, linewidth=2)
#ylims!(ax1, -.002, 0.02)
DataInspector(fig) 
display(fig)
