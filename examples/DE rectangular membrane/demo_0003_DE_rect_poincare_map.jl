using DynamicalSystems
using GLMakie
import DifferentialEquations as DE
GLMakie.closeall()
#=
Results for the paper:
Nonlinear dynamic characteristics of a dielectric elastomer membrane undergoing in-plane deformation
Sheng et al. 2014

Equation: 20

t :time
disp: displacement
vel: veclocity
=#
# Parameters


function DE!(du, u, p, t)
    (Jm, Ω , Vdc, Vac, P, c) = p 
    du[1] = u[2]
    stiffness = Jm * (u[1] - u[1]^(-5)) / (Jm - 2*u[1]^2 - u[1]^(-4) + 3)
    du[2] = -c * u[2] + P + Vdc*(1 + Vac*sin(Ω*t))^2 * u[1]^3 - stiffness
end

Jm = 100.0
Ω = 1.58
Vdc = 0.1      # (ε*Φdc^2)/(μ*L3^2)
Vac= 0.1       # (Φac)/(Φdc)
P = 0.5        # S/μ
c = 0.0
p = [Jm, Ω , Vdc, Vac, P, c]

T = 2π/Ω
x0 = [1.149, 0.0]
tmax = 300.0

### tol decrease simulation run time increae ## tolerance be fixed
diffeq = (alg = DE.Tsit5(), abstol = 1e-10, reltol = 1e-10)
# Define system
ds = ContinuousDynamicalSystem(DE!, x0, p ; diffeq)

# Compute trajectory using Δt = forcing period
Y, t = trajectory(ds, tmax; Δt=T, Ttr= 0.0)

# The Poincaré map is simply all trajectory points (sampled once per period)
poincare_points = Y  # each row corresponds to one period

# Plot Poincaré map
# scatter(poincare_points[:,1], poincare_points[:,2], color=:red, markersize=5)
#scatter(Y, markersize= 10)

fig = Figure(size = (600,400))
ax = Axis(fig[1, 1], xlabel=L"\lambda", ylabel=L"d\lambda/d t")
scatter!(ax, Y, color = :black, markersize= 10)
xlims!(ax, [0.0 5.0])
ylims!(ax,[-3.0 3.0])
display(GLMakie.Screen(), fig)