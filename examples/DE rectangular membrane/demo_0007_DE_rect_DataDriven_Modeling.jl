using DataDrivenDiffEq
using ModelingToolkit
using DifferentialEquations
using LinearAlgebra
using DataDrivenSparse
using Plots
#gr()
plotly() 
function DielectricElastomer(x, p, t)
    dx = similar(x)  
    Jm = 100.0
    Ω = 1.58
    Vdc = 0.1
    Vac = 0.1
    P = 0.
    c = 0.2
    dx[1] = x[2]
    stiffness = Jm * (x[1] - x[1]^-5) / (Jm - 2 * x[1]^2 - x[1]^-4 + 3)
    dx[2] = -c * x[2] + P + Vdc * (1 + Vac * sin(Ω * t))^2 * x[1]^3 - stiffness
    return dx
end

u0 = [1.0; 0.0]
tspan = (0.0, 10000.0)
dt = 0.01
prob = ODEProblem(DielectricElastomer, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12, saveat=dt)
plot(sol)
X = sol[:, :]
DX = similar(X)
for (i, xi) in enumerate(eachcol(X))
    DX[:, i] = DielectricElastomer(xi, [], sol.t[i])
end
t = sol.t
ddprob = ContinuousDataDrivenProblem(X, t, DX=DX)
@parameters t
@variables u[1:2]
u = collect(u)
Ψ = Basis([u; u[1]^2;u[1]^3; (0.1 + 0.1 * sin.(1.58 * t))], u, independent_variable=t)
λ = 1e-3
opt = STLSQ(λ)
res = solve(ddprob, Ψ, opt)

basis_result = get_basis(res)

pm = get_parameter_map(basis_result)

discovered_basis = get_basis(res)
p_discovered = get_parameter_values(discovered_basis)
discovered_prob = ODEProblem(discovered_basis, u0, tspan, p_discovered)
# Solve the discovered ODE
discovered_sol = solve(discovered_prob, Tsit5(),reltol = 1e-12, abstol = 1e-12, saveat= dt)

# Time-history plot
p1 = plot(sol, idxs=1, label="Original x(t)", linestyle=:dash, linewidth = 3)
plot!(p1, discovered_sol, idxs=1, label="Discovered x(t)", xlimits=(9800,10000))
xlabel!("Time")
ylabel!("x")
title!("Time History Comparison")