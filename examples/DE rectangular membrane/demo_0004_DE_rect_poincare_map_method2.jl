using Comodo.GLMakie
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
Jm = 100.0
Ω = 1.58
Vdc = 0.1      # (ε*Φdc^2)/(μ*L3^2)
Vac= 0.1       # (Φac)/(Φdc)
P = 0.5        # S/μ
c = 0.0
function duffing!(dx, x, p, t)
    dx[1] = x[2]
    stiffness = Jm*(x[1] - x[1]^-5) /(Jm - 2*x[1]^2 - x[1]^-4 + 3)
    dx[2] = -c*x[2] + P + Vdc*(1 + Vac*sin(Ω*t))^2*x[1]^3 - stiffness
end

x0 = [1.149; 0.0]
tspan = (0.0, 300.0)
prob = DE.ODEProblem(duffing!, x0, tspan)
sol = DE.solve(prob,DE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
# solution = hcat(sol.u...)
t = sol.t
disp = sol[1, :]
velo = sol[2, :]
fig = Figure(size = (600,400))


T = 2π / Ω
steady_start = 0.0
n_start = ceil(Int, steady_start / T)
n_end = floor(Int, tspan[2] / T)
sample_times = T .* (n_start:n_end)

disp_poincare = [sol(t)[1] for t in sample_times]
velo_poincare = [sol(t)[2] for t in sample_times]
ax = Axis(fig[1, 1], xlabel=L"\lambda", ylabel=L"d\lambda/d t")
scatter!(ax, disp_poincare,velo_poincare, markersize= 10)
xlims!(ax, [0.0 4.0])
ylims!(ax,[-2.0 2.0])
display(GLMakie.Screen(), fig)