# using GLMakie
# import DifferentialEquations as DE
# GLMakie.closeall()

# function DE!(dx, x, p, t)
#     (Jm, Ω , Vdc, Vac, P, c) = p 
#     dx[1] = x[2]
#     stiffness = Jm*(x[1] - x[1]^-5) /(Jm - 2*x[1]^2 - x[1]^-4 + 3)
#     dx[2] = -c*x[2] + P + Vdc*(1 + Vac*sin(Ω*t))^2*x[1]^3 - stiffness
# end

# Jm = 100.0
# Ω = 1.58
# Vdc = 0.1      # (ε*Φdc^2)/(μ*L3^2)
# Vac= 0.1       # (Φac)/(Φdc)
# P = 0.5        # S/μ
# c = 0.0
# p = [Jm, Ω , Vdc, Vac, P, c]

# x0 = [1.149; 0.0]
# tspan = (0.0, 300.0)
# prob = DE.ODEProblem(DE!, x0, tspan, p)

# T = 2π/Ω
# ts = 0.0:T:tspan[2]

# sol = DE.solve(prob,DE.Tsit5(), reltol = 1e-10, abstol = 1e-10, saveat = ts)
# t = sol.t
# disp = sol[1, :]
# velo = sol[2, :]

# # visualisation 
# ## time history
# fig = Figure(size = (600,400))
# ax = Axis(fig[1, 1], xlabel=L"\lambda", ylabel=L"d\lambda/d t")
# # lines!(ax1, disp, velo, color = :blue, linewidth = 2)
# scatter!(ax, disp, velo)
# # display(GLMakie.Screen(), fig)
# xlims!(ax, [0.0 5.0])
# ylims!(ax,[-3.0 3.0])
# display(GLMakie.Screen(), fig)




using GLMakie
import DifferentialEquations as DE
GLMakie.closeall()

function DE!(dx, x, p, t)
    (Jm, Ω , Vdc, Vac, P, c) = p 
    dx[1] = x[2]
    stiffness = Jm*(x[1] - x[1]^-5) /(Jm - 2*x[1]^2 - x[1]^-4 + 3)
    dx[2] = -c*x[2] + P + Vdc*(1 + Vac*sin(Ω*t))^2*x[1]^3 - stiffness
end

Jm = 100.0
Ω = 1.58
Vdc = 0.1      # (ε*Φdc^2)/(μ*L3^2)
Vac = 0.1     # (Φac)/(Φdc)
P = 0.5        # S/μ
c = 0.0
p = [Jm, Ω , Vdc, Vac, P, c]

x0 = [1.149; 0.0]
tspan = (0.0, 300.0)
prob = DE.ODEProblem(DE!, x0, tspan, p)

T = 2π/Ω
ts = 0.0:T:tspan[2]

sol = DE.solve(prob,DE.Tsit5(), reltol = 1e-10, abstol = 1e-10, saveat = ts)
t = sol.t
disp = sol[1, :]
velo = sol[2, :]

Ttr = 0.0                       # transient time

idx = findall(t -> t > Ttr, sol.t)

disp_p = disp[idx]
velo_p = velo[idx]
# visualisation 
## time history
fig = Figure(size = (600,400))
ax = Axis(fig[1, 1], xlabel=L"\lambda", ylabel=L"d\lambda/d t")
# lines!(ax1, disp, velo, color = :blue, linewidth = 2)
scatter!(ax, disp_p, velo_p)
# display(GLMakie.Screen(), fig)
xlims!(ax, [0.0 5.0])
ylims!(ax,[-3.0 3.0])
display(GLMakie.Screen(), fig)