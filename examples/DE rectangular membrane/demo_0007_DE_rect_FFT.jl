using GLMakie
import DifferentialEquations as DE
using FFTW
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
function DE!(dx, x, p, t)
     (Jm, Ω , Vdc, Vac, P, c) = p
    dx[1] = x[2]
    stiffness = Jm*(x[1] - x[1]^-5) /(Jm - 2*x[1]^2 - x[1]^-4 + 3)
    dx[2] = -c*x[2] + P + Vdc*(1 + Vac*sin(Ω*t))^2*x[1]^3 - stiffness
end

Jm = 100.0
Ω = 1.58
Vdc = 0.1      # (ε*Φdc^2)/(μ*L3^2)
Vac= 0.1       # (Φac)/(Φdc)
P = 0.5        # S/μ
c = 0.0
p = [Jm, Ω , Vdc, Vac, P, c]

x0 = [1.149; 0.0]
tspan = (0.0, 5000.0)
prob = DE.ODEProblem(DE!, x0, tspan, p)

Ts = 0.1 ## sampling time
sol = DE.solve(prob,DE.Tsit5(), reltol = 1e-10, abstol = 1e-10, saveat = Ts)
t = sol.t
disp = sol[1, :]
velo = sol[2, :]
#####------------------
fig = Figure(size = (800,600))
ax1 = Axis(fig[1, 1], xlabel=L"t", ylabel=L"\lambda")
lines!(ax1, t, disp, color = :black, linewidth = 1)
xlims!(ax1, 4700, 5000.0)
#####------------------
# ## Hanning Window (to avoid frequency leakage)
L = length(disp) ## length of signal
nn =1:L
W0=@. 1/2*(1-cos(2*pi*(nn-1)/L)); ## window
disp=W0.*disp;

Y = fft(disp)
P2 = abs.(Y ./ L)
P1 = P2[1:div(L, 2)+1]  # integer division ensures an integer index
P1[2:end-1] .= 2 .* P1[2:end-1]
Fs = 1/Ts;
f = Fs/L*(0:(L/2))
ω =  2*π*f # circular frequency

####------------------
ax2 = Axis(fig[2, 1], xlabel="Frequency", ylabel=L"|\lambda(\omega)|")
lines!(ax2, ω, P1, color = :black, linewidth = 1)
xlims!(ax2, 0.01, 3.6)
ylims!(ax2, 0.0, 0.3)

#####------------------
ax3 = Axis(fig[2, 2], xlabel=L"\omega", ylabel= L"$|\lambda(\omega)|$ [dB]")

lines!(ax3, ω, 20*log.(P1), color = :black, linewidth = 1)
xlims!(ax3, 0.1, 3.0)
#ylims!(ax3, 0.0, 1000.0)

DataInspector(fig) 
display(GLMakie.Screen(), fig)

