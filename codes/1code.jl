using DynamicalSystems
using LinearAlgebra
using ControlSystems

dx = Systems.lorenz([-8.0, 8.0, 27.0]; σ = 10.0, ρ = 28.0, β = 8/3)

time = 200.0
step = 0.01
xdat = trajectory(dx, time-step; Δt = step).data

x = [a for (a,b,c) in xdat];
function embed2(data, stackmax)
    H = zeros(stackmax,length(x)-stackmax);
    for k = 1:stackmax
        H[k,:] = x[k:end-stackmax-1+k]
    end  
    return H
end

stackmax = 100
H = embed2(x, stackmax);

U,S,V = svd(H)