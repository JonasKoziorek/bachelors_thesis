using DynamicalSystems
using LinearAlgebra
using ControlSystems

dx = Systems.lorenz([-8.0, 8.0, 27.0]; σ = 10.0, ρ = 28.0, β = 8/3)

time = 200.0
step = 0.01
xdat,t = trajectory(dx, time-step; Δt = step)

function embed2(data, stackmax)
    H = zeros(stackmax,length(data)-stackmax);
    for k = 1:stackmax
        H[k,:] = data[k:end-stackmax-1+k]
    end  
    return H
end

stackmax = 100
H = embed2(xdat[:,1], stackmax);

U,S,V = svd(H)