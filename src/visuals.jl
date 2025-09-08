using Plots

function calculate_pendulum_positions(pen::Pendulum)    
    n = length(pen.M)+1; x = zeros(n); y = zeros(n)
    for i in range(2, n)
        x[i] = x[i-1] + pen.L[i-1] * sin(pen.θ[i-1])
        y[i] = y[i-1] - pen.L[i-1] * cos(pen.θ[i-1])
    end
    return x,y
end

function animate_pendulum(pen::Pendulum; dt=0.01, g=9.82, frames::Int=100, fps::Int=60, filename="pendulum.gif", giftitle="Pendulum Simulation", xlabel="X positions", ylabel="Y position")
    max_len = sum(pen.L)
    anim = @animate for i in 1:frames
        pen(; dt=dt, g=g)
        x, y = calculate_pendulum_positions(p)
        plot(x, y, 
            xlims=(-max_len,max_len), ylims=(-max_len,max_len), 
            aspect_ratio=:equal,
            line=(:black, 2),
            marker=(:circle, 6),
            title=giftitle,
            xlabel=xlabel, ylabel=ylabel)
    end
    gif(anim, filename, fps=fps)
end
