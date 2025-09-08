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
        x, y = calculate_pendulum_positions(pen)
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

using GLMakie, LinearAlgebra
using GLMakie: Button, scatter!
# Simulation struct to hold multiple pendulums
struct Simulation
    pendulums::Vector{Pendulum}
    colors::Vector{Symbol}
    
    function Simulation(pendulums::Vector{<:Pendulum})
        # Auto-assign colors
        default_colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown]
        colors = [default_colors[mod1(i, length(default_colors))] for i in 1:length(pendulums)]
        new(pendulums, colors)
    end
    
    function Simulation(pendulums::Vector{<:Pendulum}, colors::Vector{Symbol})
        @assert length(pendulums) == length(colors) "Number of pendulums must match number of colors"
        new(pendulums, colors)
    end
end

# Convenience constructors
Simulation(pen::Pendulum) = Simulation([pen])
Simulation(pens::Pendulum...) = Simulation([pens...])

function animated_multi_pendulum(sim::Simulation)
    # Calculate dynamic axis limits based on all pendulums
    max_reach = maximum(sum(pen.L) for pen in sim.pendulums) + 0.5
    axis_lim = max_reach + 0.2
    
    # Create observables for each pendulum
    pendulum_data = []
    for (i, pen) in enumerate(sim.pendulums)
        x, y = calculate_pendulum_positions(pen)
        
        line_points = Observable([Point2f(x[j], y[j]) for j in 1:length(x)])
        mass_points = Observable([Point2f(x[j], y[j]) for j in 2:length(x)])
        anchor_point = Observable([Point2f(x[1], y[1])])
        
        push!(pendulum_data, (
            pendulum = pen,
            color = sim.colors[i],
            line_points = line_points,
            mass_points = mass_points,
            anchor_point = anchor_point
        ))
    end
    
    is_animating = Observable(false)
    animation_task = Ref{Union{Task, Nothing}}(nothing)
    
    # Create figure
    fig = Figure(size=(700, 750), figure_padding=20)

    # Clean axis without decorations
    ax = Axis(fig[1, 1:2], aspect=1, limits=(-axis_lim, axis_lim, -axis_lim, axis_lim),
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false,
              leftspinevisible=false, rightspinevisible=false,
              topspinevisible=false, bottomspinevisible=false)
    
    # Buttons
    button_row = GridLayout(fig[2, 1:2])
    play_btn = Button(button_row[1, 1], label="Play", width=100)
    pause_btn = Button(button_row[1, 2], label="Pause", width=100)
    
    # Layout settings
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    rowsize!(fig.layout, 1, Relative(1.0))
    rowsize!(fig.layout, 2, 50)
    colsize!(fig.layout, 1, Relative(1.0))
    colsize!(fig.layout, 2, Relative(1.0))

    # Disable interactions
    deactivate_interaction!(ax, :rectanglezoom)
    deactivate_interaction!(ax, :limitreset) 
    deactivate_interaction!(ax, :scrollzoom)
    deactivate_interaction!(ax, :dragpan)
    
    # Plot all pendulums with their colors
    for data in pendulum_data
        lines!(ax, data.line_points, color=data.color, linewidth=4)
        scatter!(ax, data.anchor_point, color=:black, markersize=12)
        scatter!(ax, data.mass_points, color=data.color, markersize=20)
    end
    
    # Dragging state - using the EXACT same logic as the working single pendulum
    dragging = Ref(false)
    dragging_pendulum_idx = Ref(0)  # Which pendulum
    dragging_mass_idx = Ref(0)      # Which mass within that pendulum
    
    # Update function for a single pendulum (same as original)
    function update_single_pendulum_visual!(p_idx)
        data = pendulum_data[p_idx]
        x_new, y_new = calculate_pendulum_positions(data.pendulum)
        data.line_points[] = [Point2f(x_new[i], y_new[i]) for i in 1:length(x_new)]
        data.mass_points[] = [Point2f(x_new[i], y_new[i]) for i in 2:length(x_new)]
        return x_new, y_new
    end
    
    # Update all pendulums
    function update_all_visuals!()
        for i in 1:length(pendulum_data)
            update_single_pendulum_visual!(i)
        end
    end
    
    # Mouse dragging - EXACT same logic as working single pendulum
    on(events(fig).mousebutton) do event
        if !is_animating[] && event.button == Mouse.left
            if event.action == Mouse.press
                mouse_pos = mouseposition(ax.scene)
                # Check which pendulum mass is being clicked
                for (p_idx, data) in enumerate(pendulum_data)
                    for (m_idx, mass_pos) in enumerate(data.mass_points[])
                        distance = norm(mouse_pos - mass_pos)
                        if distance < 0.3
                            dragging[] = true
                            dragging_pendulum_idx[] = p_idx
                            dragging_mass_idx[] = m_idx
                            return
                        end
                    end
                end
            else
                dragging[] = false
                dragging_pendulum_idx[] = 0
                dragging_mass_idx[] = 0
            end
        end
    end
    
    # EXACT same mouse position logic as the working single pendulum
    on(events(fig).mouseposition) do mouse_pos
        if !is_animating[] && dragging[] && dragging_pendulum_idx[] > 0
            world_pos = mouseposition(ax.scene)
            
            # Get the pendulum and mass index
            pen = sim.pendulums[dragging_pendulum_idx[]]
            pendulum_idx = dragging_mass_idx[]  # Which mass we're dragging
            
            # Calculate angle for the dragged pendulum segment - EXACT same logic
            x_curr, y_curr = calculate_pendulum_positions(pen)
            
            # Calculate angle relative to the appropriate anchor point
            if pendulum_idx == 1
                # First pendulum - relative to fixed anchor
                anchor_x, anchor_y = x_curr[1], y_curr[1]
            else
                # Later pendulum - relative to previous mass position
                anchor_x, anchor_y = x_curr[pendulum_idx], y_curr[pendulum_idx]
            end
            
            dx = world_pos[1] - anchor_x
            dy = world_pos[2] - anchor_y
            
            pen.θ[pendulum_idx] = atan(dx, -dy)
            pen.ω[pendulum_idx] = 0.0  # Only zero THIS mass's velocity
            
            # Update only this pendulum's visual
            x, y = update_single_pendulum_visual!(dragging_pendulum_idx[])
        end
    end
    
    # Play button - animates ALL pendulums
    on(play_btn.clicks) do n
        if is_animating[]
            is_animating[] = false
            play_btn.label = "Play"
        else
            is_animating[] = true
            play_btn.label = "Stop"
            
            animation_task[] = @async begin
                while is_animating[]
                    for data in pendulum_data
                        data.pendulum()
                    end
                    update_all_visuals!()
                    sleep(0.016)
                end
            end
        end
    end
    
    # Pause button - stops ALL pendulums
    on(pause_btn.clicks) do n
        for data in pendulum_data
            data.pendulum.ω .= 0.0
        end
        is_animating[] = false
        play_btn.label = "Play"
    end
    
    return fig
end

# Example usage:
p1 = Pendulum(1; L=[1.5], θ=[π/4], ω=[0.0])
p2 = Pendulum(2; L=[1.0, 1.0], θ=[π/6, π/3], ω=[0.0, 0.0])
p3 = Pendulum(3; L=[0.8, 0.8, 0.8], θ=[π/8, π/4, π/2], ω=[0.0, 0.0, 0.0])

sim = Simulation([p1, p2, p3])
animated_multi_pendulum(sim)