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
function animated_draggable_pendulum(n_pendulums=1)
    # Create pendulum based on number requested
    if n_pendulums == 1
        pen = Pendulum(1; L=[1.5], θ=[π/4], ω=[0.0])
    elseif n_pendulums == 2
        pen = Pendulum(2; L=[1.0, 1.0], θ=[π/4, π/6], ω=[0.0, 0.0])
    elseif n_pendulums == 3
        pen = Pendulum(3; L=[0.8, 0.8, 0.8], θ=[π/4, π/6, π/8], ω=[0.0, 0.0, 0.0])
    else
        error("Only 1, 2, or 3 pendulums supported")
    end
    
    # Observables
    x, y = calculate_pendulum_positions(pen)
    
    # Create line segments connecting all pendulum points
    line_points = Observable([Point2f(x[i], y[i]) for i in 1:length(x)])
    # Mass points (all pendulum bobs)
    mass_points = Observable([Point2f(x[i], y[i]) for i in 2:length(x)])
    # Anchor point
    anchor_point = Observable([Point2f(x[1], y[1])])
    
    is_animating = Observable(false)
    animation_task = Ref{Union{Task, Nothing}}(nothing)
    
    # Create figure with padding around the entire figure
    fig = Figure(size=(600, 650), figure_padding=20)

    # Clean axis without decorations
    ax = Axis(fig[1, 1:2], aspect=1, limits=(-2.5, 2.5, -2.5, 2.5),
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false,
              leftspinevisible=false, rightspinevisible=false,
              topspinevisible=false, bottomspinevisible=false)
    
    # Buttons in a separate row
    button_row = GridLayout(fig[2, 1:2])
    play_btn = Button(button_row[1, 1], label="Play", width=100)
    pause_btn = Button(button_row[1, 2], label="Pause", width=100)
    
    # Layout settings
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    
    # Make axis huge, buttons tiny
    rowsize!(fig.layout, 1, Relative(1.0))
    rowsize!(fig.layout, 2, 50)

    colsize!(fig.layout, 1, Relative(1.0))
    colsize!(fig.layout, 2, Relative(1.0))

    # Disable interactions
    deactivate_interaction!(ax, :rectanglezoom)
    deactivate_interaction!(ax, :limitreset) 
    deactivate_interaction!(ax, :scrollzoom)
    deactivate_interaction!(ax, :dragpan)
    
    # Plot elements - lines connecting all points
    lines!(ax, line_points, color=:black, linewidth=4)
    # Anchor point
    scatter!(ax, anchor_point, color=:red, markersize=15)
    # All pendulum masses
    scatter!(ax, mass_points, color=:blue, markersize=25)
    
    # Dragging state
    dragging = Ref(false)
    dragging_index = Ref(0)  # Which pendulum mass is being dragged
    
    # Function to update the visual
    function update_visual!()
        x_new, y_new = calculate_pendulum_positions(pen)
        line_points[] = [Point2f(x_new[i], y_new[i]) for i in 1:length(x_new)]
        mass_points[] = [Point2f(x_new[i], y_new[i]) for i in 2:length(x_new)]
        return x_new, y_new
    end
    
    # Mouse dragging (only when not animating) - can drag any pendulum mass
    on(events(fig).mousebutton) do event
        if !is_animating[] && event.button == Mouse.left
            if event.action == Mouse.press
                mouse_pos = mouseposition(ax.scene)
                # Check which pendulum mass is being clicked
                for i in 1:length(mass_points[])
                    distance = norm(mouse_pos - mass_points[][i])
                    if distance < 0.3
                        dragging[] = true
                        dragging_index[] = i  # Which pendulum mass (1-indexed)
                        break
                    end
                end
            else
                dragging[] = false
                dragging_index[] = 0
            end
        end
    end
    
    on(events(fig).mouseposition) do mouse_pos
        if !is_animating[] && dragging[] && dragging_index[] > 0
            world_pos = mouseposition(ax.scene)
            
            # Calculate angle for the dragged pendulum segment
            x_curr, y_curr = calculate_pendulum_positions(pen)
            pendulum_idx = dragging_index[]  # Which pendulum we're dragging
            
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
            pen.ω[pendulum_idx] = 0.0
            
            x, y = update_visual!()
        end
    end
    
    # Play button functionality
    on(play_btn.clicks) do n
        if is_animating[]
            # Stop animation
            is_animating[] = false
            play_btn.label = "Play"
        else
            # Start animation
            is_animating[] = true
            play_btn.label = "Stop"
            
            # Start animation loop
            animation_task[] = @async begin
                while is_animating[]
                    pen()  # Update physics using your RK4 code!
                    x, y = update_visual!()
                    sleep(0.016)  # ~60 FPS
                end
            end
        end
    end
    
    # Pause button - stops motion but keeps current position
    on(pause_btn.clicks) do n
        # Stop all velocities
        pen.ω .= 0.0
        # Stop animation if running
        is_animating[] = false
        play_btn.label = "Play"
    end
    
    return fig
end