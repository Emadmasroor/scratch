using DifferentialEquations, NLsolve, Plots

include("functions6.jl")

global ϵ = 0.01

# Set up initial condition
global switchΔt = 0.2

global lefttrapped = false
global righttrapped = false
global stepnum = 0

# Arrays which will contain the stitched-together full arrays
global rr = Array{Float64,1}(undef,0)
global θθ = Array{Float64,1}(undef,0)
global tt = Array{Float64,1}(undef,0)
global HH = Array{Float64,1}(undef,0)

# Set up callbacks for sink-eject. These won't change.
cb1b = DiscreteCallback(sinkRight,ejectRight!)
cb2b = DiscreteCallback(sinkLeft,ejectLeft!)

global u0 = [0.6,3π/4]
for i in 1:300
    # Declare all global variables
    global u0, switchΔt, stepnum, lefttrapped, righttrapped, sol1, sol2, cb1a,cb1b,cb2a,cb2b,ψvalue
    ##############---Right---####################
    stepnum+=1

    println("--->>>Integrating right-handed system, step ",stepnum)
    # Update the manifold we are on, only if we are not trapped on the right
    if ~righttrapped && ~lefttrapped
        #println("updating ψvalue on the right")
        ψvalue = ψ_right(u0[1],u0[2])
        cb1a = ManifoldProjection(energyManifold_right)
    end

    # Advance things forward for "switchΔt" time
    if ~lefttrapped
        prob1 = ODEProblem(ODErhs_right!,u0,(2*switchΔt*(i-1),2*switchΔt*(i-1)+switchΔt), callback=CallbackSet(cb1a,cb1b))
        sol1 = solve(prob1,Tsit5(),abstol=1e-8,reltol=1e-8)
        #println("Sent particle from ",short(sol1[:,1])," to ",short(sol1[:,end]))
        #println("Currently, ψvalue = ",short(ψvalue))
    else
        prob1 = ODEProblem(ODErhs_empty!,sol2[:,end],(2*switchΔt*(i-1),2*switchΔt*(i-1)+switchΔt))
        #println("!!! Particle trapped in left-pipe. Not integrating, sit and wait.")
        #println("Meanwhile, current value is (r,θ) = ( ",round.(sol2[:,end],digits=3),")")
        sol1 = solve(prob1,Euler(),dt=0.001)
    end
    # Save information into the relevant arrays
    append!(rr,sol1[1,:])
    append!(θθ,sol1[2,:])
    append!(tt,sol1.t)
    append!(HH,ψ_right.(sol1[1,:],sol1[2,:]))
    ##############---Left----####################
    stepnum+=1
    println("--->>>Integrating left-handed system, step ",stepnum)
    # Decide the new manifold to stay on:
    if ~lefttrapped && ~righttrapped
        # This code will only be executed if the system wasn't initially trapped.
        # If it was trapped, we would like to delay the execution of this snippet
        # of code until the particle has emerged.
        #println("updating ψvalue on the left")
        ψvalue = ψ_left(sol1[1,end],sol1[2,end])
        cb2a = ManifoldProjection(energyManifold_left)
    end

    # Advance things forward for "switchΔt" time
    if ~righttrapped
        prob2 = ODEProblem(ODErhs_left!,sol1[:,end],(2*switchΔt*i-switchΔt,2*switchΔt*i),callback=CallbackSet(cb2a,cb2b))
        sol2 = solve(prob2,Tsit5(),abstol=1e-8,reltol=1e-8)
        #println("Sent particle from ",short(sol2[:,1])," to ",short(sol2[:,end]))
        #println("Currently, ψvalue = ",short(ψvalue))
    else
        prob2 = ODEProblem(ODErhs_empty!,sol1[:,end],(2*switchΔt*i-switchΔt,2*switchΔt*i))
        #println("!!! Particle trapped in right-pipe. Not integrating, sit and wait.")
        #println("Meanwhile, current value is (r,θ) = ( ",short(sol1[:,end]),")")
        sol2 = solve(prob2,Euler(),dt=0.001)
    end
    # Save information into the relevant arrays
    append!(rr,sol2[1,:])
    append!(θθ,sol2[2,:])
    append!(tt,sol2.t)
    append!(HH,ψ_left.(sol2[1,:],sol2[2,:]))

    # Re-assign u0
    u0 = sol2[:,end]
end

# Plot trajectory
plot(θθ,rr,proj=:polar,
    lims=(0,1),
    linetype=:scatter,
    markersize = 2,
    legend=:none)
    savefig("Trajectory141steps.png")

# Which points are landing outside circle of radius 1?
outside=[(rr[i],tt[i], i) for i in 1:length(rr) if rr[i] > 1]

# Looks like there are four when I went back and forth 6 times.
# All these points occur right after ejection. And, since we pull ourselves back
# immediately after, due to the manifold projection, I think these points should
# not matter.

# Plot the ψ values
# These will necessarily be wrong whenever the particle is trapped
plot(tt,HH,
    linetype=:scatter,
    markersize=1,
    legend=:none)
