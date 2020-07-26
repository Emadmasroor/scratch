###############################
### Right-hand-side of ODEs ###
###############################

function ODErhs_right!(Ẋ,X,p,t)
# This function takes as input a 2-vector (r,θ) andreturns a two-vector, (rdot, θdot).
# The symbolic functions for these two speeds have been calculated from Mathematica.
    r=X[1]
    θ=X[2]
    # Equations of motion
    Ẋ[1] = (-0.8*r*sin(θ)/sqrt(0.64 - 1.6*r*cos(θ) + r^2) -  0.8*(-r*sin(θ)/sqrt(0.64 - 1.6*r*cos(θ) + r^2) + 0.8*r*sin(θ)*(0.8 - r*cos(θ))/(0.64 - 1.6*r*cos(θ) + r^2)^(3/2)) - 1.0*(-r*sin(θ)/sqrt(1.5625 - 2.5*r*cos(θ) + r^2) + 1.25*r*sin(θ)*(1.25 - r*cos(θ))/(1.5625 - 2.5*r*cos(θ) + r^2)^(3/2)))/(r^2*sin(θ))
    Ẋ[2] = -(1.0*(1 + (-1/2)*(2*r - 1.6*cos(θ))/sqrt(0.64 - 1.6*r*cos(θ) + r^2)) - 0.8*(cos(θ)/sqrt(0.64 - 1.6*r*cos(θ) + r^2) + (1/2)*(0.8 - r*cos(θ))*(2*r - 1.6*cos(θ))/(0.64 - 1.6*r*cos(θ) + r^2)^(3/2)) - 1.0*(cos(θ)/sqrt(1.5625 - 2.5*r*cos(θ) + r^2) + (1/2)*(1.25 - r*cos(θ))*(2*r - 2.5*cos(θ))/(1.5625 - 2.5*r*cos(θ) + r^2)^(3/2)))/(r*sin(θ))
end

function ODErhs_left!(Ẋ,X,p,t)
# This function takes as input a 2-vector (r,θ) andreturns a two-vector, (rdot, θdot).
# The symbolic functions for these two speeds have been calculated from Mathematica.

    r=X[1]
    θ=X[2]
    # Equations of motion
    Ẋ[1] = -(0.8*r*sin(θ)/sqrt(0.64 + 1.6*r*cos(θ) + r^2) - 0.8*(r*sin(θ)/sqrt(0.64 + 1.6*r*cos(θ) + r^2) - 0.8*r*sin(θ)*(0.8 + r*cos(θ))/(0.64 + 1.6*r*cos(θ) + r^2)^(3/2)) - 1.0*(r*sin(θ)/sqrt(1.5625 + 2.5*r*cos(θ) + r^2) - 1.25*r*sin(θ)*(1.25 + r*cos(θ))/(1.5625 + 2.5*r*cos(θ) + r^2)^(3/2)))/(r^2*sin(θ))

    Ẋ[2] = (1.0*(1 + (-1/2)*(2*r + 1.6*cos(θ))/sqrt(0.64 + 1.6*r*cos(θ) + r^2)) - 0.8*(-cos(θ)/sqrt(0.64 + 1.6*r*cos(θ) + r^2) + (1/2)*(0.8 + r*cos(θ))*(2*r + 1.6*cos(θ))/(0.64 + 1.6*r*cos(θ) + r^2)^(3/2)) - 1.0*(-cos(θ)/sqrt(1.5625 + 2.5*r*cos(θ) + r^2) + (1/2)*(1.25 + r*cos(θ))*(2*r + 2.5*cos(θ))/(1.5625 + 2.5*r*cos(θ) + r^2)^(3/2)))/(r*sin(θ))
end

function ODErhs_empty!(Ẋ,X,p,t)
# This function takes as input a 2-vector (r,θ) andreturns a two-vector, (rdot, θdot).
# The symbolic functions for these two speeds have been calculated from Mathematica.

    r=X[1]
    θ=X[2]
    # Equations of motion
    Ẋ[1] = 0
    Ẋ[2] = 0
end

###################################
### Energy manifold definitions ###
###################################

function energyManifold_right(resid,u,p,t)
	global initval1,ψvalue
    r = u[1]
    θ = u[2]
    resid[1] = ψvalue - ψ_right(r,θ)
    resid[2] = 0
end

function energyManifold_left(resid,u,p,t)
	global initval2,ψvalue
    r = u[1]
    θ = u[2]
    resid[1] = ψvalue - ψ_left(r,θ)
    resid[2] = 0
#println("--Left energy manifold called, residue = ",resid[1])
end

#############################
### Source-sink functions ###
#############################

function sinkRight(u,t,integrator)
    # This function encodes a condition which is true when we want the event to occur,
    # and false otherwise.
    r = u[1]
    θ = u[2]
    global ϵ
    if ((r *cos(θ) < 4/5 && abs(r *sin(θ)) < ϵ && r *cos(θ) > 0 ) || r < ϵ && ( θ <= -π/2 || θ >= π/2 ))
        return true
    else
        return false
    end
end

function ejectRight!(integrator)
global stepnum,righttrapped,ψvalue,ϵ
    if righttrapped
        # Forget about the most recent few time steps and start over from
        # the beginning of this integration cycle.
        println("sink-eject detected but particle was already trapped")
        r1 = integrator.sol[1,1]
        θ1 = integrator.sol[2,1]
        t1 = integrator.sol.t[1]
        L = 4/5
        tunnelspeed = 10 # can tweak this
        tunneldist = L-abs(r1*cos(θ1))
        tunneltime = tunneldist/tunnelspeed
        ejecttime = t1 + tunneltime
        # The precise location to eject them cannot be calculated starting
        # from this intermediate location. It must be found from a global
        # variable stored previously.
        sgn1 = sign(sin(θ1))
        ejectLoc = nlsolve(findEjectLocationRight!, [4/5,sgn1*ϵ,sgn1*π/4], method =:newton)
        r2 = ejectLoc.zero[1]
        θ2 = ejectLoc.zero[2]
        ψ2 = ψ_right(r2,θ2)
        integrator.u[1] = r2
        integrator.u[2] = θ2
        integrator.t = t1 + tunneltime
        righttrapped = false
        # Report
        println("... Report of sink-eject event ...")
        println("Before ejection, (r,θ)        = (",r1,",",θ1,")")
        println("Entered tube at time = ",t1)
        println("Traveled inside tube = ",tunneltime)
        println("Ejected from tube at = ",ejecttime)
        println("Value of ψ before sink-eject = ",ψvalue)
        println("Value of ψ after sink-eject  = ",ψ2)
        println("After ejection, (r,θ)        = (",r2,",",θ2,")")
    else

        # This function releases particle close to r = L+ϵ, θ = ±(π-0.1)
        # Current value of r and θ
        r1 = integrator.u[1]
        θ1 = integrator.u[2]
        x1 = r1*cos(θ1)
        y1 = r1*sin(θ1)
        t1 = integrator.t

		# declared global so that the function findEjectLocationLeft! can find it.
        ψvalue = ψ_right(r1,θ1)
        sgn1 = sign(sin(θ1))

        # Solve a nonlinear system of equations for the value of
        # r and θ which is close to the point source, and has the
        # same value of the streamfunction as we currently have.
        ejectLoc = nlsolve(findEjectLocationRight!, [4/5,sgn1*ϵ,sgn1*π/4], method =:newton)

        # The values of r and θ at which to emerge are:
        r2 = ejectLoc.zero[1]
        θ2 = ejectLoc.zero[2]
        ψ2 = ψ_right(r2,θ2)

        # Tunnel speed and time
        L = 4/5
        tunnelspeed = 10 # can tweak this
        tunneldist = L-abs(r1*cos(θ1))
        tunneltime = tunneldist/tunnelspeed
        ejecttime = t1 + tunneltime

        # Assign the values to the integrator.

        if ejecttime > stepnum*switchΔt
            # i.e. the particle will be trapped.
            # determine how far it will travel.
            Δt = stepnum*switchΔt - t1
            x2 = x1 + Δt*tunnelspeed
            y2 = y1
            r_trapped = sqrt(x2^2 + y2^2)
            θ_trapped = acos(x2/r_trapped)
            righttrapped = true
            integrator.u[1] = r_trapped
            integrator.u[2] = θ_trapped
            integrator.t = stepnum*switchΔt
            # Report
            println("... Report of sink-eject event ....")
            println("Before ejection, (x,y)        = (",x1,",",y1,")")
            println("Value of ψ before sink-eject = ",ψvalue)
            println("Value of ψ after sink-eject  = ",ψ2)
            println("After ejection, (r,θ)  should be (",r2,",",θ2,")")
            println("Entered tube at time = ",t1)
            println("Targeted ejection time was = ",ejecttime)
            println("Proposed tunnel duration exceeds total available time!")
            println("Will be ejected instead at ",ejecttime+switchΔt)
            println("Traveled only for Δt = ",Δt)
            println("Particle stopped inside pipe")
            println("Stopped at (x,y)= (",x2,",",y2,")")
            println("Stopped at (r,θ)= (",r_trapped,",",θ_trapped,")")
            println("Sending integrator to ",stepnum*switchΔt)
        else
            integrator.u[1] = r2
            integrator.u[2] = θ2
            integrator.t += tunneltime
            # Report
            println("... Report of sink-eject event ...")
            println("Before ejection, (r,θ)        = (",r1,",",θ1,")")
            println("Entered tube at time = ",t1)
            println("Traveled inside tube = ",tunneltime)
            println("Ejected from tube at = ",ejecttime)
            println("Value of ψ before sink-eject = ",ψvalue)
            println("Value of ψ after sink-eject  = ",ψ2)
            println("After ejection, (r,θ)        = (",r2,",",θ2,")")
            if ~isapprox(ψ2,ψvalue)
                println("----Warning---Too much deviation from required ψ---")
            end
        end
        println("..................................")
    end
end

function sinkLeft(u,t,integrator)
    # This function encodes a condition which is true when we want the event to occur,
    # and false otherwise.
    r = u[1]
    θ = u[2]
    global ϵ
    if ((r *cos(θ) > -4/5 && abs(r *sin(θ)) < ϵ && r *cos(θ) < 0 ) || r < ϵ && ( θ >= -π/2 || θ <= π/2 ))
        return true
    else
        return false
    end
end

function ejectLeft!(integrator)
    global stepnum,lefttrapped,ψvalue,ϵ
    if lefttrapped
        # Forget about the most recent few time steps and start over from
        # the beginning of this integration cycle.
        println("sink-eject detected but particle was already trapped")
        r1 = integrator.sol[1,1]
        θ1 = integrator.sol[2,1]
        t1 = integrator.sol.t[1]
        L = 4/5
        tunnelspeed = 10 # can tweak this
        tunneldist = L-abs(r1*cos(θ1))
        tunneltime = tunneldist/tunnelspeed
        ejecttime = t1 + tunneltime
        # The precise location to eject them cannot be calculated starting
        # from this intermediate location. It must be found from a global
        # variable stored previously.
        sgn1 = sign(sin(θ1))
        ejectLoc = nlsolve(findEjectLocationLeft!, [4/5,(π-ϵ*sgn1),π-sgn1*π/4], method =:newton)
        r2 = ejectLoc.zero[1]
        θ2 = ejectLoc.zero[2]
        ψ2 = ψ_left(r2,θ2)
        integrator.u[1] = r2
        integrator.u[2] = θ2
        integrator.t = t1 + tunneltime
        lefttrapped = false
        # Report
        println("... Report of sink-eject event ...")
        println("Before ejection, (r,θ)        = (",r1,",",θ1,")")
        println("Entered tube at time = ",t1)
        println("Traveled inside tube = ",tunneltime)
        println("Ejected from tube at = ",ejecttime)
        println("Value of ψ before sink-eject = ",ψvalue)
        println("Value of ψ after sink-eject  = ",ψ2)
        println("After ejection, (r,θ)        = (",r2,",",θ2,")")
    else

        # This function releases particle close to r = L+ϵ, θ = ±(π-0.1)
        # Current value of r and θ
        r1 = integrator.u[1]
        θ1 = integrator.u[2]
        x1 = r1*cos(θ1)
        y1 = r1*sin(θ1)
        t1 = integrator.t

        # declared global so that the function findEjectLocationLeft! can find it.
        global ψvalue = ψ_left(r1,θ1)
        sgn1 = sign(sin(θ1))

        # Solve a nonlinear system of equations for the value of
        # r and θ which is close to the point source, and has the
        # same value of the streamfunction as we currently have.
        ejectLoc = nlsolve(findEjectLocationLeft!, [4/5,(π-ϵ*sgn1),π-sgn1*π/4], method =:newton)

        # The values of r and θ at which to emerge are:
        r2 = ejectLoc.zero[1]
        θ2 = ejectLoc.zero[2]
        ψ2 = ψ_left(r2,θ2)

        # Tunnel speed and time
        L = 4/5
        tunnelspeed = 10 # can tweak this
        tunneldist = L-abs(r1*cos(θ1))
        tunneltime = tunneldist/tunnelspeed
        ejecttime = t1 + tunneltime

        # Assign the values to the integrator.

        if ejecttime > stepnum*switchΔt
            # i.e. the particle will be trapped.
            # determine how far it will travel.
            Δt = stepnum*switchΔt - t1
            x2 = x1 - Δt*tunnelspeed
            y2 = y1
            r_trapped = sqrt(x2^2 + y2^2)
            θ_trapped = acos(x2/r_trapped)
            lefttrapped = true
            integrator.u[1] = r_trapped
            integrator.u[2] = θ_trapped
            integrator.t = stepnum*switchΔt
            # Report
            println("... Report of sink-eject event ....")
            println("Before ejection, (x,y)        = (",x1,",",y1,")")
            println("Value of ψ before sink-eject = ",ψvalue)
            println("Value of ψ after sink-eject  = ",ψ2)
            println("After ejection, (r,θ)  should be (",r2,",",θ2,")")
            println("Entered tube at time = ",t1)
            println("Targeted ejection time was = ",ejecttime)
            println("Proposed tunnel duration exceeds total available time!")
            println("Will be ejected instead at ",ejecttime+switchΔt)
            println("Traveled only for Δt = ",Δt)
            println("Particle stopped inside pipe")
            println("Stopped at (x,y)= (",x2,",",y2,")")
            println("Stopped at (r,θ)= (",r_trapped,",",θ_trapped,")")
            println("Sending integrator to ",stepnum*switchΔt)
        else
            integrator.u[1] = r2
            integrator.u[2] = θ2
            integrator.t += tunneltime
            # Report
            println("... Report of sink-eject event ...")
            println("Before ejection, (r,θ)        = (",r1,",",θ1,")")
            println("Entered tube at time = ",t1)
            println("Traveled inside tube = ",tunneltime)
            println("Ejected from tube at = ",ejecttime)
            println("Value of ψ before sink-eject = ",ψvalue)
            println("Value of ψ after sink-eject  = ",ψ2)
            println("After ejection, (r,θ)        = (",r2,",",θ2,")")
            if ~isapprox(ψ2,ψvalue)
                println("----Warning---Too much deviation from required ψ---")
            end
        end
        println("..................................")
    end
end

################################
### Calculate eject location ###
################################
function findEjectLocationRight!(F,x)
    # This function is defined such that when it is equal to (0,0,0),
    # then we have found r, θ and ϕ corresponding to the location on
    # the ejection semi-circle where we would like to start off our
    # particle trajectories.

    global ψvalue # calculated before this function was called.
    global ϵ
    L = 4/5

    # The following three values of r, θ and ϕ
    # are the "guess" values passed to this function
    # from which the guess will start. They should be in
    # a "reasonable neighborhood" of the point source.
    # These are the values which will be iterated upon.

    r = x[1]
    θ = x[2]
    ϕ = x[3]

    # Define the function F(X) which should go to zero at the
    # required value of X.
    F[1] = ψ_right(r,θ) - ψvalue
    F[2] = L^2+ϵ^2+2*L*ϵ*cos(ϕ)-r^2
    F[3] = r - ϵ *(sin(ϕ))/(sin(θ))
end

function findEjectLocationLeft!(F,x)
    # This function is defined such that when it is equal to (0,0,0),
    # then we have found r, θ and ϕ corresponding to the location on
    # the ejection semi-circle where we would like to start off our
    # particle trajectories.

    global ψvalue	# calculated before this function was called.
    global ϵ
    L = 4/5

    # The following three values of r, θ and ϕ
    # are the "guess" values passed to this function
    # from which the guess will start. They should be in
    # a "reasonable neighborhood" of the point source.
    # These are the values which will be iterated upon.

    r = x[1]
    θ = x[2]
    ϕ = x[3]

    # Define the function F(X) which should go to zero at the
    # required value of X.
    F[1] = ψ_left(r,θ) - ψvalue
    F[2] = L^2+ϵ^2-2*L*ϵ*cos(ϕ)-r^2
    F[3] = r - ϵ *(sin(ϕ))/(sin(θ))

end

######################
### Streamfunction ###
######################
function ψ_right(r,θ)
    H = -0.8*(1 - (0.8 - r*cos(θ))/sqrt(0.64 - 1.6*r*cos(θ) + r^2)) - 1.0*(1 - (1.25 - r*cos(θ))/sqrt(1.5625 - 2.5*r*cos(θ) + r^2)) + 1.0*(0.8 +  r - sqrt(0.64 - 1.6*r*cos(θ) + r^2))
end

function ψ_left(r,θ)
    H = -(-0.8*(1 - (0.8 + r*cos(θ))/sqrt(0.64 + 1.6*r*cos(θ) + r^2)) - 1.0*(1 - (1.25 + r*cos(θ))/sqrt(1.5625 + 2.5*r*cos(θ) + r^2)) + 1.0*(0.8 + r - sqrt(0.64 + 1.6*r*cos(θ) + r^2)))
end

############################################
### Plotting and miscellaneous functions ###
############################################

function showTrajectory(r,θ,H,t)
    p1 = plot(θ, r, proj = :polar, color=:black, seriestype=:scatter,
        markersize=0.2,ylim=(0,1),legend=:false,size=(700,700))
    p2 = plot(t,H,ylims=(-1,1),title="Streamfunction",
        seriestype=:scatter,markersize=0.2,linestyle=:solid,legend=:none,size=(500,250))
    plot(p1,p2,dpi=600)
end

function short(variable)
    if typeof(variable)==Array{Float64,1}
        return round.(variable,digits=3)
    else
        return round(variable,digits=3)
    end
end

