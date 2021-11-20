# Example 2   : Quadrotor path planning

# Constraints :
# 				Position in box
# 				Velocity in ball
# 				Thrust in ball ∩ cone
#				Thrust in halfspace (z component lower bound) 
# 				Thrust rate in box
		
# problem definition
module prb

	using LinearAlgebra, StaticArrays

	const nx = 6												# state dim
	const nu = 3												# imput dim
	const N = 30												# horizon length
	const Δ = 0.25 												# sampling time
	const tvec = [(j-1)*Δ for j in 1:N]   						# time vector (t=0:N-1)

	# solution method settings
		# common
		const kmax = 2000										# max iterations for main algorithms (and outer loop)
		const scale_Hg = true
		
		# JuMP
		const ϵ_pd = 1e-9										# solver precision
		const ϵ_gap = 1e-9										# solver precision
	
		# PIPG
		const ϵ_pipg = 1e-4										# termination precision 
		const restrt_idx = 300									# index for restarting	

		# ADMM
		const α_admm = 1.0										# step size
		const ϵ_nest_admm = 1e-10								# Nesterov's method termination precision 
		const kmax_nest = 100									# max iterations for Nesterov's method

		# Chambolle and Pock
		const β⁰ = 2.0											# initial step size		
		const γ⁰ = 1.0											# initial step size

	# system parameters
	const nxb2 = Int64(nx/2)									# position veloity dim
	const α_tl = tand( 45 )										# tangent of max thrust pointing angle

	const m = 0.35												# quadrotor mass
	const gₑ = 9.806	 										# acceleration due to gravity at Earth sea-level (m s⁻²)
	const pmax = 3.0											# max position ∞ norm	
	const vmax = 5.0											# max velocity 2 norm 	 	
	const umax = 5.0					 						# max input (forcing) 2 norm 	
	const umin = 2.0											# min z component for input
	const dumax = 3.0										 	# max input rate ∞ norm

	# const p1 = SVector{3}([-10.0,-10.0,0.0])  				# center of obstacle #1 
	# const p2 = SVector{3}([8.0,-8.0,0.0])	  					# center of obstacle #2	
	# const p3 = SVector{3}([10.0,10.0,0.0])	  				# center of obstacle #3

	# const r1 = 7.6										  	# radius of obstacle #1
	# const r2 = 9.0						      			    # radius of obstacle #2
	# const r3 = 7.6                            			    # radius of obstacle #3 

	const p1 = SVector{3}([-1.5,-1.5,0.0])  					# center of obstacle #1 
	const p2 = SVector{3}([1.2,-1.2,0.0])	  					# center of obstacle #2	
	const p3 = SVector{3}([1.5,1.5,0.0])	  					# center of obstacle #3

	const r1 = 0.8										  	  	# radius of obstacle #1
	const r2 = 1.2						      					# radius of obstacle #2
	const r3 = 0.8                            					# radius of obstacle #3 

	# reference trajectory	
	const x0 = SVector{nx}(cat(p1 .- 
		[0.0,1.2*r1,0.0],[0.0,0.0,0.0],dims=1))					# initial condition (t=0)
	const xf = SVector{nx}(cat(p3 .+ 
		[1.2*r3,0.0,0.0],[0.0,0.0,0.0],dims=1))
	const xref = [SVector{nx}((1-t/(N-1)) .* x0 .+ 
		(t/(N-1)) .* xf) for t in 1:N-1]		        		# state reference (t=1:N-1)
	const uref = [SVector{nu}(zeros(nu)) for _ in 1:N-1]	    # input reference (t=0:N-2)

	# dynamics (t=0:N-2)
	const A = [SMatrix{nx,nx}([1.0 0.0 0.0 Δ   0.0 0.0
							   0.0 1.0 0.0 0.0 Δ   0.0
							   0.0 0.0 1.0 0.0 0.0 Δ  
							   0.0 0.0 0.0 1.0 0.0 0.0
							   0.0 0.0 0.0 0.0 1.0 0.0
							   0.0 0.0 0.0 0.0 0.0 1.0]) for _ in 1:N-1]
	const B = [SMatrix{nx,nu}([0.5*(Δ^2)/m  0.0          0.0
							   0.0          0.5*(Δ^2)/m  0.0
							   0.0          0.0          0.5*(Δ^2)/m
							   Δ/m          0.0          0.0
							   0.0          Δ/m          0.0
							   0.0          0.0          Δ/m]) for _ in 1:N-1]
	const e = [SVector{nx}([0.0, 0.0, -0.5*(Δ^2)*gₑ, 0.0, 0.0, -Δ*gₑ]) for _ in 1:N-1]

	# state and control penalty (LQR weights)

	# const q_tmp = [1.0,3.0]
	# const r_tmp = 5.0

	const q_tmp = [1.0,2.5]
	const r_tmp = 0.5

	const Q = [SMatrix{nx,nx}(Array(Diagonal(cat(q_tmp[1] .* ones(nxb2), q_tmp[2] .* ones(nxb2), dims=1)))) for _ in 1:N-1] # t=1:N-1
	const R = [SMatrix{nu,nu}(r_tmp .* Array(Diagonal(ones(nu)))) for _ in 1:N-1] 											# t=0:N-2

end

# auxiliary definitions for the the problem
module prb2
	using ..prb
	using LinearAlgebra, StaticArrays

	const nx = copy(prb.nx)
	const nu = copy(prb.nu)
	const nxb2 = copy(prb.nxb2)
	const N = copy(prb.N)

	# obstacle avoidance
	# linearization
	# initial and final point for obstacle linearization reference
	const x̄0 = SVector{nx}(cat(prb.p1 .- [0.0,prb.r1,0.0],[0.0,0.0,0.0],dims=1))
	const x̄f = SVector{nx}(cat(prb.p3 .+ [prb.r3,0.0,0.0],[0.0,0.0,0.0],dims=1))
	# obstructed straight-line reference for obstacle linearization
	const x̄_obs = [SVector{nx}( (1-t/(N-1)) .* x̄0 .+ (t/(N-1)) .* x̄f ) for t=1:(N-1)]

	# move state to nearest obstacle boundary if penetrated
	function path_unobstructed(z)
		# globals used p1, p2, p3
		yy = zeros(prb.nx)
		yy .= z

		tol_fac = 1.0

		nrm_y_p1 = norm(yy[1:3]-prb.p1)
		if nrm_y_p1 ≤ tol_fac*prb.r1
			yy[1:3] .= (yy[1:3] .- prb.p1) .* (tol_fac*prb.r1/nrm_y_p1) .+ prb.p1 	
		end

		nrm_y_p2 = norm(yy[1:3] .- prb.p2)
		if nrm_y_p2 ≤ tol_fac*prb.r2
			yy[1:3] .= (yy[1:3] .- prb.p2) .* (tol_fac*prb.r2/nrm_y_p2) .+ prb.p2 	
		end

		nrm_y_p3 = norm(yy[1:3]-prb.p3)
		if nrm_y_p3 ≤ tol_fac*prb.r3
			yy[1:3] .= (yy[1:3] .- prb.p3) .* (tol_fac*prb.r3/nrm_y_p3) .+ prb.p3 	
		end
		return yy
	end

	# unobstructed straight-line reference for obstacle linearization
	const x̄_unobs = [SVector{nx}(path_unobstructed(x̄_obs[t])) for t=1:(N-1)]

	# terms in the obstacle avoidance constraint linearization
	#  - ϕᵢ⋅x + θᵢ ≈ norm(x[1:3]-pᵢ) - rᵢ ≥ 0
	const θ₁ = SVector{N-1}([norm(x̄_unobs[t][1:3] .- prb.p1) - prb.r1 + (-1 * norm(x̄_unobs[t][1:3])^2 + prb.p1⋅x̄_unobs[t][1:3]) / norm(x̄_unobs[t][1:3] .- prb.p1) for t=1:(N-1)])
	const θ₂ = SVector{N-1}([norm(x̄_unobs[t][1:3] .- prb.p2) - prb.r2 + (-1 * norm(x̄_unobs[t][1:3])^2 + prb.p2⋅x̄_unobs[t][1:3]) / norm(x̄_unobs[t][1:3] .- prb.p2) for t=1:(N-1)])
	const θ₃ = SVector{N-1}([norm(x̄_unobs[t][1:3] .- prb.p3) - prb.r3 + (-1 * norm(x̄_unobs[t][1:3])^2 + prb.p3⋅x̄_unobs[t][1:3]) / norm(x̄_unobs[t][1:3] .- prb.p3) for t=1:(N-1)])			

	const ϕ₁ = [SVector{nxb2}(-1 .* (x̄_unobs[t][1:3] .- prb.p1) ./ norm(x̄_unobs[t][1:3] .- prb.p1) ) for t=1:N-1]
	const ϕ₂ = [SVector{nxb2}(-1 .* (x̄_unobs[t][1:3] .- prb.p2) ./ norm(x̄_unobs[t][1:3] .- prb.p2) ) for t=1:N-1]
	const ϕ₃ = [SVector{nxb2}(-1 .* (x̄_unobs[t][1:3] .- prb.p3) ./ norm(x̄_unobs[t][1:3] .- prb.p3) ) for t=1:N-1]
end


include("proj_funcs.jl") 					# projection library

# assemble quantities for vectorized PIPG+ implementation
module asm
	
	using ..prb
	using ..prb2
	using ..proj 											# projection onto simple sets
	using LinearAlgebra, StaticArrays, JuMP

	const scale_Hg = prb.scale_Hg
	const nH = 2*prb.nu*(prb.N-2)+(prb.N-1)+3*(prb.N-1)	    # no. of inequality constraints

	function construct_H(nx,nu,N,A,B,Inx,Inu)
		nxnu = nx+nu
		nxb2 = Int64(nx/2)
		H = zeros(nx*(N-1) + nH, nxnu*(N-1))
		H[1:nx,1:nxnu] .= [-B[1] Inx]
		for j=1:N-2
			H[j*nx+1:(j+1)*nx, nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+2*nx+nu] .= [-A[j] -B[j] Inx]

			H[nx*(N-1)+(j-1)*nu+1:nx*(N-1)+j*nu, (j-1)*nxnu+1:(j-1)*nxnu+nu+nxnu] .= [Inu zeros(nu,nx) -Inu]

			H[nx*(N-1)+nu*(N-2)+(j-1)*nu+1:nx*(N-1)+nu*(N-2)+j*nu, (j-1)*nxnu+1:(j-1)*nxnu+nu+nxnu] .= [-Inu zeros(nu,nx) Inu]
		end
		for j=1:N-1
			H[nx*(N-1)+2*nu*(N-2)+j,(j-1)*nxnu+1:(j-1)*nxnu+nu] .= [0.0,0.0,1.0]

			H[nx*(N-1)+2*nu*(N-2)+(N-1)+j,nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+nxb2] .= -1 .* prb2.ϕ₁[j]

			H[nx*(N-1)+2*nu*(N-2)+2*(N-1)+j,nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+nxb2] .= -1 .* prb2.ϕ₂[j]

			H[nx*(N-1)+2*nu*(N-2)+3*(N-1)+j,nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+nxb2] .= -1 .* prb2.ϕ₃[j]
		end

		norm_row_H = zeros(nx*(N-1)+nH)
		for j=1:nx*(N-1)+nH
			norm_row_H[j] = norm(H[j,1:end],Inf)
		end

		if scale_Hg
			print("Scaling rows of H\n")
			for j=1:length(norm_row_H)
				H[j,1:end] .= H[j,1:end] ./ norm_row_H[j]
			end
		end

		return H,norm_row_H
	end

	function construct_G(nx,nu,N,Q,R)
		nxnu = nx+nu
		G = zeros(nxnu*(N-1),nxnu*(N-1))
		for j=1:N-1
			G[(j-1)*nxnu+1:j*nxnu,(j-1)*nxnu+1:j*nxnu] .= [R[j] 		zeros(nu,nx)
														   zeros(nx,nu) Q[j]]
		end
		return G
	end

	H_scl,norm_row_H = construct_H(prb.nx,prb.nu,prb.N,prb.A,prb.B,Diagonal(ones(prb.nx)),Diagonal(ones(prb.nu)))
	const H = deepcopy(H_scl)
	const HT = transpose(H)

	const zref = cat([cat(prb.uref[t],prb.xref[t],dims=1) for t=1:prb.N-1]...,dims=1)
	const h = cat([cat(-prb.R[t]*prb.uref[t], -prb.Q[t]*prb.xref[t], dims=1) for t=1:prb.N-1]...,dims=1)
	
	if scale_Hg
		print("Scaling components of g\n")
		const scl_val_Hg = norm_row_H	
	else
		const scl_val_Hg = 1.0
	end

	const g = cat(prb.A[1]*prb.x0 .+ prb.e[1],cat([prb.e[t] for t=2:prb.N-1]...,dims=1), -prb.dumax .* ones(2*prb.nu*(prb.N-2)), 
			prb.umin .* ones(prb.N-1),-1 .* prb2.θ₁...,-1 .* prb2.θ₂...,-1 .* prb2.θ₃...,dims=1) ./ scl_val_Hg

	@assert length(h) == (prb.nx+prb.nu)*(prb.N-1) "Length of h is incorrect."
	@assert length(g) == 2*prb.nu*(prb.N-2)+(prb.N-1)+3*(prb.N-1)+prb.nx*(prb.N-1) "Length of g is incorrect."

	const Heq = [H -1 .* Array(Diagonal(ones(prb.nx*(prb.N-1) + nH)))]
	const HeqT = transpose(Heq)

	const G = construct_G(prb.nx,prb.nu,prb.N,prb.Q,prb.R)

	const svd_HTH = svd(HT*H)
	const svd_HeqTHeq = svd(HeqT*Heq)
	
	const svd_G = svd(G)

	const σ = max(svd_HTH.S...)
	const σ_eq = max(svd_HeqTHeq.S...)
	const λ = max(svd_G.S...)
	const μ = min(svd_G.S...)

	function construct_xu(z::Array{Float64,1})
		return cat([prb.x0],[z[nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+nx] for j=1:N-1],dims=1), [z[(j-1)*nxnu+1:(j-1)*nxnu+nu] for j=1:N-1] 
	end

	# useful constants
	const nx = copy(prb.nx)
	const nu = copy(prb.nu)
	const nxnu = nx+nu
	const nxb2 = Int64(nx/2)
	const N = copy(prb.N)

	const κ_x = MVector{nx}(randn(nx))
	const κ_xb2 = MVector{nxb2}(randn(nxb2))
	const κ_u = MVector{nu}(randn(nu))

	# projection
	function project_U!(u)
		κ_u .= u
		γ_u = randn()
		proj.soc_ball!(u,κ_u,prb.α_tl,prb.umax,nu,γ_u)
	end
	function project_X!(x)
		κ_x .= x
		proj.box!(κ_xb2,κ_x[1:nxb2],-prb.pmax,prb.pmax,nxb2)
		x[1:nxb2] .= κ_xb2 
		proj.ball!(κ_xb2,κ_x[nxb2+1:nx],prb.vmax,nxb2)
		x[nxb2+1:nx] .= κ_xb2
	end

	# set state and input constraints in JuMP
	function set_constraint_jump_vectorized!(model,zz)
		# acceleration, position and velocity
		@constraint(model,[j=1:N-1],[prb.umax,zz[(j-1)*nxnu+1:(j-1)*nxnu+nu]...] in SecondOrderCone())
		@constraint(model,[j=1:N-1],[prb.α_tl*zz[(j-1)*nxnu+nu],zz[(j-1)*nxnu+1:(j-1)*nxnu+nu-1]...] in SecondOrderCone())
		@constraint(model,[j=1:N-1],-prb.pmax .≤ zz[nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+nxb2] .≤  prb.pmax)
		@constraint(model,[j=1:N-1],[prb.vmax,zz[nu+(j-1)*nxnu+nxb2+1:nu+(j-1)*nxnu+nx]...] in SecondOrderCone())
	end

	# set state and input constraints (everything other than dynamics and initial condition) in JuMP (decomposed version)
	function set_constraint_jump_decomposed!(model,x,u)

		@constraint(model,[t=1:N-1], -prb.pmax .≤ x[1:nxb2,t+1] .≤ prb.pmax)					
		@constraint(model,[t=1:N-1], [prb.vmax,x[nxb2+1:nx,t+1]...] in SecondOrderCone())
		@constraint(model,[t=1:N-1], [prb.umax,u[1:nu,t]...] in SecondOrderCone())
		@constraint(model,[t=1:N-1], [prb.α_tl*u[nu,t],u[1:nu-1,t]...] in SecondOrderCone())

		# linear inequality constraints
		@constraint(model,[t=1:N-2], -prb.dumax .≤ u[:,t+1] .- u[:,t] .≤ prb.dumax)
		@constraint(model,[t=1:N-1], u[3,t] ≥ prb.umin)
		@constraint(model,[t=1:N-1], prb2.ϕ₁[t]⋅x[1:nxb2,t+1] - prb2.θ₁[t] ≤ 0)
		@constraint(model,[t=1:N-1], prb2.ϕ₂[t]⋅x[1:nxb2,t+1] - prb2.θ₂[t] ≤ 0)
		@constraint(model,[t=1:N-1], prb2.ϕ₃[t]⋅x[1:nxb2,t+1] - prb2.θ₃[t] ≤ 0)
	end

end

module plotter
using LinearAlgebra
using Plots, LaTeXStrings
gr()
# This ensures that a legend entry is not created by default
default(lab="",markersize=2,markerstrokewidth=0.1,xtickfontsize=12,ytickfontsize=12,
	ztickfontsize=12,legendfontsize=9)
using ..prb
using ..prb2

function trajectory2D(x1,u1,x2,u2,slver=:ecos,pltflg=true)

	y_unscl = cat([prb.x0],prb.xref,dims=1)

	# Position

	# obstacle regions
	x = prb.p1[1] .+ collect(-prb.r1:(2*prb.r1/(100-1)):prb.r1)
	y = prb.p1[2] .+ sqrt.( abs.(prb.r1^2 .- (x .- prb.p1[1]) .^ 2) )
	p1 = plot(x,y,color=:black,line=:solid)
	y = prb.p1[2] .- sqrt.( abs.(prb.r1^2 .- (x .- prb.p1[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	x = prb.p2[1] .+ collect(-prb.r2:(2*prb.r2/(100-1)):prb.r2)
	y = prb.p2[2] .+ sqrt.( abs.(prb.r2^2 .- (x .- prb.p2[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)
	y = prb.p2[2] .- sqrt.( abs.(prb.r2^2 .- (x .- prb.p2[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	x = prb.p3[1] .+ collect(-prb.r3:(2*prb.r3/(100-1)):prb.r3)
	y = prb.p3[2] .+ sqrt.( abs.(prb.r3^2 .- (x .- prb.p3[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)
	y = prb.p3[2] .- sqrt.( abs.(prb.r3^2 .- (x .- prb.p3[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	plot!([prb.x0[1]],[prb.x0[2]],color=:orange,marker=:circle,linealpha=0.01,lab=L"x_0[1:3]")
	plot!([prb.xf[1]],[prb.xf[2]],color=:orange,marker=:square,linealpha=0.01,lab=L"x_f[1:3]")

	plot!([y_unscl[t][1] for t=1:prb.N],[y_unscl[t][2] for t=1:prb.N],color=:red,line=:solid,lw=0.7,lab="Cost Ref.")
	plot!([prb2.x̄_unobs[t][1] for t=1:prb.N-1],[prb2.x̄_unobs[t][2] for t=1:prb.N-1],color=:magenta,line=:solid,lw=1.3,lab="Obstcl. Ref.")
	# lim_set = [-max(prb.r1,prb.r2,prb.r3),max(prb.r1,prb.r2,prb.r3)]
	lim_set = [-1.0*prb.pmax,prb.pmax*1.0]
	plot!(xlim=lim_set,ylim=lim_set)

	if pltflg
		plot!([x1[t][1] for t=1:prb.N],[x1[t][2] for t=1:prb.N],color=:blue,line=:solid,lab="PIPG")
		plot!([x2[t][1] for t=1:prb.N],[x2[t][2] for t=1:prb.N],color=:cyan,line=:dash,lab=string("JuMP - ",String(slver)))
	end

	plot!(size=[500,500],legend=:topleft,title="Position")

	vmax = prb.vmax
	umax = prb.umax
	umin = prb.umin

	# Velocity
	p2 = plot(prb.tvec,vmax .* ones(prb.N),color=:black,line=:solid)

	if pltflg
		plot!(prb.tvec,[norm(x1[t][4:6]) for t=1:prb.N],color=:blue,line=:solid,lab="PIPG")
		plot!(prb.tvec,[norm(x2[t][4:6]) for t=1:prb.N],color=:cyan,line=:dash,lab=string("JuMP - ",String(slver)))
	end

	plot!(size=[500,500],legend=:none,title="Velocity")

	# Acceleration

	p3 = plot(prb.tvec[1:prb.N-1],umax .* ones(prb.N-1),color=:black,line=:solid,label=L"u_{max}")
	plot!(prb.tvec[1:prb.N-1],umin .* ones(prb.N-1),color=:black,line=:solid,label=L"u_{min}")

	if pltflg	
		plot!(prb.tvec[1:prb.N-1],[norm(u1[t]) for t=1:prb.N-1],color=:blue,line=:solid,lab="PIPG")		
		plot!(prb.tvec[1:prb.N-1],[norm(u2[t]) for t=1:prb.N-1],color=:cyan,line=:dash,lab=string("JuMP - ",String(slver)))
	end

	plot!(size=[500,500],legend=:topright,title="Acceleration",ylim=[0.8*umin,1.1*umax])

	λ = @layout [a b c]
	display(plot(p1,p2,p3,layout=λ,size=[999,333]))

	if pltflg
		p0 = plot(prb.tvec,[x1[t][3] for t=1:prb.N],color=:orange,line=:solid,lab="PIPG",title="Altitude",size=[400,400])
		plot!(prb.tvec,[x2[t][3] for t=1:prb.N],color=:green,line=:dash,lab=string("JuMP - ",String(slver)))
		display(p0)
	end
end

end