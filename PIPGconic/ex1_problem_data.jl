# Example 1 : 	Oscillating masses

# Constraints :
# 				Position in box
# 				Velocity in box
#				Acceleration in box
#       		Acceleration rate in box
		
# problem definition
module prb

	using LinearAlgebra, StaticArrays

	const nx = 8											# state dim
	const nu = 4											# imput dim
	const N = 30											# horizon length
	const Δ = 0.25 											# sampling time
	const tvec = [(j-1)*Δ for j in 1:N]   					# time vector (t=0:N-1)

	# solution method settings
		# common
		const kmax = 500										# max iterations for main algorithms (and outer loop)
		const scale_Hg = true
		
		# JuMP
		const ϵ_pd = 1e-9										# solver precision
		const ϵ_gap = 1e-9										# solver precision
	
		# PIPG
		const ϵ_pipg = 1e-4										# termination precision 
		const restrt_idx = 40									# index for restarting	

		# ADMM
		const α_admm = 1.0										# step size
		const ϵ_nest_admm = 1e-10								# Nesterov's method termination precision 
		const kmax_nest = 100									# max iterations for Nesterov's method

		# Chambolle and Pock (PDHG)
		const β⁰ = 2.0											# initial step size		
		const γ⁰ = 1.0											# initial step size

	# system parameters
	const nxb2 = Int64(nx/2)									# position veloity dim
	const pmax = 2.0											# position	
	const vmax = 2.0											# velocity		
	const umax = 2.0					 						# input (forcing)	
	const dumax = 0.5											# input rate

	# reference trajectory	
	const xref = [SVector{nx}(cat(ones(nxb2)...,zeros(nxb2)...,dims=1)) for _ in 1:N-1]	# state reference (t=1:N-1)
	const uref = [SVector{nu}(zeros(nu)) for _ in 1:N-1]							# input reference (t=0:N-2)
	const x0 = 1 .* ones(nx)								    					# initial condition (t=0)

	# dynamics (t=0:N-2)
	const Ac = Array([zeros(nxb2,nxb2) Diagonal(ones(nxb2))
						       Tridiagonal(ones(nxb2-1),-2 .* ones(nxb2),ones(nxb2-1)) zeros(nxb2,nxb2)])
	const Ad = SMatrix{nx,nx}(exp(Ac*Δ))
	const Bc = [zeros(nxb2,nxb2) 
				Diagonal(ones(nxb2,nxb2))]
	const Bd = SMatrix{nx,nu}(inv(Ac)*(Ad-Diagonal(ones(nx)))*Bc)

	const A = [Ad for _ in 1:N-1]
	const B = [Bd for _ in 1:N-1]
	const e = [SVector{nx}(zeros(nx)) for _ in 1:N-1]
	
	# state and control penalty (LQR weights)

	const q_tmp = [1.0,1.0]
	const r_tmp = 1.0

	# const q_tmp = [0.5,1.0]
	# const r_tmp = 1.0

	const Q = [SMatrix{nx,nx}(Array(Diagonal(cat(q_tmp[1] .* ones(nxb2), q_tmp[2] .* ones(nxb2), dims=1)))) for _ in 1:N-1] # t=1:N-1
	const R = [SMatrix{nu,nu}(r_tmp .* Array(Diagonal(ones(nu)))) for _ in 1:N-1] 											# t=0:N-2

end

include("proj_funcs.jl") 		# projection library

# assemble quantities for vectorized PIPG+ implementation
module asm
	
	using ..prb
	using ..proj 								# projection onto simple sets
	using LinearAlgebra, StaticArrays, JuMP

	const scale_Hg = prb.scale_Hg
	const nH = 2*prb.nu*(prb.N-2)				# no. of inequality constraints

	function construct_H(nx,nu,N,A,B,Inx,Inu)
		nxnu = nx+nu
		H = zeros(nx*(N-1) + nH, nxnu*(N-1))
		H[1:nx,1:nxnu] .= [-B[1] Inx]
		for j=1:N-2
			H[j*nx+1:(j+1)*nx, nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+2*nx+nu] .= [-A[j] -B[j] Inx]

			H[nx*(N-1)+(j-1)*nu+1:nx*(N-1)+j*nu, (j-1)*nxnu+1:(j-1)*nxnu+nu+nxnu] .= [Inu zeros(nu,nx) -Inu]

			H[nx*(N-1)+nu*(N-2)+(j-1)*nu+1:nx*(N-1)+nu*(N-2)+j*nu, (j-1)*nxnu+1:(j-1)*nxnu+nu+nxnu] .= [-Inu zeros(nu,nx) Inu]
		end

		norm_row_H = zeros(nx*(N-1)+nH)
		for j=1:nx*(N-1)+nH
			norm_row_H[j] = norm(H[j,1:end],Inf)
		end

		if scale_Hg
			# print("Scaling rows of H\n")
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
		# print("Scaling components of g\n")
		const scl_val_Hg = norm_row_H	
	else
		const scl_val_Hg = 1.0
	end

	const g = cat(prb.A[1]*prb.x0+prb.e[1],cat([prb.e[t] for t=2:prb.N-1]...,dims=1), -prb.dumax .* ones(2*prb.nu*(prb.N-2)),dims=1)

	@assert length(h) == (prb.nx+prb.nu)*(prb.N-1) "Length of h is incorrect."
	@assert length(g) == 2*prb.nu*(prb.N-2)+prb.nx*(prb.N-1) "Length of g is incorrect."

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
		return [z[nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+nx] for j=1:N-1], [z[(j-1)*nxnu+1:(j-1)*nxnu+nu] for j=1:N-1] 
	end

	# useful constants
	const nx = copy(prb.nx)
	const nu = copy(prb.nu)
	const nxnu = nx+nu
	const nxb2 = Int64(nx/2)
	const N = copy(prb.N)

	const κ_x = MVector{nx}(randn(nx))
	const κ_u = MVector{nu}(randn(nu))

	# projection
	function project_U!(u)
		κ_u .= u
		proj.box!(u,κ_u,-prb.umax,prb.umax,nu)
	end
	function project_X!(x)
		κ_x .= x
		proj.box!(x[1:nxb2],κ_x[1:nxb2],-prb.pmax,prb.pmax,nxb2)
		proj.box!(x[nxb2+1:nx],κ_x[nxb2+1:nx],-prb.vmax,prb.vmax,nxb2)
	end

	# set state and input constraints in JuMP (vectorized version)
	function set_constraint_jump_vectorized!(model,zz)
		# acceleration, position and velocity in box
		@constraint(model,[j=1:N-1],-prb.umax .≤ zz[(j-1)*nxnu+1:(j-1)*nxnu+nu] .≤  prb.umax)
		@constraint(model,[j=1:N-1],-prb.pmax .≤ zz[nu+(j-1)*nxnu+1:nu+(j-1)*nxnu+nxb2] .≤  prb.pmax)
		@constraint(model,[j=1:N-1],-prb.vmax .≤ zz[nu+(j-1)*nxnu+nxb2+1:nu+(j-1)*nxnu+nx] .≤  prb.vmax)
	end

	# set state and input constraints (everything other than dynamics and initial condition) in JuMP (decomposed version)
	function set_constraint_jump_decomposed!(model,x,u)
		@constraint(model,[t=1:N-1], -prb.pmax .≤ x[1:nxb2,t+1] .≤ prb.pmax)
		@constraint(model,[t=1:N-1], -prb.vmax .≤ x[nxb2+1:nx,t+1] .≤ prb.vmax)
		@constraint(model,[t=1:N-1], -prb.umax .≤ u[1:nu,t] .≤ prb.umax)
		@constraint(model,[t=1:N-2],-prb.dumax .≤ u[:,t+1] .- u[:,t] .≤ prb.dumax)
	end

end