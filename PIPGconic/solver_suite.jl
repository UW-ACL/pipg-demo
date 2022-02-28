# projection functions for the vectorized problem
module prj
	using ..prb
	using ..asm
	using LinearAlgebra

	# useful constants
	const nx = copy(prb.nx)
	const nu = copy(prb.nu)
	const nxnu = nx+nu
	const N = copy(prb.N)

	function project_Kdual!(w,κK)
		κK .= w[nx*(N-1)+1:end]
		w[nx*(N-1)+1:end] .= max.(κK,0.0)
	end
	function project_Kdual(w,κK)
		ww = deepcopy(w)
		project_Kdual!(ww,κK)
		return ww
	end
	function project_K!(w,κK)
		w[1:nx*(N-1)] .= 0.0
		κK .= w[nx*(N-1)+1:end]
		w[nx*(N-1)+1:end] .= max.(κK,0.0)		 
	end
	function project_K(w,κK)
		ww = deepcopy(w)
		project_K!(ww,κK)
		return ww
	end
	function project_Kpolar!(w,κK)
		w .= w .- project_K(w,κK)
	end
	function project_Kpolar(w,κK)
		ww = deepcopy(w)
		project_Kpolar!(ww,κK)
		return ww
	end
	function project_Z!(z,u,x)
		for j = 1:N-1
			u .= z[(j-1)*nxnu+1:(j-1)*nxnu+nu]
			asm.project_U!(u)
			z[(j-1)*nxnu+1:(j-1)*nxnu+nu] .= u

			x .= z[(j-1)*nxnu+nu+1:nu+(j-1)*nxnu+nx]
			asm.project_X!(x)
			z[(j-1)*nxnu+nu+1:nu+(j-1)*nxnu+nx] .= x
		end
	end
	function project_Z(z,u,x)
		zz = deepcopy(z)
		project_Z!(zz,u,x)
		return zz 
	end

	const compute_relerr(z,zopt,nrm2_zopt) = (norm(z .- zopt)^2)/nrm2_zopt
	const compute_relfval(z,fopt) = abs( 0.5 * sum(z .* (asm.G*z)) + asm.h⋅z - fopt )/abs( fopt )
	const compute_rd2K(z,nrm_z_jump,κK) = (norm( ((asm.H*z .- asm.g) .* asm.scl_val_Hg) .- project_K(((asm.H*z .- asm.g) .* asm.scl_val_Hg),κK))^2)/nrm_z_jump

end

module jump
	
	using ..prb
	using ..asm
	using JuMP, MosekTools, ECOS, Gurobi, COSMO, SCS
	using LinearAlgebra

	# useful constants
	const nx = copy(prb.nx)
	const nu = copy(prb.nu)
	const nxnu = nx+nu
	const nxb2 = Int64(nx/2)
	const N = copy(prb.N)

	# initialize main variable
	z = randn(nxnu*(N-1))
	z_xu = randn(nxnu*(N-1))

	function solver_proxy!(z,slvr=:ecos,verbosity=true)
		if slvr == :ecos
			model= Model(ECOS.Optimizer)
			set_optimizer_attribute(model,"printlevel",0)
			set_optimizer_attribute(model,"feastol",prb.ϵ_pd)
			set_optimizer_attribute(model,"abstol",prb.ϵ_gap)
			set_optimizer_attribute(model,"reltol",prb.ϵ_gap)
		elseif slvr == :gurobi
			model= Model(Gurobi.Optimizer)
			set_optimizer_attribute(model,"Presolve",0)
			set_optimizer_attribute(model,"FeasibilityTol",prb.ϵ_pd)
			set_optimizer_attribute(model,"OptimalityTol",prb.ϵ_pd)
			set_optimizer_attribute(model,"BarConvTol",prb.ϵ_gap)
			set_optimizer_attribute(model,"BarQCPConvTol",prb.ϵ_gap)
		elseif slvr == :mosek
			model= Model(Mosek.Optimizer)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_DFEAS",prb.ϵ_pd)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_PFEAS",prb.ϵ_pd)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_REL_GAP",prb.ϵ_gap)
		elseif slvr == :cosmo
			model= Model(COSMO.Optimizer)
			set_optimizer_attribute(model,"eps_abs",prb.ϵ_pd)
			set_optimizer_attribute(model,"eps_rel",prb.ϵ_pd)
			set_optimizer_attribute(model,"eps_prim_inf",prb.ϵ_pd)
			set_optimizer_attribute(model,"eps_dual_inf",prb.ϵ_pd)
		elseif slvr == :scs
			model= Model(SCS.Optimizer)
			set_optimizer_attribute(model,"eps",prb.ϵ_pd)
		else
			error("Invalid solver choice for JuMP.")
		end
		if ~verbosity
			set_silent(model)
		end

		@variable(model,zz[1:nxnu*(N-1)])
		
		asm.set_constraint_jump_vectorized!(model,zz)

		# dynamics
		@constraint(model,asm.H[1:nx*(N-1),1:end]*zz .- asm.g[1:nx*(N-1)] .== 0)

		# affine inequality constraints
		@constraint(model,asm.H[nx*(N-1)+1:end,1:end]*zz .- asm.g[nx*(N-1)+1:end] .≥ 0)

		# cost function
		@objective(model,Min,0.5 * sum(zz .* (asm.G*zz)) + asm.h⋅zz)

		optimize!(model)

		exit_status = string(termination_status(model))
		print("$(string(slvr)) termination status: $(exit_status)\n")
		@assert (exit_status ∈ ("OPTIMAL","SLOW_PROGRESS","ALMOST_OPTIMAL")) "Problem not solved correctly."

		z .= value.(zz)[1:end]
	end

	solver!(solver_name=:ecos,verbosity=true) = solver_proxy!(z,solver_name,verbosity)

	function solver_xu_proxy!(z,slvr=:ecos,verbosity=true)
		if slvr == :ecos
			model= Model(ECOS.Optimizer)
			set_optimizer_attribute(model,"printlevel",0)
			set_optimizer_attribute(model,"feastol",prb.ϵ_pd)
			set_optimizer_attribute(model,"abstol",prb.ϵ_gap)
			set_optimizer_attribute(model,"reltol",prb.ϵ_gap)
		elseif slvr == :gurobi
			model= Model(Gurobi.Optimizer)
			set_optimizer_attribute(model,"Presolve",0)
			set_optimizer_attribute(model,"FeasibilityTol",prb.ϵ_pd)
			set_optimizer_attribute(model,"OptimalityTol",prb.ϵ_pd)
			set_optimizer_attribute(model,"BarConvTol",prb.ϵ_gap)
			set_optimizer_attribute(model,"BarQCPConvTol",prb.ϵ_gap)
		elseif slvr == :mosek
			model= Model(Mosek.Optimizer)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_DFEAS",prb.ϵ_pd)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_PFEAS",prb.ϵ_pd)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_REL_GAP",prb.ϵ_gap)
		elseif slvr == :cosmo
			model= Model(COSMO.Optimizer)
			set_optimizer_attribute(model,"eps_abs",prb.ϵ_pd)
			set_optimizer_attribute(model,"eps_rel",prb.ϵ_pd)
			set_optimizer_attribute(model,"eps_prim_inf",prb.ϵ_pd)
			set_optimizer_attribute(model,"eps_dual_inf",prb.ϵ_pd)
		elseif slvr == :scs
			model= Model(SCS.Optimizer)
			set_optimizer_attribute(model,"eps",prb.ϵ_pd)
		else
			error("Invalid solver choice for JuMP.")
		end
		if ~verbosity
			set_silent(model)
		end

		@variable(model,x[1:nx,1:N])
		@variable(model,u[1:nu,1:N-1])
		@variable(model,β1[1:N-1] ≥ 0)        			# aux variable for defining objective
		@variable(model,β2[1:N-1] ≥ 0)					# aux variable for defining objective

		# dynamics constraint
		Inx = Array(Diagonal(ones(nx)))
		@constraint(model,[t=1:N-1],[Inx -prb.A[t] -prb.B[t] -Inx]*cat(x[:,t+1],x[:,t],u[:,t],prb.e[t],dims=1) .== 0)

		# initial condition
		@constraint(model,x[:,1] .== prb.x0)

		# set problem constraints
		asm.set_constraint_jump_decomposed!(model,x,u) 

		# define objective terms
		@constraint(model,[t=1:N-1],[β1[t],(sqrt(prb.R[t])*(u[:,t] .- prb.uref[t]))...] in SecondOrderCone())
		@constraint(model,[t=1:N-1],[β2[t],(sqrt(prb.Q[t])*(x[:,t+1] .- prb.xref[t]))...] in SecondOrderCone())

		@objective(model,Min,0.5*sum(β1 .^ 2) + 0.5*sum(β2 .^ 2))

		optimize!(model)

		exit_status = string(termination_status(model))
		print("$(string(slvr)) termination status: $(exit_status)\n")
		@assert (exit_status ∈ ("OPTIMAL","SLOW_PROGRESS","ALMOST_OPTIMAL")) "Problem not solved correctly."

		z .= cat([cat(value.(u)[:,t],value.(x)[:,t+1],dims=1) for t=1:N-1]...,dims=1)
	end

	solver_xu!(solver_name=:ecos,verbosity=true) = solver_xu_proxy!(z_xu,solver_name,verbosity)

	compute_nrm_z() = norm(z)^2
	compute_fopt() = 0.5 * sum(z .* (asm.G*z)) + asm.h⋅z

end


# pipg solver and variable containers
module pipg

	using ..prb									# problem definition
	using ..prj
	using ..asm									# assembly of auxiliary parameters
	using ..jump
	using LinearAlgebra, StaticArrays, Plots
	pgfplotsx()

	# useful constants
	const nx = copy(prb.nx)
	const nu = copy(prb.nu)
	const nxnu = nx+nu
	const nH = copy(asm.nH)
	const N = copy(prb.N)

	const μb2σ = 0.5*asm.μ/asm.σ
	const μp2λ = asm.μ+2*asm.λ		

	const nrm_z_jump = jump.compute_nrm_z()
	const fopt = jump.compute_fopt()

	# initialize main variables
	z = randn(nxnu*(N-1))
	v = randn(nx*(N-1)+nH)
	w = randn(nx*(N-1)+nH)	

	zeq = randn(nxnu*(N-1))
	yeq = randn(nx*(N-1)+nH)
	veq = randn(nx*(N-1)+nH)
	weq = randn(nx*(N-1)+nH)	

	# solution stats
	rd2K = randn(prb.kmax)
	rd2o = randn(prb.kmax)
	rfval = randn(prb.kmax)

	rd2K_restrt = randn(prb.kmax)
	rd2o_restrt = randn(prb.kmax)
	rfval_restrt = randn(prb.kmax)
	
	rd2K_eq = randn(prb.kmax)
	rd2o_eq = randn(prb.kmax)
	rfval_eq = randn(prb.kmax)

	# temporary variables
	κK = randn(nH) 
	const κx = MVector{nx}(randn(nx))
	const κu = MVector{nu}(randn(nu))

	# pipg+ solver (version 1)
	function solver_v1_proxy!(z,v,w,κK,κu,κx,rd2o,rd2K,rfval,restart_idx=prb.kmax)

		for k = 1:prb.kmax
			# βᵏ = (k+1)*μb2σ
			βᵏ = (mod(k,restart_idx)+1)*μb2σ
			
			# feasibility measure
			rd2K[k] = prj.compute_rd2K(z,nrm_z_jump,κK)

			# relative objective value decrease
			# rfval[k] = prj.compute_relfval(z,fopt)
			rfval[k] = abs(z⋅(asm.G*z)+z⋅asm.h - jump.z⋅(asm.G*jump.z)-jump.z⋅asm.h)/abs(jump.z⋅(asm.G*jump.z)+jump.z⋅asm.h)

			v .= w .- βᵏ .* (asm.H*z .- asm.g)
			prj.project_Kdual!(v,κK)
			
			# αᵏ = 2/(k*asm.μ+μp2λ)
			αᵏ = 2/(mod(k,restart_idx)*asm.μ+μp2λ)

			z .= z .- αᵏ .* (asm.G*z .+ asm.h .- asm.HT*v)
			prj.project_Z!(z,κu,κx)

			w .= w .- βᵏ .* (asm.H*z .- asm.g)
			# prj.project_Kdual!(w,κK)

			rd2o[k] = prj.compute_relerr(z,jump.z,nrm_z_jump)

		end
		if restart_idx == prb.kmax
			print("PIPG+ relative distance to optimum: $(prj.compute_relerr(z,jump.z,nrm_z_jump))\n")
		else
			print("PIPG+ restarted w/ $(restart_idx) relative distance to optimum: $(prj.compute_relerr(z,jump.z,nrm_z_jump))\n")
		end	

	end

	# pipg+ solver
	function solver_proxy!(z,v,w,κK,κu,κx,rd2o,rd2K,rfval,restart_idx=prb.kmax)

		zkm1 = randn(length(z))
		for k = 1:prb.kmax
			# βᵏ = (k+1)*μb2σ
			βᵏ = (mod(k,restart_idx)+1)*μb2σ
			
			# feasibility measure
			rd2K[k] = prj.compute_rd2K(z,nrm_z_jump,κK)

			# relative objective value decrease
			# rfval[k] = prj.compute_relfval(z,fopt)
			rfval[k] = abs(z⋅(asm.G*z)+z⋅asm.h - jump.z⋅(asm.G*jump.z)-jump.z⋅asm.h)/abs(jump.z⋅(asm.G*jump.z)+jump.z⋅asm.h)

			w .= v .+ βᵏ .* (asm.H*z .- asm.g)
			prj.project_Kpolar!(w,κK)
			
			# αᵏ = 2/(k*asm.μ+μp2λ)
			αᵏ = 2/(mod(k,restart_idx)*asm.μ+μp2λ)

			zkm1 .= z
			z .= z .- αᵏ .* (asm.G*z .+ asm.h .+ asm.HT*w)
			prj.project_Z!(z,κu,κx)

			v .= w .+ βᵏ .* (asm.H)*(z .- zkm1)

			rd2o[k] = prj.compute_relerr(z,jump.z,nrm_z_jump)

		end
		if restart_idx == prb.kmax
			print("PIPG+ relative distance to optimum: $(prj.compute_relerr(z,jump.z,nrm_z_jump))\n")
		else
			print("PIPG+ restarted w/ $(restart_idx) relative distance to optimum: $(prj.compute_relerr(z,jump.z,nrm_z_jump))\n")
		end	

	end

	solver!(restart_idx=prb.kmax) = solver_proxy!(z,v,w,κK,κu,κx,rd2o,rd2K,rfval,restart_idx)

	# pipg solver (non-strongly convex formulation; affine equality constraint replaces the affine conic constraint)
	function solver_eq_proxy!(z,y,v,w,κK,κu,κx,rd2o,rd2K,rfval)
		β = 2.0
		α = 1/(β*asm.σ_eq+asm.λ)
		for k = 1:prb.kmax

			# feasibility measure
			rd2K[k] = prj.compute_rd2K(z,nrm_z_jump,κK)

			# relative objective value decrease
			# rfval[k] = prj.compute_relfval(z,fopt)
			rfval[k] = abs(z⋅(asm.G*z)+z⋅asm.h - jump.z⋅(asm.G*jump.z)-jump.z⋅asm.h)/abs(jump.z⋅(asm.G*jump.z)+jump.z⋅asm.h)

			v .= w .+ β .* (asm.H*z .- y .- asm.g)

			z .= z .- α .* (asm.G*z .+ asm.h .+ asm.HT*v)	
			prj.project_Z!(z,κu,κx)

			y .= y .+ α .* v
			prj.project_K!(y,κK)

			w .= w .+ β .* (asm.H*z .- y .- asm.g)

			rd2o[k] = prj.compute_relerr(z,jump.z,nrm_z_jump)
		end

		print("PIPGeq relative distance to optimum: $(prj.compute_relerr(z,jump.z,nrm_z_jump))\n")

	end

	solver_eq!() = solver_eq_proxy!(zeq,yeq,veq,weq,κK,κu,κx,rd2o_eq,rd2K_eq,rfval_eq)

	plot_solstat() = begin plt_d2o = plot(1:prb.kmax,rd2o,line=:solid,color=:magenta,title="Relative distance to optimum",lab="PIPG",lw=1.1,yaxis=:log,xlabel="Iterations") 
			plot!(1:prb.kmax,rd2o_restrt,line=:solid,color=:red,lab="PIPG restarted",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rd2o_eq,line=:solid,color=:purple,lab="PIPGeq",lw=1.1,yaxis=:log)
			plot!(size=(400,500),legend=:none)
			plt_fval = plot(1:prb.kmax,rfval,line=:solid,color=:magenta,title="Relative objective value",lab="PIPG",lw=1.1,yaxis=:log,xlabel="Iterations") 
			plot!(1:prb.kmax,rfval_restrt,line=:solid,color=:red,lab="PIPG restarted",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rfval_eq,line=:solid,color=:purple,lab="PIPGeq",lw=1.1,yaxis=:log)
			plot!(size=(400,500),legend=:none)
			plt_feas = plot(1:prb.kmax,rd2K,line=:solid,color=:magenta,title="Relative feasibility",lab="PIPG",lw=1.1,xlabel="Iterations",yaxis=:log)
			plot!(1:prb.kmax,rd2K_restrt,line=:solid,color=:red,lab="PIPG restarted",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rd2K_eq,line=:solid,color=:purple,lab="PIPGeq",lw=1.1,yaxis=:log)
			plot!(size=(400,500))
			λ = @layout [a b c]
			display(plot(plt_d2o,plt_fval,plt_feas,layout=λ,size=(4*1200/5,4*500/5)))
	end

	function reset_var!(flg=:both)
		if flg == :plus
			z .= randn(nxnu*(N-1))
			v .= randn(nx*(N-1)+nH)
			w .= randn(nx*(N-1)+nH)
		elseif flg == :eq
			zeq .= randn(nxnu*(N-1))
			yeq .= randn(nx*(N-1)+nH)
			veq .= randn(nx*(N-1)+nH)
			weq .= randn(nx*(N-1)+nH)
		elseif flg == :both
			z .= randn(nxnu*(N-1))
			v .= randn(nx*(N-1)+nH)
			w .= randn(nx*(N-1)+nH)

			zeq .= randn(nxnu*(N-1))
			yeq .= randn(nx*(N-1)+nH)
			veq .= randn(nx*(N-1)+nH)
			weq .= randn(nx*(N-1)+nH)
		else
			error("Invalid reset flag.")
		end
	end

end

# rival solvers
module rival

	using ..prb									# problem definition
	using ..prj
	using ..asm									# assembly of auxiliary parameters
	using ..jump
	using LinearAlgebra, StaticArrays, Plots, LaTeXStrings
	gr()
	default(xtickfontsize=11,ytickfontsize=11,ztickfontsize=11,legendfontsize=9,labelfontsize=11)

	# useful constants
	const nx = copy(prb.nx)
	const nu = copy(prb.nu)
	const nxnu = nx+nu
	const nH = copy(asm.nH)
	const N = copy(prb.N)

	const nrm_z_jump = jump.compute_nrm_z()
	const fopt = jump.compute_fopt()

	# initialize main variables
	z = randn(nxnu*(N-1))
	x = randn(nxnu*(N-1))
	y = randn(nx*(N-1)+nH)
	v = randn(nx*(N-1)+nH)	
	w = randn(nx*(N-1)+nH)	

	# solution stats
	rd2K = randn(prb.kmax)
	rd2o = randn(prb.kmax)
	rfval = randn(prb.kmax)
	nest_iter_count = zeros(prb.kmax)

	# temporary variables
	# temporary variables
	κK = randn(nH) 
	const κx = MVector{nx}(randn(nx))
	const κu = MVector{nu}(randn(nu))
	κz = randn(nxnu*(N-1))

	# ADMM
	function solver_admm_proxy!(z,x,y,v,w,κz,κK,κx,κu,rd2o,rd2K,rfval,nest_iter_count)
		α = copy(prb.α_admm)
		λ̃ = asm.λ + α*asm.σ
		λ̃inv = 1/λ̃
		β = (sqrt(λ̃) - sqrt(asm.μ))/(sqrt(λ̃) + sqrt(asm.μ))
		
		inner_term_crit = 100.0
		k_nest = 0
		for k=1:prb.kmax

			# feasibility measure
			rd2K[k] = prj.compute_rd2K(z,nrm_z_jump,κK)

			# relative objective value decrease
			# rfval[k] = prj.compute_relfval(z,fopt)
			rfval[k] = abs(z⋅(asm.G*z)+z⋅asm.h - jump.z⋅(asm.G*jump.z)-jump.z⋅asm.h)/abs(jump.z⋅(asm.G*jump.z)+jump.z⋅asm.h)

			while inner_term_crit ≥ prb.ϵ_nest_admm
				v .= asm.H*x .- y .- asm.g .+ w

				κz .= z 
				z .= x .- λ̃inv .* (asm.G*x .+ asm.h .+ α .* asm.HT*v)
				prj.project_Z!(z,κu,κx)	


				x .= z .+ β .* (z .- κz)

				k_nest > prb.kmax_nest && break

				inner_term_crit = norm(z .- x)/norm(z)
				k_nest+=1
				nest_iter_count[k] += 1
			end	
			k_nest = 0
			inner_term_crit = 100.0

			y .= asm.H*z .- asm.g .+ w
			prj.project_K!(y,κK)
			w .= w .+ asm.H*z .- y .- asm.g

			rd2o[k] = prj.compute_relerr(z,jump.z,nrm_z_jump)
		end
		print("ADMM relative distance to optimum: $(prj.compute_relerr(z,jump.z,nrm_z_jump))\n")
	end

	solver_admm!() = solver_admm_proxy!(z,x,y,v,w,κz,κK,κx,κu,rd2o,rd2K,rfval,nest_iter_count)

	z_2 = randn(nxnu*(N-1))
	y_2 = randn(nxnu*(N-1))
	w_2 = randn(nx*(N-1)+nH)	

	# solution stats
	rd2K_2 = randn(prb.kmax)
	rd2o_2 = randn(prb.kmax)
	rfval_2 = randn(prb.kmax)

	# Chambolle and Pock variable step size
	function solver_cp_proxy!(z,y,w,κK,κx,κu,rd2o,rd2K,rfval)

	    βᵏ = prb.β⁰
	    γᵏ = prb.γ⁰
	    αᵏ = 1/(asm.λ+asm.σ*βᵏ)

		for k=1:prb.kmax

			# feasibility measure
			rd2K[k] = prj.compute_rd2K(z,nrm_z_jump,κK)

			# relative objective value decrease
			# rfval[k] = prj.compute_relfval(z,fopt)
			rfval[k] = abs(z⋅(asm.G*z)+z⋅asm.h - jump.z⋅(asm.G*jump.z)-jump.z⋅asm.h)/abs(jump.z⋅(asm.G*jump.z)+jump.z⋅asm.h)
        
	        w .= -1 .* prj.project_Kdual( -1 .* w .- βᵏ .* (asm.H*(z .+ γᵏ .* (z .- y)) .- asm.g),κK)
	        y .= z
	        α_step = αᵏ/(asm.μ*αᵏ+1)
	        z .= prj.project_Z(z .- α_step .* (asm.G*z .+ asm.h .+ asm.HT*w),κu,κx)
	        
	        # step-size update
	        term_1 = (0.5*asm.λ*αᵏ + sqrt((1+asm.μ*αᵏ)*(1-asm.λ*αᵏ) + (0.5*asm.λ*αᵏ)^2 ))
	        γᵏ = term_1/(1+asm.μ*αᵏ)
	        αᵏ = αᵏ/term_1
	        βᵏ = (-asm.λ + 1/αᵏ)/asm.σ

			rd2o[k] = prj.compute_relerr(z,jump.z,nrm_z_jump)
		end
		print("Chambolle and Pock relative distance to optimum: $(prj.compute_relerr(z,jump.z,nrm_z_jump))\n")
	end

	solver_cp!() = solver_cp_proxy!(z_2,y_2,w_2,κK,κx,κu,rd2o_2,rd2K_2,rfval_2)

	function reset_var!(flg::Symbol)
		if flg == :admm
			z .= randn(nxnu*(N-1))
			x .= randn(nxnu*(N-1))
			y .= randn(nx*(N-1)+nH)
			v .= randn(nx*(N-1)+nH)	
			w .= randn(nx*(N-1)+nH)	
		elseif flg == :cp
			z_2 .= randn(nxnu*(N-1))
			y_2 .= randn(nxnu*(N-1))
			w_2 .= randn(nx*(N-1)+nH)
		else
			error("Invalid reset flag.")
		end
	end

	using ..pipg

	plot_solstat() = begin plt_d2o = plot(1:prb.kmax,pipg.rd2o,line=:solid,color=:magenta,title=L"\frac{\|z-z^{\star}\|_2^2}{\|z^{\star}\|_2^2}",lab="PIPG",lw=1.1,yaxis=:log,xlabel="Iterations") 
			plot!(1:prb.kmax,pipg.rd2o_restrt,line=:solid,color=:red,lab="PIPG w/ restart",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,pipg.rd2o_eq,line=:solid,color=:purple,lab="PIPGeq",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rd2o,line=:solid,color=:green,lab="ADMM",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rd2o_2,line=:solid,color=:orange,lab="PDHG",lw=1.1,yaxis=:log)
			plot!(size=(400,500),legend=:none,xlim=[1,prb.kmax],ylim=[1e-15,1e1])

			plt_rfval = plot(1:prb.kmax,pipg.rfval,line=:solid,color=:magenta,title=L"\frac{|f(z)-f(z^{\star})|}{|f(z^{\star})|}",lab="PIPG",lw=1.1,yaxis=:log,xlabel="Iterations") 
			plot!(1:prb.kmax,pipg.rfval_restrt,line=:solid,color=:red,lab="PIPG w/ restart",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,pipg.rfval_eq,line=:solid,color=:purple,lab="PIPGeq",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rfval,line=:solid,color=:green,lab="ADMM",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rfval_2,line=:solid,color=:orange,lab="PDHG",lw=1.1,yaxis=:log)
			plot!(size=(400,500),legend=:none,xlim=[1,prb.kmax],ylim=[1e-15,1e1])

			plt_nest_iter = plot(1:prb.kmax,nest_iter_count,line=:solid,color=:purple,title="Nesterov iteration count in ADMM",lw=1.1,xlabel="Outer iterations",lab="")
			plot!(1:prb.kmax,prb.kmax_nest .* ones(prb.kmax),line=:dash,lw=1.5,lab="Max iterations",color=:black,size=(1200,500))

			plt_feas = plot(1:prb.kmax,pipg.rd2K,line=:solid,color=:magenta,title=L"\frac{d_{K}(Hz-g)}{\|z^{\star}\|_2^2}",lab="PIPG",lw=1.1,xlabel="Iterations",yaxis=:log)
			plot!(1:prb.kmax,pipg.rd2K_restrt,line=:solid,color=:red,lab="PIPG w/ restart",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,pipg.rd2K_eq,line=:solid,color=:purple,lab="PIPGeq",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rd2K,line=:solid,color=:green,lab="ADMM",lw=1.1,yaxis=:log)
			plot!(1:prb.kmax,rd2K_2,line=:solid,color=:orange,lab="PDHG",lw=1.1,yaxis=:log)	
			plot!(size=(400,500),legend=:topright,xlim=[1,prb.kmax],ylim=[1e-15,1e1])

			λ = @layout [a b c ; d]
			display(plot(plt_d2o,plt_rfval,plt_feas,plt_nest_iter,layout=λ,size=(4*1200/5,4*1000/5)))
	end

end