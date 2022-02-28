module proj
using LinearAlgebra
# ! indicates that the variable is modified

function infnorm(x,n,γ)
# compute ∞-norm of x
# ---
# x : input vector 				(real n-vector)
# n : length of x
# γ : temp storage variable		(real scalar)    !

	γ = 0.0
	for j=1:n
		γ = max(γ,abs(x[j]))
	end
	return γ
end

function hlfspace!(x,z,u,ζ,n,β)
# projection onto halfspace { x | u⋅x ≤ ζ }
# ---
# x : projected result					(real n-vector) !
# z : vector to be projected			(real n-vector) 
# u : normal to supporting hyperplane 	(real n-vector)
# ζ : scalar on RHS 					(real scalar)
# β : temporary variable                (real scalar)   !

	@assert BLAS.nrm2(n,u,1)≈1 "Normal vector u should have unit norm."

	β = BLAS.dot(n,z,1,u,1) - ζ

	if β > 0 
		BLAS.blascopy!(n,z,1,x,1)
		axpy!(-β,u,x)
	else
		BLAS.blascopy!(n,z,1,x,1)
	end

end

function hlfspace2!(x,z,u1,u2,ζ1,ζ2,n,β,α1,α2)
# projection onto intersection of two halfspaces { x | u1⋅x ≤ ζ1 } ∩ { x | u2⋅x ≤ ζ2 }
# ---
# x  : projected result					(real n-vector) !
# z  : vector to be projected			(real n-vector) 
# u1 : normal to supporting hyperplane 	(real n-vector)
# u2 : normal to supporting hyperplane 	(real n-vector)
# ζ1 : scalar on RHS 					(real scalar)
# ζ1 : scalar on RHS 					(real scalar)
# β  : temporary variable               (real scalar)   !
# α1 : temporary variable               (real scalar)   !
# α2 : temporary variable               (real scalar)   !

	@assert BLAS.nrm2(n,u1,1)≈1 "Normal vector u1 should have unit norm." 
	@assert BLAS.nrm2(n,u2,1)≈1 "Normal vector u2 should have unit norm." 

	β = BLAS.dot(n,u1,1,u2,1)
	α1 = BLAS.dot(n,z,1,u1,1)
	α2 = BLAS.dot(n,z,1,u2,1)

	if β ≈ 1
		β = min(ζ1,ζ2)
		if α1 > β
			BLAS.blascopy!(n,z,1,x,1)
			axpy!(β-α1,u1,x)
		else
			BLAS.blascopy!(n,z,1,x,1)
		end
	elseif β ≈ -1
		@assert ζ1+ζ2>0 "The two halfspaces do not intersect."
		if α1 < -ζ2
			BLAS.blascopy!(n,z,1,x,1)
			axpy!(-α1-ζ2,u1,x)		
		elseif α1 > ζ1
			BLAS.blascopy!(n,z,1,x,1)
			axpy!(-α1+ζ1,u1,x)
		else
			BLAS.blascopy!(n,z,1,x,1)
		end
	else
		if α1-ζ1+β*(ζ2-α2)>0 && α2-ζ2+β*(ζ1-α1)>0
			BLAS.blascopy!(n,z,1,x,1)
			axpy!((α1-ζ1+β*(ζ2-α2))/(β*β-1),u1,x)
			axpy!((α2-ζ2+β*(ζ1-α1))/(β*β-1),u2,x)
		elseif α2>ζ2 && α1-ζ1+β*(ζ2-α2)≤0
			BLAS.blascopy!(n,z,1,x,1)
			axpy!(ζ2-α2,u2,x)
		elseif α1>ζ1 || α2>ζ2
			BLAS.blascopy!(n,z,1,x,1)
			axpy!(ζ1-α1,u1,x)
		else
			BLAS.blascopy!(n,z,1,x,1)
		end
	end

end

function box!(x,z,l,u,n)
# projection onto l∞-norm ball
# ---
# x : projected result			(real n-vector) !
# z : vector to be projected	(real n-vector) 
# l : lower bound			 	(real scalar)
# u : upper bound 				(real scalar)
# n : length of z 				(integer)

	@assert u>l "Upper bound should be greater than lower bound."

	for j = 1:n
		x[j] = min(u,max(z[j],l))
	end

end

function ball!(x,z,r,n)
# projection onto l2-norm ball
# ---
# x : projected result			(real n-vector) !
# z : vector to be projected	(real n-vector) 
# r : radius of l2-norm ball 	(real scalar)
# n : length of z 				(integer)
	
	@assert r>0 "Radius should be positive."

	BLAS.blascopy!(n,z,1,x,1)
	BLAS.scal!(n,r/max(BLAS.nrm2(n,x,1),r),x,1)

end

function soc!(x,z,α,n,β)
# projection onto second-order cone (aka soc, ice-cream, lorentz cone) 
# type 1 : cone axis along last coordinate and arctan(α) is the cone angle
# ---
# x : projected result 		     (real n-vector) !
# z : vector to be projected	 (real n-vector)
# α : tan(cone angle)			 (real scalar)
# n : length of z (integer)		 (integer)	
# β : temp storage variable  	 (real scalar)   !

	@assert α>0 "Cone angle arctan(α) must be positive."

	β = BLAS.nrm2(n-1,z,1)

	if β ≤ α*z[n]
		BLAS.blascopy!(n,z,1,x,1)
	elseif α*β + z[n] ≤ 0
		fill!(x,0.0)
	else
		BLAS.blascopy!(n-1,z,1,x,1)
		x[n] = (β*α+z[n])/(α*α+1)
		BLAS.scal!(n-1,α*x[n]/β,x,1)
	end

end

function soc2!(x,z,θ1,θ2,u,n,β,γ,κ)
# projection onto second-order cone (aka soc, ice-cream, lorentz cone) 
# type 2 : cone axis along unit vector u and arccos(θ) is the cone angle
# ---
# x  : projected result 			(real n-vector) !
# z  : vector to be projected		(real n-vector)
# θ1 : cos(cone angle) 				(real scalar)
# θ2 : sin(cone angle) 				(real scalar)
# u  : unit vector along cone axis	(real n-vector)
# n  : length of z 					(integer)
# β  : temp storage variable 		(real scalar)   !
# γ  : temp storage variable 		(real scalar)   !
# κ  : temp storage variable 		(real n-vector) !

	@assert min(θ1,1-θ1)>0 && min(θ2,1-θ2)>0 "θ1,θ2 ∈ [0,1]"
	@assert BLAS.nrm2(n,u,1)≈1 "u must be a unit vector"

	β = BLAS.nrm2(n,z,1)
	γ = BLAS.dot(n,z,1,u,1)

	if θ1*β ≤ γ
		BLAS.blascopy!(n,z,1,x,1)
	elseif β*θ2 + γ ≤ 0
		fill!(x,0.0)
	else
		BLAS.blascopy!(n,u,1,x,1)
		BLAS.blascopy!(n,z,1,κ,1)
		axpy!(-γ,u,κ)
		axpby!(θ2/BLAS.nrm2(n,κ,1),κ,θ1,x)
		BLAS.scal!(n,BLAS.dot(n,z,1,x,1),x,1)
	end

end

function soc_ball!(x,z,α,r,n,β)
# projection onto intersection of soc cone and l2-norm ball
# ---
# x : projected result 			(real n-vector) !
# z : vector to be projected	(real n-vector)
# α : tan(cone angle) 			(real scalar)
# r : radius of l2-norm ball 	(real scalar)
# n : length of z 				(integer)
# β : temp storage variable 	(real scalar)   !

	soc!(x,z,α,n,β)
	ball!(x,x,r,n)

end

function soc2_ball!(x,z,θ1,θ2,u,r,n,β,γ,κ)
# projection onto intersection of soc cone and l2-norm ball
# ---
# x  : projected result 			(real n-vector) !
# z  : vector to be projected		(real n-vector)
# θ1 : cos(cone angle) 				(real scalar)
# θ2 : cos(cone angle) 				(real scalar)
# u  : unit vector along cone axis	(real n-vector)
# r  : radius of l2-norm ball 		(real scalar)
# n  : length of z 					(integer)
# β  : temp storage variable 		(real scalar)   !
# γ  : temp storage variable 		(real scalar)   !
# κ  : temp storage variable 		(real n-vector) !

	soc2!(x,z,θ1,θ2,u,n,β,γ,κ)
	ball!(x,x,r,n)

end

function admm!(x,y,u,z,f1,f2,α,β,ϵ,n,γ1,γ2,κ1,κ2)
# projection onto intersection of two sets with respective projections defined by f1 and f2
# ---
# x  : projected result (copy 1)	(real n-vector) !
# y  : projected result (copy 2)	(real n-vector) !
# u  : dual variable 				(real n-vector) !
# z  : vector to be projected		(real n-vector)
# f1 : projection function			(function: f!(θ,ψ) where θ is modified)
# f2 : projection function 			(function: f!(θ,ψ) where θ is modified)
# α  : ADMM step size				(real scalar)
# β  : 1/(1+α) 						(real scalar)
# ϵ  : ADMM termination criteria  	(real scalar)
# n  : length of z					(integer)
# γ  : temp storage variable        (real scalar)   !
# γ2 : temp storage variable 		(real scalar)   !
# κ1 : temp storage variable 		(real n-vector) !
# κ2 : temp storage variable 		(real n-vector) !
	
# x, y and u may hold values from previous call to admm! (warm-starting)


	f1(κ1,z)
	f2(κ2,z)
	κ1 .= κ1 .- z
	κ2 .= κ2 .- z
	γ1 = BLAS.nrm2(n,κ2,1)
	γ1 = max(γ1,BLAS.nrm2(n,κ1,1))
	if γ1 ≤ ϵ
	    BLAS.blascopy!(n,z,1,x,1)
		BLAS.blascopy!(n,z,1,y,1)
		# κ2[1] will be zero fortunately
		κ2[1] = 0.0 # inserted for type-stability
	else
		κ1 .= z .- u
		axpy!(α,y,κ1)
		BLAS.scal!(n,β,κ1,1)
		f1(x,κ1)
		κ2 .= z .+ u
		axpy!(α,x,κ2)
		BLAS.scal!(n,β,κ2,1)
		f2(y,κ2)
		κ1 .= x .- y
		axpy!(α,κ1,u)

		κ2[1] = infnorm(κ1,n,γ2)
		γ1 = 1.0
		while κ2[1]≥ϵ
			κ1 .= z .- u
			axpy!(α,y,κ1)
			BLAS.scal!(n,β,κ1,1)
			f1(x,κ1)
			κ2 .= z .+ u
			axpy!(α,x,κ2)
			BLAS.scal!(n,β,κ2,1)
			f2(y,κ2)
			κ1 .= x .- y
			axpy!(α,κ1,u)

			γ1 = γ1 + 1.0 # iteration counter

			κ2[1] = infnorm(κ1,n,γ2)
		end
		κ2[1] = γ1 # store the iteration count
	end

end

function multi_admm!(x,y,u,z,f,α,β,ϵ,n,N,γ1,γ2,ψ,κ)
# projection onto intersection of N sets with projections defined by f[i], i=1,...,N
# ---
# x  : projected results (copy 1)	(N-vector of real n-vectors) !
# y  : projected results (copy 2)	(N-vector of real n-vectors) !
# u  : dual variable 				(N-vector of real n-vectors) !
# z  : vector to be projected		(real n-vector)
# f  : N projection functions   	(N-vector of functions: f[i]!(θ,ψ) where θ is modified, i=1,...,N)
# α  : ADMM step-size				(positive scalar)
# β  : 1/(2α(N-1)+1)				(positive scalar)
# n  : length of z 					(integer)
# N  : no. of sets					(integer)
# γ1 : temp storage variable  		(real scalar)                !
# γ2 : temp storage variable  		(real scalar)                !
# ψ  : temp storage variable  		(real scalar)                !
# κ  : temp storage variable        (real n-vector)              !

	@assert length(f)==N "No. of projection functions is not N."
	# @assert N ≥ 3 "Use admm! for fewer than 3 sets."

	# γ1 = 100.0
	# ψ = 0.0 # iteration counter

	for i = 1:N
		κ .= x[i] # BLAS.blascopy!(n,x[i],1,κ,1)
		lmul!(N-2,κ)
		for j = 1:N
			κ .= κ .+ x[j]  
		end
		lmul!(α,κ)
		κ .= κ .- u[i] .+ z
		lmul!(β,κ) 
		f[i](y[i],κ)
	end

	γ1 = 0.0 # holds termination criteria value
	x[1] .= y[1]
	for i = 2:N
		x[i] .= y[i]
		κ .= x[1] .- x[i]

		γ1 = max(γ1,infnorm(κ,n,γ2))
	end

	for i = 1:N
		κ .= x[i]
		lmul!(N,κ)
		for j = 1:N
			κ .= κ .- x[j]
		end
		axpy!(α,κ,u[i])
	end
	
	ψ = 1.0

	while γ1≥ϵ 
		for i = 1:N
			κ .= x[i] # BLAS.blascopy!(n,x[i],1,κ,1)
			lmul!(N-2,κ)
			for j = 1:N
				κ .= κ .+ x[j]  
			end
			lmul!(α,κ)
			κ .= κ .- u[i] .+ z
			lmul!(β,κ) 
			f[i](y[i],κ)
		end

		γ1 = 0.0 # holds termination criteria value
		x[1] .= y[1]
		for i = 2:N
			x[i] .= y[i]
			κ .= x[1] .- x[i]

			γ1 = max(γ1,infnorm(κ,n,γ2))
		end

		for i = 1:N
			κ .= x[i]
			lmul!(N,κ)
			for j = 1:N
				κ .= κ .- x[j]
			end
			axpy!(α,κ,u[i])
		end
		
		ψ = ψ + 1.0
	end 
	κ[1] = ψ

end

end