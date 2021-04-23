"""
	W1_1D(a, b)

Compute the 1-Wasserstein distance between measures `a` and `b`, assuming
that both are supported on the points 1:N.
"""
function W1_1D(a, b)
    length(a) == length(b) || throw(DimensionMismatch("lengths are not compatible"))
    W1_1D_unsafe(a, b)
end

function W1_1D_unsafe(a, b)
	s = 0.0
    w = 0.0
    for i in eachindex(a)
        s += @inbounds a[i] - b[i]
        w += abs(s)
    end
    return w
end

# TODO: Doesn't give the same result as W1 for unnormalized data, why?
function solve_convex_OT_plan(μ, ν, x, y, c)
	_solve_convex_OT(μ, ν, x, y, c, Val(true))
end

function solve_convex_OT_cost(μ, ν, x, y, c)
	_solve_convex_OT(μ, ν, x, y, c, Val(false))
end

#TODO: make this return a sparse matrix
function _solve_convex_OT(μ, ν, x, y, c, return_plan::Val{T}) where T
	m = length(x)
	n = length(y)
	massνj = 0.0

	i = 1
	massμi = μ[1]
	if T
		plan = zeros(n, m)
	end
	cost = 0.0

	for j in eachindex(ν)
		massνj = ν[j]
		while massνj > 0
			mass = min(massμi, massνj)
			massμi -= mass
			massνj -= mass
			if T
				plan[i,j] += mass
			end
			cost += mass*c(x[i], y[j])
			if massμi==0
				if i < m
					i+=1
					massμi = μ[i]
				else
					T ? (return (cost, plan)) : (return cost)
				end
			end
		end
	end
	T ? (return (cost, plan)) : (return cost)
end
