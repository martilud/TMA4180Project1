using Plots

mutable struct instance
  n::Int
  lengths::Array{Float64}
  angles::Array{Float64}
  point::Array{Float64}
end

mutable struct point
  x::Float64
  y::Float64
end

ex = instance(3, [3,2,2], [0,0,0], [3,2]) 

function draw(instance)
  n = instance.n
  lengths = instance.lengths
  angles = instance.angles
  curr = point(0,0)
  for i = 1:n
    next = point(sin(angles[i])*lengths[i], cos(angles[i])*lengths[i])
    plot([curr.x, next.x], [curr.y, next.y], show = true)
    curr = next
  end
end
draw(ex)
