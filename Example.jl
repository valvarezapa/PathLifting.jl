include("PathLifting.jl")

using PyPlot
using Images
using Polynomials
using BenchmarkTools
using FileIO
using Polynomials
using DelimitedFiles
using Plots

# Polygonal path
#=
points=complex([0.02852374436962157,0.9583268689378289,0.040248228476352135+0.06865542334894827*im,0.028521558015962945,0.040248228476352135-0.06865542334894827*im,-0.0035557874067469165+0.003299155231960405*im,-0.0035557874067469165-0.003299155231960405*im])
radius=PL.min_dist_max(points)/3.0

dist=PL.dist_max(complex(0.9583290552914876),points[1]+complex(radius))
print("Distancia que recorre el camino: ");println(dist)
step=minimum([PL.dist_max(complex(0.9583290552914876),points),PL.dist_max(points[1]+complex(radius),points),PL.min_dist_max(points)])/3.0
print("Tamaño de los pasos: ");println(step)
print("Relacion: ");println(step/dist)

PL.polygonal_path(points,complex(0.9583290552914876),points[1]+complex(radius))
=#

# Plot polygonal path
#PL.plot_polygonal_path([complex(1.0)],complex(-5.0),complex(6.0))
#gcf()
#clf()

# Plot loops
#=
points=complex([-2.0,2.0])
loops=PL.discrete_loops(points)
l=length(loops)
@inbounds for i in 1:l
    PL.plot_discrete_loop(points,loops[i])
end
gcf()
=#
#clf()


# Plot lifting of a path
#=
pol=complex([-1.0,0.0,1.0])
path=PL.polygonal_path([complex(-1.0)],complex(2.0),complex(-2.0))
starting_point=3.0*im
PL.precise_lifting_plot(pol,path,starting_point;preci=2.0)
gcf()
=#

# Generators initial point
#=
co=complex([0.0,-3.0,0.0,1.0])
a=PL.generators_initial_point(co,complex(0.0),complex([-2.0,2.0]))
print(a)
=#

# new N-solve aux
#=
co=complex([0.02222222222222222,-0.04444444444444444,-0.06666666666666667,0.08888888888888888,0.11111111111111109,-0.13333333333333333,0.15555555555555553,-0.17777777777777776,-0.1999999999999996])
gens=PL.new_nsolve_aux(co,complex(0.0),complex([0.02852374436962157,0.9583268689378289,0.040248228476352135+0.06865542334894827*im,0.028521558015962945,0.040248228476352135-0.06865542334894827*im,-0.0035557874067469165+0.003299155231960405*im,-0.0035557874067469165-0.003299155231960405*im]))
=#

# Non critical N-solve aux
#=
co=complex([0.02222222222222222,-0.04444444444444444,-0.06666666666666667,0.08888888888888888,0.11111111111111109,-0.13333333333333333,0.15555555555555553,-0.17777777777777776,-0.1999999999999996])
a=PL.noncritical_nsolve_aux(co,complex(0.0),complex([0.02852374436962157,0.9583268689378289,0.040248228476352135+0.06865542334894827*im,0.028521558015962945,0.040248228476352135-0.06865542334894827*im,-0.0035557874067469165+0.003299155231960405*im,-0.0035557874067469165-0.003299155231960405*im]))
print(a)
=#

# N-solve aux
#=
@inbounds init_pol=[complex(1.0) for i in 1:9+1]
# 112, 273, 174
mod=PL.binary(517)
PL.sameLength(mod,init_pol)
pol=init_pol
@inbounds for i in 1:9+1
    if mod[i]==true
     pol[i]=complex(-1.0)
    end
end
=#

#co=complex([0.02222222222222222,-0.04444444444444444,-0.06666666666666667,0.08888888888888888,0.11111111111111109,-0.13333333333333333,0.15555555555555553,-0.17777777777777776,-0.1999999999999996])
#a=PL.discrete_loop(complex([-0.07934558977667738-0.015617398369028439*im,4.198591731440788,-0.07934558977667738+0.01561739836902844*im,-0.11567222702033735+0.1287753983441095*im,-0.25476495925382336,-0.046331037899165456+0.00941611102652034*im,-0.046331037899165456-0.009416111026520341*im,-0.11567222702033737-0.1287753983441095*im]),complex(-0.25476495925382336),0.00627740735101356)
#print(a)


# N-solve
#=
co=complex([0.0,0.0,-1.0,1.0])
PL.nsolve(co,complex(0.0);pre=8.0)
=#

# Zeroes mult
#=
co=complex([0.0,0.0,0.0,1.0])
a=PL.zeroes_mult(co;precision=14.0)
=#

# Compute the roots of every Littlewood polynomial of degree n, and save them into a file.
#PL.littlewood_roots(9)

# Plot every computed root of Littlewood polynomials
#PL.plot_littlewood((-1.5,1.5),(-1.5,1.5))
#PyPlot.savefig(fname="plot.png")
#println("Figure saved succesfully!")

# Julia set of a polynomial
#PL.julia_set(complex([1.618033988749,0.0,1.0]),complex([0.0,20,0]),16)

# Lifting of the fundamental group
# pol 3 roots: complex([-2.0,0.0,1.0,1.0])
# pol 3 roots through an affine map: complex([119.0-85.0*im,-151.0-217.0*im,-130.0+90.0*im,18.0+26.0*im])
# pol 5 roots: complex([0.0,-1.0,0.0,0.0,0.0,1.0])
# pol PCF depth>1: complex([1.0,0.0,-4.5,3.0])
PL.lift_fundamental_group(complex([-1.0,0.0,0.0,1.0]),1)
gcf()

#clf()