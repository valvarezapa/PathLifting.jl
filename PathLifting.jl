module PL # Path Lifting.

using Colors
using PyPlot
using Polynomials
using DelimitedFiles
using PyCall
using Plots

np=pyimport("numpy")


# Computes the norm of the maximum of a complex number.
function norm_max(z::Complex{Float64})::Float64
    # Pre:
    # Post: Returns the maximum between the real part and the imaginary part of the given complex number.
    return max(abs(real(z)),abs(imag(z)))
  end
  
  
  #Example
  #norm_max(1.0-3.0*im)


# Computes the distance (considering the norm of the maximum) between two complex numbers.
function dist_max(z1::Complex{Float64},z2::Complex{Float64})::Float64
  # Pre:
  # Post: Returns the distance between two given complex numbers.
  return norm_max(z2-z1)
end


#Example
#dist_max(1.0-3.0*im,complex(25.0))


# Computes the distance (considering the norm of the maximum) between two complex numbers.
function dist_max(points::Tuple{Complex{Float64},Complex{Float64}})::Float64
  # Pre:
  # Post: Returns the distance between two given complex numbers.
  return norm_max(points[2]-points[1])
end


#Example
#dist_max((1.0-3.0*im,complex(25.0)))


# Computes the distance (considering the norm of the maximum) between a complex number and a list of points.
function dist_max(z::Complex{Float64},points::Array{Complex{Float64},1})::Float64
  # Pre:
  # Post: Returns the distance between the given point and the list of complex numbers.
  l=length(points)
  if l>0
    dis=dist_max(z,points[1])
    @inbounds for i in 2:l
      if dist_max(z,points[i])<dis
        dis=dist_max(z,points[i])
      end
    end
  else
    dis=0
  end
  return dis
end


#Example
#dist_max((1.0-3.0*im,complex(25.0)))


# Computes the minimum distance (considering the norm of the maximum) in a list of complex numbers.
function min_dist_max(points::Array{Complex{Float64},1})::Float64
  # Pre:
  # Post: Returns the minimum distance between points of the given list.
  l=length(points)
  if l>1
    min=dist_max(points[1],points[2])
    @inbounds for i in 1:l
      @inbounds for j in 1:l
        if i!=j && dist_max(points[i],points[j])<min
          min=dist_max(points[i],points[j])
        end
      end
    end
  else
    min=0
  end
  return min
end


#Example
#min_dist_max((1.0-3.0*im,complex(25.0),3.05+18*im))


# Computes the direction of a discrete path joining two complex numbers (z1 and z2), in two steps.
# The first component of the resulting vector is the direction which z1 has to follow in order to get to the intermediate point m1.
# The second component of the resulting vector is the direction which m1 has to follow in order to get to the next point m2 of the discrete path.
# In case the two points are vertically or horizontally aligned, the resulting vector only has one element.
function direction(z1::Complex{Float64},z2::Complex{Float64})::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the direction vector.
  r1=real(z1)
  i1=imag(z1)
  r2=real(z2)
  i2=imag(z2)
  if abs(r2-r1)<abs(i2-i1) 
    if r2!=r1 
      return [((i2-i1)/abs(i2-i1))*im,complex((r2-r1)/abs(r2-r1))]
    else 
      return [((i2-i1)/abs(i2-i1))*im]
    end
  else
    if i2!=i1
      return [complex((r2-r1)/abs(r2-r1)),((i2-i1)/abs(i2-i1))*im]
    else
      return [complex((r2-r1)/abs(r2-r1))]
    end
  end
end


#Example
#direction(1.0-3.0*im,complex(25.0))


# Constructs a discrete corner joining two complex numbers.
function corner(z1::Complex{Float64},z2::Complex{Float64})::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the list of points that make up the corner.
  r1=real(z1)
  i1=imag(z1)
  r2=real(z2)
  i2=imag(z2)
  if abs(r2-r1)<abs(i2-i1)
    if r2!=r1 
      return [z1,r1+i2*im,z2]
    else 
      return [z1,z2]
    end
  else
    if i2!=i1
      return [z1,r2+i1*im,z2]
    else 
      return [z1,z2]
    end
  end
end


#Example
#corner(1.0-3.0*im,complex(25.0))


# Constructs a discrete path joining the two given complex numbers z1 and z2, avoiding to intersect a neighborhood of each point in the given list.
function polygonal_path(points::Array{Complex{Float64},1},z1::Complex{Float64},z2::Complex{Float64})::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the list of points that make up the polygonal path.
  if z1==z2
    pol_path=vcat(polygonal_path(points,z1,z1+1),polygonal_path(points,z1+1,z1))
    return pol_path
  else
    l=length(points)
    if l==0
      return corner(z1,z2)
    else
      if l==1
        mult=minimum([dist_max(z1,points[1]),dist_max(z2,points[1])])/3.0
      else
        mult=minimum([dist_max(z1,points),dist_max(z2,points),min_dist_max(points)])/3.0 # REVISAR FORMULA
      end
      if dist_max(z1,z2)<=mult
        return corner(z1,z2)
      else
        pol_path=[z1]
        dir=direction(last(pol_path),z2)
        while (dist_max(last(pol_path),z2)>mult)
          if (dist_max(last(pol_path),points)>=mult)
            dir=direction(last(pol_path),z2)
            @inbounds for i in 1:length(dir)
              push!(pol_path,last(pol_path)+mult*dir[i])
            end
          else
            pop!(pol_path)
            if length(dir)==1
              dir=[dir[1]*im]
            else
              dir=[dir[1]]
            end
            push!(pol_path,last(pol_path)+1.1*mult*dir[1])
          end
        end
        cor=corner(last(pol_path),z2)
        len_cor=length(cor)
        @inbounds for i in 1:len_cor
          push!(pol_path,cor[i])
        end
        return pol_path
      end
    end
  end
end


#Example
#polygonal_path([complex(10.0),15.0+3*im],1.0-3.0*im,complex(25.0))


# Plots a discrete path joining the two given complex numbers, avoiding to intersect a neighborhood of each point in the given list.
function plot_polygonal_path(points::Array{Complex{Float64},1},z1::Complex{Float64},z2::Complex{Float64};plotcolor::AbstractString="blue")
  # Pre:
  # Post: Plots the polygonal path.
  len_points=length(points)
  pol_path=polygonal_path(points,z1,z2)
  len_path=length(pol_path)
  if len_points==0
    @inbounds for i in 1:len_path-1
      x1=real(pol_path[i])
      y1=imag(pol_path[i])
      x2=real(pol_path[i+1])
      y2=imag(pol_path[i+1])
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
    end
  else
    if len_points==1
      mult=minimum([1.0,dist_max(z1,points[1]),dist_max(z2,points[1])])
    else
      mult=minimum([dist_max(z1,points),dist_max(z2,points),min_dist_max(points)/3.0])
    end
    @inbounds squares=[((real(points[k])-mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])+mult),(real(points[k])-mult,imag(points[k])+mult)) for k=1:len_points]
    @inbounds for i in 1:len_points
      @inbounds for j in 1:4
        if j==4
          x1=squares[i][j][1]
          y1=squares[i][j][2]
          x2=squares[i][1][1]
          y2=squares[i][1][2]
        else
          x1=squares[i][j][1]
          y1=squares[i][j][2]
          x2=squares[i][j+1][1]
          y2=squares[i][j+1][2]
        end
        x_values=[x1,x2]
        y_values=[y1,y2]
        PyPlot.plot(x_values,y_values,color="gray",linestyle="-")
      end
    end
    @inbounds for i in 1:len_path-1
      x1=real(pol_path[i])
      y1=imag(pol_path[i])
      x2=real(pol_path[i+1])
      y2=imag(pol_path[i+1])
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
    end
  end
end


#Example
#clf() # Can be useful to clear the current plot.

#plot_polygonal_path([complex(10.0)],1.0-3.0*im,complex(25.0))
#gcf()


# Plots a given discrete path.
function plot_polygonal_path(pol_path::Array{Complex{Float64},1};plotcolor::AbstractString="blue")
  # Pre:
  # Post: Plots the polygonal path.
  len_path=length(pol_path)
  @inbounds for i in 1:len_path-1
    x1=real(pol_path[i])
    y1=imag(pol_path[i])
    x2=real(pol_path[i+1])
    y2=imag(pol_path[i+1])
    x_values=[x1,x2]
    y_values=[y1,y2]
    PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
  end
end


# Constructs a discrete segment joining two complex numbers.
# If the value of the variable 'string' is 'c', the whole segment will be returned as an array of complex points.
# If 'string' is 'r', the last point of the segment (z2) will be removed, whereas if 'string' is 'l', the first
# point (z1) will be removed.
# If the value of the variable 'reverse' is -1, the program returns the corresponding reverse segment between z2 and z1.
function discrete_segment(z1::Complex{Float64},z2::Complex{Float64},num_points::Int64=20;string::Char='c',reverse::Int64=1)::Array{Complex{Float64},1} # ¿PUEDE REDEFINIRSE POLYGONAL PATH CON MULTIPLE DISPATCH PARA EL CASO EN QUE NO HAYA QUE ESQUIVAR NADA?
  # Pre:
  # Post: Returns the list of points that make up the segment.
  if reverse==1
    @inbounds seg=[z1+(i/num_points)*(z2-z1) for i in 0:num_points]
    if string=='c'
      return seg
    elseif string=='r'
      pop!(seg)
      return seg
    elseif string=='l'
      popfirst!(seg)
      return seg
    end
  elseif reverse==-1
    if string=='c'
      return discrete_segment(z2,z1,num_points;string='c',reverse=1)
    elseif string=='r' 
      return discrete_segment(z2,z1,num_points;string='l',reverse=1)
    elseif string=='l'
      return discrete_segment(z2,z1,num_points;string='r',reverse=1)
    end
  end
end


#Example
#discrete_segment(1.0-3.0*im,complex(25.0),50;string='r',reverse=1)


# Constructs a discrete square loop of given 'radius' around a complex number.
# If the value of the variable 'string' is 'c', the whole loop will be returned as an array of complex points.
# If 'string' is 'r', the last point of the loop (center) will be removed, whereas if 'string' is 'l', the first
# point (center) will be removed.
# If the value of the variable 'reverse' is -1, the program returns the corresponding reversed loop.
function discrete_loop(center::Complex{Float64},radius::Float64,num_points::Int64=20;string::Char='c',reverse::Int64=1)::Array{Complex{Float64},1} # ¿PUEDE REDEFINIRSE POLYGONAL PATH CON MULTIPLE DISPATCH PARA EL CASO EN QUE NO HAYA QUE ESQUIVAR NADA?
  # Pre:
  # Post: Returns the list of points that make up the loop.
  first_half_right=discrete_segment(center+complex(radius), center+complex(radius)+radius*im, num_points;string='r',reverse=1)
  up=discrete_segment(center+complex(radius)+radius*im, center-complex(radius)+radius*im, 2*num_points;string='r',reverse=1)
  left=discrete_segment(center-complex(radius)+radius*im, center-complex(radius)-radius*im, 2*num_points;string='r',reverse=1)
  down=discrete_segment(center-complex(radius)-radius*im, center+complex(radius)-radius*im, 2*num_points;string='r',reverse=1)
  second_half_right=discrete_segment(center+complex(radius)-radius*im, center+complex(radius), num_points;string='c',reverse=1)
  loop=[first_half_right; up; left; down; second_half_right]
  if reverse==1
    if string=='c'
      return loop
    elseif string=='r'
      pop!(loop)
      return loop
    elseif string=='l'
      popfirst!(loop)
      return loop
    elseif string=='o'
      popfirst!(loop)
      pop!(loop)
      return loop
    end
  elseif reverse==-1
    if string=='c'
      return reverse!(loop)
    elseif string=='r' 
      pop!(loop)
      return reverse!(loop)
    elseif string=='l'
      popfirst!(loop)
      return reverse!(loop)
    elseif string=='o'
      popfirst!(loop)
      pop!(loop)
      return reverse!(loop)
    end
  end
end


#Example
#discrete_loop(1.0-3.0*im,2.5,50;string='r',reverse=1)


# Computes the base point of the given list of points.
# The base point is considered to be the real point with value equal to the maximum real part of the points of the list plus the minimum distance (using the norm of the maximum) between the points.
function basepoint(points::Array{Complex{Float64},1})::Complex{Float64}
  # Pre:
  # Post: Returns the base point.
  l=length(points)
  if l==1
    return complex(real(points[1])+1.0/3.0)
  elseif l==0
    return complex(0.0)
  else
    @inbounds points_real=[real(points[i]) for i in 1:l]
    xmax=maximum(points_real)
    xmax1=xmax+min_dist_max(points)
    return complex(xmax1)
  end
end


#Example
#basepoint([1.0-3.0*im,complex(2.5),50.0+6.2*im])


# Computes a discrete loop that starts and ends at the base point of the given list of points, and that loops around the 'loop_center' 
# avoiding to intersect a neighborhood of radius 'radius' of the points in the given list.
function discrete_loop(points::Array{Complex{Float64},1},loop_center::Complex{Float64},radius::Float64,num_points::Int64=20)::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the discrete loop starting and ending at 'loop_center'.
  bp=basepoint(points)
  segment=polygonal_path(points, bp, loop_center+complex(radius))
  loop_around=discrete_loop(loop_center,radius,num_points;string='o',reverse=1)
  return [segment; loop_around; reverse(segment)]
end


#Example
#discrete_loop([1.0-3.0*im,complex(2.5),50.0+6.2*im],complex(0.0),0.5)


# Computes a list of discrete loops that start and end at the base point of the given list of points, and that loops around the each point of the list, 
# avoiding to intersect a neighborhood of the  rest of the points in the given list.
function discrete_loops(points::Array{Complex{Float64},1},num_points::Int64=20)::Array{Array{Complex{Float64},1},1}
  # Pre:
  # Post: Returns the list of discrete loops.
  loops=[]
  l=length(points)
  if l==1
    radius=1.0/3.0
  elseif l==0
    return loops
  else
    radius=min_dist_max(points)/3.0
  end
  @inbounds for i in 1:l
    new_loop=discrete_loop(points,points[i],radius,num_points)
    push!(loops,new_loop)
  end
  return loops
end


#Example
#discrete_loops([1.0-3.0*im,complex(2.5),50.0+6.2*im])


# Plots a given loop, avoiding to intersect a neighborhood of each point in the given list.
function plot_discrete_loop(points::Array{Complex{Float64},1},loop::Array{Complex{Float64},1};plotcolor::AbstractString="blue",enable_plot::Bool=true)::Array{Any,1}
  # Pre: The loop does not intersect a fixed neighborhood of any of the points in the given list.
  # Post: Plots the loop and the obstacles, and returns their elements joined in a list.
  l=length(points)
  len_loop=length(loop)
  if l==0
    if enable_plot==true
      @inbounds for i in 1:len_loop-1
        x1=real(loop[i])
        y1=imag(loop[i])
        x2=real(loop[i+1])
        y2=imag(loop[i+1])
        x_values=[x1,x2]
        y_values=[y1,y2]
        PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
      end
    end
    return loop
  else
    mult=min_dist_max(points)/3.0
    @inbounds squares=[((real(points[k])-mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])+mult),(real(points[k])-mult,imag(points[k])+mult)) for k=1:l]
    if enable_plot==true
      @inbounds for i in 1:l
        @inbounds for j in 1:4
          if j==4
            x1=squares[i][j][1]
            y1=squares[i][j][2]
            x2=squares[i][1][1]
            y2=squares[i][1][2]
          else
            x1=squares[i][j][1]
            y1=squares[i][j][2]
            x2=squares[i][j+1][1]
            y2=squares[i][j+1][2]
          end
          x_values=[x1,x2]
          y_values=[y1,y2]
          PyPlot.plot(x_values,y_values,color="gray",linestyle="-")
        end
      end
      @inbounds for i in 1:len_loop-1
        x1=real(loop[i])
        y1=imag(loop[i])
        x2=real(loop[i+1])
        y2=imag(loop[i+1])
        x_values=[x1,x2]
        y_values=[y1,y2]
        PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
      end
    end
    return [loop; squares]
  end
end


#Example
#plot_discrete_loop([1.0-3.0*im,complex(2.5),50.0+6.2*im],loop)


# Normalizes a given polynomial, assuming that the norm of the polynomial
# is defined as the sum in absolute value of its coefficients.
function normalize_pol(coeffs::Array{Complex{Float64},1})::Array{Complex{Float64},1}
  # Pre:
  # Post: Retuns the normalized polynomial.
  norm=sum(abs,coeffs)
  return coeffs*(1.0/norm)
end


#Example
#=
pol=complex([1.0,0.0,0.0,2.0]) # f(z)=z^3-1
norm_pol=normalize_pol(pol)
=#


# Computes the next point.
function new_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64})::Complex{Float64}
  # Pre: The point 'x' is not a critical point of the given polynomial.
  # Post: Computes the next point.
  p=Polynomial(coeffs)
  der_p=derivative(p)
  @assert der_p(x)!=0
  return ((new_y-y)/der_p(x))+x
end


#Example
#x_new(Polynomial([1.0-3.0*im,complex(2.5),50.0+6.2*im]),y,new_y,x)

#=

# Computes the next point (for rational maps).
function new_x(pair::Tuple{Array{Complex{Float64},1},Array{Complex{Float64},1}},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64})::Complex{Float64}
  # Pre: The point 'x' is not a critical point of the rational map given by its numerator and denominator 'pair'.
  # Post: Computes the next point.
  p=Polynomial(pair[1])
  q=Polynomial(pair[2])
  der_p=derivative(p)
  der_q=derivative(q)
  num_der=der_p*q-p*der_q
  den_der=q*q
  @assert num_der(x)!=0
  return complex((new_y-y)*den_der(x)/num_der(x)+x)
end


#Example
#x_new((Polynomial([1.0-3.0*im,complex(2.5),50.0+6.2*im]),Polynomial([1.0,im])),y,new_y,x)

=#


# Computes the last point.
function new_precise_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64};precision::Float64=8.0)::Complex{Float64}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post: Computes the last point.
  p=Polynomial(coeffs)
  while abs(p(new_x(coeffs,y,new_y,x))-new_y)>1.0/10^precision
    x=new_x(coeffs,y,new_y,x)
    y=p(x)
  end
  return new_x(coeffs,y,new_y,x)
end


#Example
#new_precise_x([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function new_precise_x_list(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64};pre::Float64=8.0)::Array{Complex{Float64},1}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post:
  return [x,new_precise_x(coeffs,y,new_y,x;precision=pre)]
end


#Example
#new_precise_x_list([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function precise_lifting(coeffs::Array{Complex{Float64},1},y_points::Array{Complex{Float64},1},x::Complex{Float64};pre::Float64=8.0)::Complex{Float64}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post:
  p=Polynomial(coeffs)
  l=length(y_points)
  @inbounds y_list=[y_points[i] for i in 1:l]
  @assert l>=2
  if l==2 
    return new_precise_x(coeffs,y_points[1],y_points[2],x;precision=pre)
  else
    y=y_list[1]
    new_y=y_list[2]
    new_x=new_precise_x(coeffs,y,new_y,x;precision=pre)
    while length(y_list)>2
      popfirst!(y_list)
      y=y_list[1]
      new_y=y_list[2]
      new_x=new_precise_x(coeffs,y,new_y,new_x;precision=pre)
    end
    return new_precise_x(coeffs,y,new_y,new_x;precision=pre)
  end
end


#Example
#precise_lifting([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function precise_lifting_list(coeffs::Array{Complex{Float64},1},y_points::Array{Complex{Float64},1},x::Complex{Float64};prec::Float64=8.0)::Array{Complex{Float64},1}
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post:
  p=Polynomial(coeffs)
  l=length(y_points)
  @inbounds y_list=[y_points[i] for i in 1:l]
  @assert l>=2
  if l==2 
    return new_precise_x_list(coeffs,y_points[1],y_points[2],x;pre=prec)
  else
    y=y_list[1]
    new_y=y_list[2]
    new_x_list=new_precise_x_list(coeffs,y,new_y,x;pre=prec)
    while length(y_list)>2
      popfirst!(y_list)
      y=y_list[1]
      new_y=y_list[2]
      new_x=new_precise_x(coeffs,y,new_y,last(new_x_list);precision=prec)
      push!(new_x_list,new_x)
    end
    push!(new_x_list,new_precise_x(coeffs,y,new_y,last(new_x_list);precision=prec))
    return new_x_list
  end
end


#Example
#precise_lifting([1.0-3.0*im,complex(2.5),50.0+6.2*im],y,new_y,x;precision=8)


# 
function precise_lifting_plot(coeffs::Array{Complex{Float64},1},y_points::Array{Complex{Float64},1},x::Complex{Float64};preci::Float64=8.0,plotcolor::AbstractString="blue")
  # Pre: The point 'x' is not a critical point of a the polynomial given by its coefficients 'coeffs'.
  # Post: 
  p=Polynomial(coeffs)
  l=length(y_points)
  coveringlist=precise_lifting_list(coeffs,y_points,x;prec=preci)
  @inbounds for i in 1:l
    x1=real(coveringlist[i])
    y1=imag(coveringlist[i])
    x2=real(coveringlist[i+1])
    y2=imag(coveringlist[i+1])
    x_values=[x1,x2]
    y_values=[y1,y2]
    PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
  end
  @inbounds for i in 1:l
    x1=real(y_points[i])
    y1=imag(y_points[i])
    x2=real(y_points[i+1])
    y2=imag(y_points[i+1])
    x_values=[x1,x2]
    y_values=[y1,y2]
    PyPlot.plot(x_values,y_values,color="red",linestyle="-")
  end
end


#Example
#precise_lifting_plot(Polynomial(complex([1.0,0.0,1.0])),polygonal_path([1.0+1.0*im],complex(0.0),2.0+2.0*im),3.0*im)


# Computes a vector containing the positions in which 'element' appears in the given list 'list'.
function positions(list,element)
  # Pre:
  # Post: Returns the vector of positions.
  positions_vector=[]
  l=length(list)
  @inbounds for i in 1:l
    if list[i]==element
      push!(positions_vector,i)
    end
  end
  return positions_vector
end


#Example
#positions([1,2,4,6,1,25,3,8,4,2,4],4)


# 
function new_taylor_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64})
  # Pre: The polynomial given by 'coeffs' is not zero.
  # Post: 
  degree=length(coeffs)-1
  p=Polynomial(coeffs)
  derivatives_p=[p]
  @inbounds for k in 1:degree
    push!(derivatives_p,derivative(last(derivatives_p)))
  end
  popfirst!(derivatives_p)
  @inbounds taylor_coeffs=[derivatives_p[k](x)/factorial(k) for k in 1:degree]
  @inbounds potential=[abs(taylor_coeffs[k])^(1.0/k) for k in 1:degree]
  @inbounds potenitialy=[abs(new_y-y)^(1.0/k) for k in 1:degree]
  @inbounds comparation=[potential[k]>potenitialy[k] for k in 1:degree]
  first=positions(comparation, true)[1]
  @assert abs(taylor_coeffs[first])!=0.0
  return ((new_y-y)/taylor_coeffs[first])^(1.0/first)+x
end


#Example
#new_taylor_x()


# 
function new_precise_taylor_x(coeffs::Array{Complex{Float64},1},y::Complex{Float64},new_y::Complex{Float64},x::Complex{Float64};precision::Float64=8.0)
  # Pre: The polynomial given by 'coeffs' is not zero.
  # Post: 
  p=Polynomial(coeffs)
  while abs(p(new_taylor_x(coeffs,y,new_y,x))-new_y)>1.0/10^precision
    x=new_taylor_x(coeffs,y,new_y,x)
    y=p(x)
  end
  return new_taylor_x(coeffs,y,new_y,x)
end


#Example
#new_precise_taylor_x()


# 
function base_change(path::Array{Complex{Float64},1},loops::Array{Array{Complex{Float64},1},1})::Array{Array{Complex{Float64},1},1}
  # Pre:
  # Post: 
  len_loops=length(loops)
  len_path=length(path)
  @inbounds path_aux=[path[k] for k in 1:len_path]
  pop!(path_aux)
  new_loops=[]
  @inbounds for i in 1:len_loops
    push!(new_loops,[path_aux; loops[i]; reverse(path_aux)])
  end
  return new_loops
end


#Example
#base_change()


# 
function new_solution_list(coeffs::Array{Complex{Float64},1},initial_points::Array{Complex{Float64},1},generator::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  p=Polynomial(coeffs)
  der_p=derivative(p)
  l=length(initial_points)
  @inbounds solution_list=[initial_points[i] for i in 1:l]
  len_sol=length(solution_list)
  @inbounds for k in 1:len_sol
    next_sol=precise_lifting(coeffs,generator,solution_list[k];pre=precision)
    @inbounds control=[abs(der_p(next_sol))*dist_max(next_sol,solution_list[i]) for i in 1:len_sol]
    if minimum(control)>=1.0/10.0^precision # If the new solution is neither a critical point nor an already computed solution on the list:
      push!(solution_list,next_sol)
      len_sol=length(solution_list)
    end
  end
  return solution_list
end


#Example
#new_solution_list()


# 
function new_solution_list(coeffs::Array{Complex{Float64},1},initial_points::Array{Complex{Float64},1},generators::Array{Array{Complex{Float64},1},1};pre::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  len_gens=length(generators)
  l=length(initial_points)
  @inbounds solution_list=[initial_points[k] for k in 1:l]
  i=1
  while i<=len_gens
    solution_list=new_solution_list(coeffs,solution_list,generators[i];precision=pre)
    i=i+1
  end
  return solution_list
end


#Example
#new_solution_list()


# 
function generators_initial_point(coeffs::Array{Complex{Float64},1},initial_point::Complex{Float64},critical_values::Array{Complex{Float64},1};pre::Float64=8.0)::Array{Array{Complex{Float64},1},1}
  # Pre:
  # Post:
  p=Polynomial(coeffs)
  generators=discrete_loops(critical_values)
  bp=generators[1][1]
  new_bp=p(initial_point)
  path=polygonal_path(critical_values,new_bp,bp)
  new_gens=base_change(path,generators)
  return new_gens
end


#Example
#generators_initial_point()


# 
function new_nsolve_aux(coeffs::Array{Complex{Float64},1},initial_point::Complex{Float64},critical_values::Array{Complex{Float64},1};preci::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  deg=length(coeffs)-1
  new_gens=generators_initial_point(coeffs,initial_point,critical_values;pre=preci)
  solution_list=[initial_point]
  len_sol=length(solution_list)
  solution_list=new_solution_list(coeffs,solution_list,new_gens;pre=preci)
  mod_len_sol=length(solution_list)
  while len_sol<mod_len_sol && mod_len_sol!=deg
    len_sol=length(solution_list)
    solution_list=new_solution_list(coeffs,solution_list,new_gens;pre=preci)
    mod_len_sol=length(solution_list)
  end    
  return solution_list
end


#Example
#new_nsolve_aux()


# 
function search_initial_point(coeffs::Array{Complex{Float64},1},value_list::Array{Complex{Float64},1},distance::Float64)::Complex{Float64}
  # Pre:
  # Post:
  initial_point=complex(0.0)
  p=Polynomial(coeffs)
  next_point=p(initial_point)
  while dist_max(next_point,value_list)<distance
    initial_point=initial_point+complex(1.0)
    next_point=p(initial_point)
  end
  return initial_point
end


#Example
#search_initial_point()


# 
function non_quasicritical_value(value::Complex{Float64},critical_values::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Any,1}
  # Pre:
  # Post:
  l=length(critical_values)
  @inbounds distances=[dist_max(value,critical_values[k]) for k in 1:l]
  if l==0
    return [true] # 'value' is not among the given critical values.
  else
    min=minimum(distances)
    if min<1.0/10.0^precision
      positions_list=positions(distances,min)
      critical_values_list=[false,positions_list[1],critical_values[positions_list[1]]] # 'value' is the critical value in the position positions_list[1].
      return critical_values_list
    else
      return [true] # 'value' is not among the given critical values.
    end
  end
end


#Example
#non_quasicritical_value()


# CALCULA p^{-1}(noncrit_value)
function noncritical_nsolve_aux(coeffs::Array{Complex{Float64},1},noncrit_value::Complex{Float64},critical_values::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  p=Polynomial(coeffs)
  if length(critical_values)==1
    dist=1.0/3.0
  else
    dist=min_dist_max(critical_values)/3.0
  end
  initial_point=search_initial_point(coeffs,critical_values,dist)
  path=polygonal_path(critical_values,p(initial_point),noncrit_value)
  initial_sol=new_nsolve_aux(coeffs,initial_point,critical_values;preci=precision)
  l=length(initial_sol)
  @inbounds final_sol=[precise_lifting(coeffs,path,initial_sol[k];pre=precision)  for k in 1:l]
  return final_sol
end


#Example
#noncritical_nsolve_aux()


# 
function nsolve_aux(coeffs::Array{Complex{Float64},1},value::Complex{Float64},critical_values::Array{Complex{Float64},1};precis::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  res=non_quasicritical_value(value,critical_values;precision=precis)
  if res[1]==true # In case the value is not among the critical values given on the list.
    sol=noncritical_nsolve_aux(coeffs,value,critical_values;precision=precis)
  else # In case the value is one of the critical values given on the list.
    if length(critical_values)==1
      min1=1.0/3.0
    else
      min1=min_dist_max(critical_values)/3.0
    end
    new_value=non_quasicritical_value(value,critical_values;precision=precis)[3]+complex(min1)
    provisionalsol=noncritical_nsolve_aux(coeffs,new_value,critical_values;precision=precis)
    l=length(provisionalsol)
    @inbounds sol=[new_precise_taylor_x(coeffs,new_value,value,provisionalsol[k];precision=precis) for k in 1:l]
  end
  return sol
end


#Example
#nsolve_aux()


#=

# 
function delete_repited_aux(list::Array{Any,1})::Array{Any,1}
  # Pre:
  # Post:
  sol_list=list[1]
  point_list=list[2]
  prec=list[3]
  len_points=length(point_list)
  @inbounds final_point_list=[point_list[k] for k in 2:len_points]
  if dist_max(point_list[1],final_point_list)>1.0/(10.0^prec)
    sol_list=[sol_list; point_list[1]]
  end
  return [sol_list,final_point_list,prec]
end


#Example
#delete_repited_aux()


# 
function delete_repited(points::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  l=length(points)
  @inbounds mod_points=[points[k] for k in 1:l]
  len_mod_points=length(mod_points)
  list=[[],points,precision]
  while len_mod_points>0
    list=delete_repited_aux(list)
    len_mod_points=length(list[2])
  end
  return list[1]
end


#Example
#delete_repited()

=#


# 
function is_in(el::Complex{Float64},list::Array{Complex{Float64},1};pre::Float64=8.0)::Bool
  # Pre:
  # Post:
  l=length(list)
  for i in 1:l
    if abs(el-list[i])<1.0/(10.0^pre)
      return true
    end
  end
  return false
end


#Example
#is_in()


# 
function delete_repited(list::Array{ComplexF64,1};precision::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  l=length(list)
  new_list=[first(list)]
  for i in 2:l
    if is_in(list[i],new_list;pre=precision)==false
      push!(new_list,list[i])
    end
  end
  return new_list
end


#Example
#delete_repited()


# AVOIDS UNDERFLOWS
function clean_list(list::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  l=length(list)
  new_list=[]
  for i in 1:l
    x=real(list[i])
    y=imag(list[i])
    if abs(x)<1.0/(10.0^precision)
      x=0
    end
    if abs(y)<1.0/(10.0^precision)
      y=0
    end
    push!(new_list,x+y*im)
  end
  return new_list
end


#Example
#clean_list()


# 
function nsolve(coefs::Array{Complex{Float64},1},value::Complex{Float64};pre::Float64=8.0)::Array{Any,1}
  # Pre: The polynomial is not constant.
  # Post:
  coefs[1]=coefs[1]-value
  coefs=normalize_pol(coefs)
  deg=length(coefs)-1
  @assert deg>0
  if deg==1
    return [[0,(-coefs[1])/coefs[2]]]
  else
    p=Polynomial(coefs)
    p_mod=p
    der_p_list=[coefs]
    println("Polynomial and its derivatives:")
    println(coefs)
    @inbounds for k in 1:deg-1
      p_mod=derivative(p_mod)
      print(k);print("-th derivative: ")
      println(p_mod)
      push!(der_p_list,normalize_pol(coeffs(p_mod)))
    end
    reverse!(der_p_list)
    coef=der_p_list[1]
    criticalpointlastpol=[-coef[1]/coef[2]]
    criticalvaluelist=[Polynomial(der_p_list[2])(criticalpointlastpol[1])]
    println("Critical points and values:")
    println("1- ")
    print("C.points: ");println(criticalpointlastpol)
    print("C.values: ");println(criticalvaluelist)
    counter=2
    sol=[[deg-counter+1,criticalpointlastpol]]
    value1=complex(0.0)
    while counter<deg
      criticalpointlastpol=nsolve_aux(der_p_list[counter],value1,criticalvaluelist;precis=pre)
      criticalpointlastpol=clean_list(criticalpointlastpol;precision=pre)
      criticalpointlastpol=delete_repited(criticalpointlastpol;precision=pre)
      print(counter);println("- ")
      print("C.points: ");println(criticalpointlastpol)
      counter=counter+1
      sol=[sol;[[deg-counter+1,criticalpointlastpol]]]
      l=length(criticalpointlastpol)
      criticalvaluelist=[Polynomial(der_p_list[counter])(criticalpointlastpol[i]) for i in 1:l]
      criticalvaluelist=clean_list(criticalvaluelist;precision=pre)
      criticalvaluelist=delete_repited(criticalvaluelist;precision=pre)
      print("C.values: ");println(criticalvaluelist)
    end
    println("Last step:")
    print("Polynomial ");print(der_p_list[counter]);print(" in value ");print(value1);print(" avoiding ");println(criticalvaluelist)
    criticalpointlastpol=nsolve_aux(der_p_list[counter],value1,criticalvaluelist;precis=pre)
    criticalpointlastpol=clean_list(criticalpointlastpol;precision=pre)
    criticalpointlastpol=delete_repited(criticalpointlastpol;precision=pre)
    println("Solutions: ");println(criticalpointlastpol)
    return [sol;[[0,criticalpointlastpol]]]
  end
end


#Example
#nsolve()


# 
function zeroes_mult(coefs::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Any,1}
  # Pre:
  # Post:
  sol=nsolve(coefs,complex(0.0);pre=precision)
  l=length(sol)
  reverse!(sol)
  zeroes=first(sol)[2]
  num_zeroes=length(zeroes)
  mult_zero_pair=[]
  for i in 1:num_zeroes
    mult=1
    cond=true
    j=2
    while j<=l && cond==true
      if is_in(zeroes[i],sol[j][2];pre=precision)==true
        mult=mult+1
      else
        cond=false
      end
      j=j+1
    end
    push!(mult_zero_pair,(zeroes[i],mult))
  end
  return mult_zero_pair
end


#Example
#zeroes_mult()



# ---------------- LITTLEWOOD ----------------#



# 
function binary(n::Int64)::Array{Bool,1}
  # Pre:
  # Post:
  if n==0
    return [false]
  else
    list=[]
    while n>1
      r=n%2
      n=(n-r)/2
      if r==1
        push!(list,true)
      else
        push!(list,false)
      end
    end
    push!(list,true)
    return list
  end
end


#Example
#binary()


# Ensures that the two given lists have the same length, adding complex zeros if necessary.
function sameLength(list1,list2)
  # Pre: 
  # Post: The lists have been modified if necessary to have the same length.
  l1=length(list1)
  l2=length(list2)
  if l1!=l2
      if l1<l2
       while l1<l2
           push!(list1,complex(0.0))
             l1=l1+1
          end
      else
          while l2<l1
             push!(list2,complex(0.0))
             l2=l2+1
          end
      end
  end
  le=l1 # Note that at this point, it is ensured that l1=l2.
  while le>1 && list1[le]==complex(0.0) && list2[le]==complex(0.0)
      pop!(list1)
      pop!(list2)
      le=le-1
  end
end


#Example
#=
list1=[1.0,2.0]
list2=[1.0,2.0,3.0,4.0,5.0,0.0]
sameLength(list1,list2)
=#


# Computes the roots of every Littlewood polynomial of degree deg>0 and writes them in a text file.
function littlewood_roots(deg::Int64;pre::Float64=8.0)
  # Pre:
  # Post:
  @assert deg>0
  touch("roots.txt")
  @inbounds init_pol=[complex(1.0) for i in 1:deg+1]
  zeroes=[]
  @inbounds for i in 0:2^(deg+1)-1
    mod=binary(i)
    sameLength(mod,init_pol)
    pol::Array{Complex{Float64},1}=[]
    @inbounds for i in 1:deg+1
      if mod[i]==true
        push!(pol,complex(-1.0))
      else
        push!(pol,complex(1.0))
      end
    end
    print(i);print("- Polynomial ");print(pol);print(": ")
    zeroes_pol=zeroes_mult(pol;precision=pre)
    println("done!")
    len_zeroes=length(zeroes_pol)
    @inbounds for i in 1:len_zeroes
      mult=zeroes_pol[i][2]
      @inbounds for j in 1:mult
        push!(zeroes,zeroes_pol[i][1])
      end
    end
  end
  open("roots.txt","w") do io
    writedlm(io,zeroes)
  end
end


#Example
#littlewood_roots()


# Builds a collection of points of a rectangle determined by the given intervals and precision,
# respecting the aspect ratio between the given intervals. The given precision is used to
# determine the number of points in the grid.
# Note that the collection contains P^1+(C) points [z:t] that, considering the P^1+(C) to C*
# (that is the Alexandroff compactification of C) bijection, correspond to z in C*.
function rectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
  precision=6)
  # Pre:
  # Post: Returns the collection of points.
  resolution=9.443329*(10.0^(precision)) # We take this 9-million-pixels-resolution as a reference.
  a=xinterval[1]
  b=xinterval[2]
  c=yinterval[1]
  d=yinterval[2]
  p1=floor(sqrt(resolution*((b-a)/(d-c)))) # Note that resolution=p1*p2
  p2=floor(p1*(d-c)/(b-a))
  collection=[complex(r,i) for i=d:-(d-c)/p2:c, r=a:(b-a)/p1:b]
  return collection
end


# Plots the Littlewood fractal.
function plot_littlewood(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
  precision=6)
  # Pre:
  # Post: 
  a=readdlm("roots1.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  pol=[]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots1.txt")
  a=readdlm("roots2.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots2.txt")
  a=readdlm("roots3.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots3.txt")
  a=readdlm("roots4.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots4.txt")
  a=readdlm("roots5.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots5.txt")
  a=readdlm("roots6.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots6.txt")
  a=readdlm("roots7.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots7.txt")
  a=readdlm("roots8.txt",'\t',Complex{Float64},'\n')
  l=size(a)[1]
  for i in 1:l
    push!(pol,a[i][1])
  end
  close("roots8.txt")
  pol=delete_repited(pol)
  l=length(pol)

  x=[]
  for i in 1:l
    push!(x,real(pol[i]))
  end
  y=[]
  for i in 1:l
    push!(y,imag(pol[i]))
  end
  PyPlot.scatter(x,y)

  #=
  rect=rectangle(xinterval,yinterval,precision)
  a=size(rect)[1]
  b=size(rect)[2]
  colorMatrix=Array{Float64}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      colorMatrix[i,j]=0.0
    end
  end
  for i in 1:l
    point_x = a*(real(pol[i])-xinterval[1])/(xinterval[2]-xinterval[1])   
    point_y = b*(imag(pol[i])-yinterval[1])/(yinterval[2]-yinterval[1])  
    if point_x>=0 && point_x<=a-1 && point_y>=0 && point_y<=b-1
      colorMatrix[floor(Int,point_x),floor(Int,point_y)] += 1
    end
  end
  log_max=log(maximum(colorMatrix))
  ind::Array{ComplexF64,1}=[]
  @inbounds for j in 1:b
    @inbounds for i in 1:a
      if colorMatrix[i,j]>0.0
        push!(ind,colorMatrix[i,j])
      end
    end
  end
  @inbounds for j in 1:b
    @inbounds for i in 1:a
      if is_in(complex(colorMatrix[i,j]),ind)
        colorMatrix[i,j]=log(colorMatrix[i,j])/log_max
      end
    end
  end
  fig = plt.figure()
  ax = fig.add_axes([0, 0, 1, 1], aspect=1)
  ax.axis("off")
  ax.imshow(colorMatrix, cmap="afmhot")
  =#
end


# Computes a numerical approximation to the Julia set of the given polynomial.
function julia_set(coefs::Array{Complex{Float64},1},points::Array{Complex{Float64},1},num_preiter::Int64=8;precision::Float64=8.0)
  # Pre:
  # Post:
  @assert length(coefs)>0
  new_points=[]
  @inbounds for i in 1:num_preiter
    l=length(points)
    @inbounds for j in 1:l
      sol=nsolve(coefs,points[j];pre=precision)
      reverse!(sol)
      fiber=first(sol)[2]
      new_points=[new_points;fiber]
    end
    points=new_points
    points=delete_repited(points)
    new_points=[]
  end
  return points
end


#Example
#julia_set()


# Computes and plots the k-th iterated lifting of some representatives of the fundamental group induced by the given polynomial.
function lift_fundamental_group(coefs::Array{Complex{Float64},1},k::Int64=8;pre::Float64=8.0)
  # Pre:
  # Post:
  deg=length(coefs)-1
  @assert deg>0
  p=Polynomial(coefs)
  p_mod=p
  der_p_list=[coefs]
  @inbounds for k in 1:deg-1
    p_mod=derivative(p_mod)
    push!(der_p_list,coeffs(p_mod))
  end
  reverse!(der_p_list)
  coef=der_p_list[1]
  criticalpointlastpol=[-coef[1]/coef[2]]
  criticalvaluelist=[Polynomial(der_p_list[2])(criticalpointlastpol[1])]
  counter=2
  while counter<deg
    criticalpointlastpol=nsolve_aux(der_p_list[counter],complex(0.0),criticalvaluelist;precis=pre)
    criticalpointlastpol=clean_list(criticalpointlastpol;precision=pre)
    criticalpointlastpol=delete_repited(criticalpointlastpol;precision=pre)
    counter=counter+1
    l=length(criticalpointlastpol)
    criticalvaluelist=[Polynomial(der_p_list[counter])(criticalpointlastpol[i]) for i in 1:l]
    criticalvaluelist=clean_list(criticalvaluelist;precision=pre)
    criticalvaluelist=delete_repited(criticalvaluelist;precision=pre)
  end
  print("Critical points of the polynomial: ");println(criticalpointlastpol)
  criticalpointlastpol=nsolve_aux(der_p_list[counter],complex(0.0),criticalvaluelist;precis=pre)
  criticalpointlastpol=clean_list(criticalpointlastpol;precision=pre)
  criticalpointlastpol=delete_repited(criticalpointlastpol;precision=pre)
  print("Roots of the polynomial: ");println(criticalpointlastpol)
  print("Critical values: ");println(criticalvaluelist)
  if length(criticalvaluelist)==1
    dist=1.0/3.0
  else
    dist=min_dist_max(criticalvaluelist)/3.0
  end
  #initial_point=search_initial_point(coefs,criticalvaluelist,dist)
  #print("Initial point of the generators: ");println(initial_point)
  @assert !(is_in(complex(0.0),criticalvaluelist)) # TRATAR CASO A PARTE
  gens=generators_initial_point(coefs,criticalpointlastpol[1],criticalvaluelist)
  len_gens=length(gens)
  #=
  for i in 1:len_gens
    print("Gen ");print(i);print(" starts at ");print(first(gens[i]));print(" and ends at ");println(last(gens[i]))
  end
  =#
  colors=["brown","blue","red","green","yellow"]
  num_colors=length(colors)
  lifted_points=criticalpointlastpol
  len_sol=length(lifted_points)
  for i in 1:k
    print("Depth of lifting ");print(i);println(":")
    new_lifted_points=[]
    len_sol=length(lifted_points)
    new_gens=[]
    len_gens=length(gens)
    for n in 1:len_gens
      print("  Action of Gen: ");print(n);print(" corresponding to critical value: ");println(criticalvaluelist[n])
      for m in 1:len_sol
        lift_gen=precise_lifting_list(coefs,gens[n],lifted_points[m])
        push!(new_lifted_points,last(lift_gen))
        if abs(first(lift_gen)-last(lift_gen))<1.0/10.0^pre
          push!(new_gens,lift_gen)
        end
        print("    Root ");print(lifted_points[m]);print(" to Root ");println(last(lift_gen))
        plot_polygonal_path(lift_gen;plotcolor=colors[n%num_colors+1])
      end
    end
    lifted_points=new_lifted_points
    #lifted_points=clean_list(lifted_points;precision=pre)
    #lifted_points=delete_repited(lifted_points;precision=pre)
    gens=new_gens
  end
  if k==0
    for n in 1:len_gens
      #print("Gen ");print(n);print(" starts at ");print(first(gens[n]));print(" and ends at ");println(last(gens[n]))
      plot_polygonal_path(gens[n];plotcolor=colors[n%num_colors+1])
    end
  end
end


#Example
#


end