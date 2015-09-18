include("../../Toolbox/Julia_Tools/toolbox.jl")


println(MarkovAR_95(4, .9, .05))

x=[5,7]
y=zeros(3,2)
y[2,:]=x[[1,1]]
println(y)

println("Sergio",x)

n=4
x=[4.0, 1.0, 9.7, 5.0]
y=[0.0, 0.0, 0.0, 0.0]
z=[5, 6, 7,0]

ccall((:__toolbox_MOD_sort, "../../Toolbox/Fortran_Tools/Toolbox.so"), Void, 
      (Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}), &n, x, y, z)

println(y)
println(z)

x1 = 7
a1 = [0]
b1 = [0]

r1 = ccall((:__simplemodule_MOD_foo, "./simplemodule.so"), Int64, 
            (Ptr{Int64},), &x1)

println(r1)
println()

ccall((:__simplemodule_MOD_bar, "./simplemodule.so"), Void, 
      (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}), &x1, a1, b1)
      
println(a1[1])
println(b1[1])
println()

x2 = 7.0
a2 = Cdouble[1.0]
b2 = Cdouble[1.0]

ccall((:__simplemodule_MOD_keg, "./simplemodule.so"), Void, 
      (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), &x2, a2, b2)
      
println(a2[1])
println(b2[1])
println()

x3 = [1.0, 2.0, 3.0]
y3 = [0.0, 0.0, 0.0]
z3 = [0, 0, 0,88,93]
ccall((:__simplemodule_MOD_ruf, "./simplemodule.so"), Void, 
      (Ptr{Float64}, Ptr{Float64}, Ptr{Int64}), x3, y3,z3)
      
println(y3)
println(z3)



# n_z=3
# rho=0.9
# sigma=0.05
# z_grid = [1.0,1.0,1.0] #Array(Float64, n_z)
# P_z    = [1.0 1.0 1.0;1.0 1.0 1.0;1.0 1.0 1.0]#Array(Float64, n_z,n_z)
# ccall((:__simplemodule_MOD_markov_95, "./simplemodule.so"), Void, 
#       (Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), n_z,rho,sigma,z_grid,P_z)



x_min=1.0
x_max=2.0
n_grid=10
c_grid=3.0
x_grid = ccall((:__simplemodule_MOD_grid, "./simplemodule.so"), Void, 
        (Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Float64}), x_min,x_max,n_grid,c_grid)

