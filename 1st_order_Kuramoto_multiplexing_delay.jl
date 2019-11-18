using DifferentialEquations
using DelimitedFiles
using Random
using Printf
using StatsBase


# setting the parameters
const N = 200;
const nearest_neighbour = 61.0;
const inter = 0.5;
const intra = 0.5;
const number_of_solitary_states = 2;




#delPos = [10, 30, 50, 70, 90]
#tau = [8.0, 6.0, 2.0, 7.0, 4.0]

delPos = zeros(Int64, number_of_solitary_states*2)
xx = sample(10:80, number_of_solitary_states, replace=false) # Depending on the number of solitary states 'k' we require, '10:80' in this line of code picks 'k' random position for delay between 10th to 80th node.
for i=1:number_of_solitary_states
    delPos[i] = xx[i]
    delPos[i+number_of_solitary_states] = xx[i]+100
end

rand_tau = rand(number_of_solitary_states)*5
tau = zeros(number_of_solitary_states*2)

for i=1:number_of_solitary_states
    tau[i] = rand_tau[i]
    tau[i+number_of_solitary_states] = rand_tau[i]
end
println(delPos)
println(tau)

omega = ones(N)*1.0;

adj = readdlm("network_float.dat") # Importing the network file
#A = transpose(reshape(adj,N,N))
A = zeros(N, N)

for i=1:N
    for j=1:N
	global A, adj, inter, intra
        if i<=N/2 && j<=N/2
	    A[i, j] = adj[i, j] * intra
	elseif i>N/2 && j>N/2
    	    A[i, j] = adj[i, j] * intra
	elseif i<=N/2 && j>N/2
    	    A[i, j] = adj[i, j] * inter
	elseif i>N/2 && j<=N/2
	    A[i, j] = adj[i, j] * inter
	end
    end
end

net_file = @sprintf("network.dat") # Confirming the network imported by creating another network file
fout = open(net_file, "w")
writedlm(fout, A)
close(fout)

function kuramoto_del(du,u,h,p,t)
    global omega, A, nearest_neighbour, inter
    n = floor(Int64, length(u))
    for i = 1:n
        sum_coup = 0.0
        for j = 1:n
            if i in delPos
		if j in delPos
                    l = findall(l->l==i, delPos)
                    sum_coup += A[i, j]*sin(h(p,t-tau[l[1]])[j] - u[i])
		end
            else
                sum_coup += A[i,j]*sin(u[j] - u[i])
            end
        end
        du[i] = omega[i] + sum_coup/nearest_neighbour
    end
end

ti = 0.0; tf = 500.0
dt = 0.01; tr = 495.0
tot_nts = floor(Int64,tf/dt)
tspan = (ti, tf)


fname2 = @sprintf("dn=solit_%dn_hetundel_pht.dat", number_of_solitary_states)
fout2 = open(fname2, "w")


# initial phases
rng = MersenneTwister(10)
u0 = 2pi*rand(rng, N)
h(p,t) = ones(N)


# system integration
prob = DDEProblem(kuramoto_del,u0,h,tspan; constant_lags=tau)
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-6, saveat=dt)

ph = sol[tot_nts]
for i = 1:size(ph,1)
	ph[i] = mod2pi(ph[i])
end

for ts = tr+dt:dt:tf
    global dt, tr, tf
    nts = floor(Int64,(ts-ti)/dt)
    y = sol[nts]
    yp = sol[nts-1]
    ef = (y - yp)/dt
    #push!(mf, (y-yp)/dt);
    #= for i = 1:size(y,1)
        push!(mf, (y[i] - yp[i]) / dt);
    end =#
	
    for i = 1:size(y,1)
	        
	y[i] = mod2pi(y[i])
        #yp[i] = mod2pi(yp[i])
    end

    writedlm(fout2, transpose(y), "\t")
end
close(fout2)



