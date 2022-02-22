from sys import argv

P = int(argv[1])
N = int(argv[2])
N_threads = int(argv[3])

Q = P**N
L = int(Q/N_threads)

filename = "start_P"+repr(P)+"_N"+repr(N)+"_threads"+repr(N_threads)+".sh"
file = open(filename, "w")
file.write('#!/bin/bash\n\n')
file.write('mkdir -p P' + repr(P) + '_N' + repr(N) + '\n')

for i in range(N_threads-1):
	file.write('./planar_monomials P' + repr(P) + '_N' + repr(N) + '/solutions_' + repr(i+1) + '.sol_' + repr(P) + '_' + repr(N) + ' ' + repr(P) + ' ' + repr(N) + ' ' +repr(i*L + 1) + ' ' + repr((i+1)*L) + ' > P' + repr(P) + '_N' + repr(N) + '/planar_monomials_' + repr(P) + '_' + repr(N) + '_' +repr(i+1) + '.out &\n')

file.write('./planar_monomials P' + repr(P) + '_N' + repr(N) + '/solutions_' + repr(N_threads) + '.sol_' + repr(P) + '_' + repr(N) + ' ' + repr(P) + ' ' + repr(N) + ' ' +repr((N_threads-1)*L + 1) + ' ' + repr(Q-1) + ' > P' + repr(P) + '_N' + repr(N) + '/planar_monomials_' + repr(P) + '_' + repr(N) + '_' +repr(N_threads) + '.out &\n')
	
file.close()
