1. Compile the program by running "./compile.sh"

2. Generate the script for starting one search over P, N, and the number of threads N_threads you want to use. This can be done by running "python ./gen_instance.py P N N_threads" with the particular parameters

3. Start the search by running "./start_Px_Nx_threadsx.sh" (make sure the progam is executable. If not, run "chmod u+x start_Px_Nx_threadsx.sh" first)

4. If all programs terminate, you can obtain the solution file by running "cat Px_Nx/*.sol_P_N > Px_Nx/solutions.sol_P_N"
