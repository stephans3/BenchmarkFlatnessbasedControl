using BellBruno

# Write Bell polynomials
N_der = 10;
path_2_folder = "results/bell_results/"
bp = bell_poly( N_der;
                save_on_disk=true, 
                path_to_folder=path_2_folder, 
                print_iteration=true)
                