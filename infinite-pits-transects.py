from landlab import RasterModelGrid, imshow_grid
import numpy as np
from landlab.components import LakeMapperBarnes, FlowAccumulator
import os, csv, copy
import pickle

# Define the grid dimensions
nc = 300    #nodes in x-direction
nr = 300    #nodes in y-direction
dx = 0.025   #horizontal node spacing in meters

# Define constants
h_c = 0.08            		#characteristic debris thickness in meters
D_list = [0.001, 0.005, 0.01]   #debris diffusivity constant in meters per day, less than 0.01 m/day for grains smaller than 0.01 m diameter
s_c = 1.15                		#Critical hillslope steepness in m/m

# Define basic model parameters
initial_ice_thickness = 100.     #initial ice thickness in meters (arbitrary)
total_time = 30                 #simulated runtime in melt-days (2-3 melt seasons)
dt = 0.01                     #time step interval in days
tsteps = int(total_time/dt)

b_list = [0.04]       #bare ice melt rate in meters per day
debris_depth = 100. 	#pit depth in meters
radius = 0.25    #pit radius in meters

pit = 'yes' #Type 'yes' if the debris patch exists in a pit, otherwise type 'no'

# Create a directory to hold the experiments
dirname = "infinite-pit-experiments-trans/"

if not os.path.exists(dirname):
    os.mkdir(dirname)

# Create a text file to hold the parameters from each experiment
indx_filename = "infinite-pits-trans_index.txt"

if not os.path.exists(os.path.join(dirname, indx_filename)):
        with open(os.path.join(dirname, indx_filename), 'w') as fn:
            fn.write("simulation_number\tD\th_c\tS_c\tb_0\tdebris_radius\tdebris_depth\tcone_height\tcone_width\tdays\n")

sim = 1


for D in D_list:
    for b_0 in b_list:

        simFileName = "infinite-pit-trans_"

        """
        RUN THE MODEL
        """
            # Create the model grid
        mg = RasterModelGrid((nr, nc), xy_spacing=dx)
        s = mg.add_zeros('ice_thickness', at='node')
        h = mg.add_zeros('debris_thickness', at='node')
        z = mg.add_zeros('topographic__elevation', at='node')
        qs = mg.add_zeros('debris_flux', at='link')
        melt = mg.add_zeros('melt', at='node')

        # Set the boundary conditions
        mg.set_closed_boundaries_at_grid_edges(False, False, False, False)

        """
        DEFINE THE MOULIN COORDINATES
        """
        # Assign initial ice thickness
        s += initial_ice_thickness

        # Center the debris patch in the model grid
        x_m = (nc) * dx / 2
        y_m = (nr) * dx / 2

        # Find the grid node ids that lie within the circular moulin
        dists = mg.calc_distances_of_nodes_to_point([x_m, y_m])

        debris_patch_nodes = []
        i = 0
        for distance in dists:
            if distance <= radius:
                debris_patch_nodes.append(i)
            i += 1

        debris_area = len(debris_patch_nodes) * dx**2

        # Define the moulin depth, remove the ice, fill with debris
        for ids in debris_patch_nodes:
            if pit == "yes":
                s[ids] -= debris_depth
            h[ids] += debris_depth 

        # Update initial topographic elevation
        z = s + h
        mg.at_node['topographic__elevation'] = z

        # Create empty lists to store cone geometry
        cone_height = ['cone_height_m']
        ice_height = ['ice_height_m']
        cone_width = ['cone_width_m']
        times = ['time_days']

        z_transects = {}
        s_transects = {}

        z_transects['x'] = mg.x_of_node
        s_transects['x'] = mg.x_of_node
        z_transects[0] = mg.node_vector_to_raster(z)[nr//2, :]
        s_transects[0] = copy.deepcopy(mg.node_vector_to_raster(s)[nr//2, :])

        for ts in range(tsteps):
            # Calculate sub-debris melt with hyperbolic model
            melt = b_0 * dt * h_c / (h_c + h)

            # Update the ice elevation
            s -= melt

            # Update the topographic elevation
            z = h + s

            # Calculate the topographic gradient
            grad = mg.calc_grad_at_link(z)

            # Calculate slope magnitudes
            slp = mg.calc_slope_at_node(z)
            slope = mg.map_mean_of_link_nodes_to_link(slp)

            # Calculate the "transportable" debris
            s_fill = s.copy() #create a duplicate ice surface
            fa = FlowAccumulator(mg)
            lmb = LakeMapperBarnes(mg, method = "Steepest", fill_flat = True, surface = s, fill_surface = s_fill, track_lakes = True)
            lmb.run_one_step()

            s_fill = s + lmb.lake_depths #calculate the elevation of the filled ice surface
            h_trans = z - s_fill #calculate transportable debris thickness
            h_link = mg.map_mean_of_link_nodes_to_link(h_trans) #transfer transportable debris to links

            # Calculate debris flux (nonlinear approximation)
            qs[mg.active_links] = D * (1 - np.exp(-h_link[mg.active_links] / h_c)) * grad[mg.active_links] / (1 - (grad[mg.active_links] / s_c)**2)
            flux = qs[mg.active_links]

            # Break the loop if flux "blows up", write warning to file
            if qs.max() > 2.:
                print("Flux getting big!")

            dhdt = mg.calc_flux_div_at_node(qs) * dt

            # Update the debris thickness
            h = h.copy() + dhdt

            # Update the topographic elevation
            z[mg.core_nodes] = s[mg.core_nodes] + h[mg.core_nodes]

            #Ensure that all fields are updated
            mg.at_node['ice_thickness'] = s
            mg.at_node['debris_thickness'] = h
            mg.at_node['topographic__elevation'] = z

            #Save data at model start, each half-day, when apex thickness less than 1 cm, or final timestep

            if (ts+1)*dt % 0.5 == 0 or ts == 0 or mg.node_vector_to_raster(h)[nr//2, nc//2] <= 0.01 or ts+1 == tsteps:
                np.savetxt(os.path.join(dirname, simFileName + str(sim) +'_flux-ts=' + str(ts+1) + '.txt'), flux, fmt='%.4e')

                max_height = z.max()-z[0] #Calculate dirt cone height (ice + debris) relative to bare ice

                #Append cone and ice heights to lists
                cone_height.append(np.format_float_positional(max_height, precision = 3))
                ice_height.append(np.format_float_positional(mg.node_vector_to_raster(s)[nr//2,nc//2]-s[0], precision = 3)) #Calculate ice cone height relative to bare ice

                #Calculate time in days
                time = (ts+1)*dt
                times.append(time) #Append model time to list

                #Calculate the cone width
                debris = mg.node_vector_to_raster(h)[nr//2, :]
                w = np.array(np.where(debris > 0.001)) #Assume the cone exists where debris is thicker than 1 mm
                c_width = (int(w[:,-1]) - int(w[:,0]))*dx
                cone_width.append(np.format_float_positional(c_width, precision = 3)) #Append width of the cone to lists

                #Create a transect

                z_transects[time] = mg.node_vector_to_raster(z)[nr//2, :]
                s_transects[time] = copy.deepcopy(mg.node_vector_to_raster(s)[nr//2, :])

                #Save the data when the model terminates    
                if mg.node_vector_to_raster(h)[nr//2, nc//2] <= 0.01 or ts+1 == tsteps:
                    time = (ts+1)*dt

                    #Create csv of daily data from this simulation
                    sim_filename = os.path.join(dirname, simFileName + str(sim) +'.csv')
                    data = [times, cone_height, ice_height, cone_width]
                    with open(sim_filename, "w", newline = "") as f:
                        writer = csv.writer(f)
                        writer.writerows(data)

                    #Create topographic transects csv    
                    z_trans_fname = os.path.join(dirname, simFileName + str(sim) + "_z_trans.csv")
                    with open(z_trans_fname, "w") as f:
                        writer = csv.writer(f)
                        key_list_z = list(z_transects.keys())
                        limit = len(mg.node_vector_to_raster(z)[nr//2, :])
                        writer.writerow(z_transects.keys())
                        for i in range(limit):
                            writer.writerow([z_transects[x][i] for x in key_list_z])
                    #Create ice surface transects csv
                    s_trans_fname = os.path.join(dirname, simFileName + str(sim) + "_s_trans.csv")
                    with open(s_trans_fname, "w") as f:
                        writer = csv.writer(f)
                        key_list_s = list(s_transects.keys())
                        limit = len(mg.node_vector_to_raster(s)[nr//2, :])
                        writer.writerow(s_transects.keys())
                        for j in range(limit):
                            writer.writerow([s_transects[y][j] for y in key_list_s])


                    # Write simulation parameters, cone height, and cone width to the index file            
                    with open(os.path.join(dirname, indx_filename), 'a') as fn:
                        #fn.write("simulation_number\tD\th_c\tS_c\tb_0\tdebris_radius\tdebris_depth\tcone_height\tcone_width\tdays\n")
                        fn.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sim, D, h_c, s_c, b_0, radius, debris_depth, np.round(max_height, 3), np.round(c_width, 3), np.round(ts*dt, 2)))

                    sim += 1
                    break