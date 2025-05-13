import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import uproot
import cygno as cy
import pandas as pd


# pkl maker process a run to create a pkl file for faster processing it takes as variables the number of the run you want to process and the dataframe name (look it up later)

def pkl_maker(run_number, df_name=''):
    try:
        import numpy as np
        import pandas as pd
        import uproot
        import os
        
        file_url = '/jupyter-workspace/private/PSF_test1/reco_run{:5d}_3D.root'.format(run_number)
        file_out_name = '/jupyter-workspace/private/run_data/reco_run{:5d}_3D.pkl.gz'.format(run_number)
        
        tf = uproot.open(file_url)
        names = tf["Events;1"].keys()
        branch_data = {}         # python dictionary
        for i, name in enumerate(names):
            print(i,name)
            var = tf["Events;1/"+name].array(library='np')        # this loads the branch with that name as a numpy array
            if var[0].ndim == 0:         # check if the variable is a scalar instead of an array
                branch_data[name] = np.hstack(var)       # this flattens the array itself if the first element is a scalar
            else:
                branch_data[name] = var
            #print('branch imported correctly')
        tf.close()
        
        df_all = pd.DataFrame(branch_data)
        df_all.to_pickle(file_out_name, compression={'method': 'gzip', 'compresslevel': 1})
        del df_all
        
        if(os.path.isfile(file_out_name)):    # this checks if the file was correctly created
            if(df_name != ''):
                dfk = pd.read_json(df_name)
                dfk.at[run_number, 'pkl'] = 1
                dfk.to_json(df_name)
                del dfk
            file = open('logbook.txt', 'a')
            file.write('made pkl of run {:5d}\n'.format(run_number))
            file.close()
            
    except Exception as e:
        file = open('logbook.txt', 'a')
        file.write('run {:5d} ERROR >>> {}\n'.format(run_number, e))
        file.close()


# this class takes the data from the pickle file and imports it as numpy arrays for processing and analysis


class ExtractData:
    def __init__(self, run_number, run_path):
        self.run_number = run_number
        self.run_path = run_path
        self.max_roundness = 0.7
        self.data = self._set_data()
        self._muon_mask = self._muon_filter()
        self.redpix = self._set_redpix()

    def _set_data(self):
        path = f'{self.run_path}/reco_run{self.run_number:05d}_3D.pkl.gz'
        return pd.read_pickle(path, compression = 'gzip')

    

    # This function retrieves the redpix data, indexes it correctly and retrunes the tracks that pass the filter

    def _set_redpix(self):
        rx = self.data['redpix_ix'].to_numpy() - 896
        ry = self.data['redpix_iy'].to_numpy()
        rz = self.data['redpix_iz'].to_numpy()
        Id = self.data['sc_redpixIdx'].to_numpy()

        redx, redy, redz = [], [], []
        for n_image in range(len(Id)):  # Loop over images
            cluster_x, cluster_y, cluster_z = [], [], []
    
            for cluster in range(len(Id[n_image])):
                start_index = int(Id[n_image][cluster])

                # Ignore removed clusters
                if start_index != -1 and self._muon_mask[n_image][cluster]:
                    if cluster < len(Id[n_image]) -1:
                        end_index = int(Id[n_image][cluster +1])-1
                    else:
                        end_index = len(rx[n_image])
                    cluster_x.append(rx[n_image][start_index:end_index])
                    cluster_y.append(ry[n_image][start_index:end_index])
                    cluster_z.append(rz[n_image][start_index:end_index])
    
            redx.append(cluster_x)
            redy.append(cluster_y)
            redz.append(cluster_z)
            del cluster_x, cluster_y, cluster_z
            
        return redx, redy, redz

    

    # this is the method that computes the boolean mask I use to select the muons

    def _muon_filter(self):
        boolean_mask = [(self.data['sc_integral'][n_image] > 20000) & 
                        #(self.data['sc_integral'][n_image] < 50000) & 
                        #(self.data['sc_xmax'][n_image] - self.data['sc_xmin'][n_image] < 500) &
                        (self.data['sc_width'][n_image] / self.data['sc_length'][n_image] < 0.7) &
                        (self.data['sc_width'][n_image] < 70) &
                        (self.data['sc_length'][n_image] > 1200)
                        #(self.data['sc_ymax'][n_image] - self.data['sc_ymin'][n_image] > 500)
                        for n_image in range(len(self.data['sc_integral']))]
        #boolean_mask = [self.data['sc_integral'][n_image] < 2000 for n_image in range (len(self.data['sc_integral']) ) ]
        return boolean_mask


    # This function retrieves the cluster position and integral

    def get_sc_variables(self):
        x = np.concatenate(self.data['sc_xmean'].to_numpy(), axis=0) - 896
        y = np.concatenate(self.data['sc_ymean'].to_numpy(), axis=0)
        z = np.concatenate(self.data['sc_integral'].to_numpy(), axis=0)
        
        return x, y, z



    #This is a method to retrieve the full array whenever needed
    
    def get_raw_data(self):
        return self.data



    def get_redpix(self):
        return self.redpix



    # This is a method to plot the downsampled images of selected cluster tracks
    # n_start is the starting run number
    # n_end is the number of the last run, note that the range function will actually stop at n_end-1
    # downsample_value instead tells the function how much you want to downsample the tracks, the value must be an integer >= 1
    
    def plot_data(self, n_start, n_end, downsample_value):
        fig, ax = plt.subplots(figsize=(9,8))

        for image in range(n_start, n_end):
            #print('processing image number', image)
    
            for cluster in range(len(self.redpix[0][image])):  
                #print('cluster number', cluster)
        
                n = int(max(1, len(self.redpix[0][image][cluster]) / downsample_value))  # Ensure at least 1 point
       
                downsample_indices = np.random.choice(len(self.redpix[0][image][cluster]), n, replace=False)
        
                x_sampled = self.redpix[0][image][cluster][downsample_indices]
                y_sampled = self.redpix[1][image][cluster][downsample_indices]
                z_sampled = self.redpix[2][image][cluster][downsample_indices]

                pl = ax.scatter(x_sampled, y_sampled, s=0.2, c=z_sampled)

        ax.set_xlim(0, 2304)
        ax.set_ylim(0, 2304)
        ax.set_xlabel('x (pixel)')
        ax.set_ylabel('y (pixel)')
        ax.set_title(f'Filtered tracks in runs: {n_start} to {n_end-1}')

        colorbar = fig.colorbar(pl, ax=ax)
        colorbar.set_label('Intensity')
        plt.show()
        del x_sampled, y_sampled, z_sampled