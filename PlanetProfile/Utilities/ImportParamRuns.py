import scipy.io
import pandas as pd
import numpy as np

class Entry:
    def __init__(self, frequency, amplitude, phase, color, D_ocean_km, D_seafloor_km,
                 k_ocean_Sm, rho_kgm3, k_ions_Sm, Bfield, filename):
        self.frequency = frequency
        self.amplitude = amplitude
        self.phase = phase
        self.color = color
        self.D_ocean_km = D_ocean_km
        self.D_seafloor_km = D_seafloor_km
        self.k_ocean_Sm = k_ocean_Sm
        self.rho_kgm3 = rho_kgm3
        self.k_ions_Sm = k_ions_Sm
        self.Bfield = Bfield
        self.filename = filename

class Run:
    def __init__(self):
        self.entries = []

    def add_entry(self, entry):
        self.entries.append(entry)

filename = 'JointCallistoFixedCond100kmIonospheres_10000S_withZero.mat'
# Load the .mat file
mat = scipy.io.loadmat(filename)

# Extract the 'thisrun' structure array
thisrun = mat['thisrun']
print(thisrun.dtype)
# Initialize a list to store all runs
runs = []
# Iterate over the runs
for field_name in thisrun.dtype.names:
    run_data = thisrun[field_name][0, 0]  # Access the nested structure
    run = Run()

    # Iterate over the 480 entries in each run
    for i in range(len(run_data['frequency'][0])):
        entry = Entry(
            frequency=run_data['frequency'][0][i],
            amplitude=run_data['amplitude'][0][i],
            phase=run_data['phase'][0][i],
            color=run_data['color'][0][i],
            D_ocean_km=run_data['D_ocean_km'][0][i][0],
            D_seafloor_km=run_data['D_seafloor_km'][0][i][0],
            k_ocean_Sm=run_data['k_ocean_Sm'][0][i][0],
            rho_kgm3=run_data['rho_kgm3'][0][i][0],
            k_ions_Sm=run_data['k_ions_Sm'][0][i][0],
            Bfield=run_data['Bfield'][0][i],
            filename=run_data['filename'][0][i]
        )
        run.add_entry(entry)

    runs.append(run)



# Convert a single run's entries to a DataFrame
def entries_to_dataframe(entries):
    data = {
        'frequency': [e.frequency for e in entries],
        'amplitude': [e.amplitude for e in entries],
        'phase': [e.phase for e in entries],
        'color': [e.color for e in entries],
        'D_ocean_km': [e.D_ocean_km for e in entries],
        'D_seafloor_km': [e.D_seafloor_km for e in entries],
        'k_ocean_Sm': [e.k_ocean_Sm for e in entries],
        'rho_kgm3': [e.rho_kgm3 for e in entries],
        'k_ions_Sm': [e.k_ions_Sm for e in entries],
        'Bfield': [e.Bfield for e in entries],
        'filename': [e.filename for e in entries]
    }
    return pd.DataFrame(data)


# Convert all runs to a list of DataFrames
runs_dfs = [entries_to_dataframe(run.entries) for run in runs]
# Now 'runs' is a list of Run objects, each containing 480 Entry objects
