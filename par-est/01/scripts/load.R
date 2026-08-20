# Load input data and emission measurements

tanks <- fread('../data/tanks.csv')
temp <- fread('../data/tank_temp_interp.csv')
emis <- fread('../data/meas_emis_ave_samp.csv')
level <- fread('../data/level.csv')
