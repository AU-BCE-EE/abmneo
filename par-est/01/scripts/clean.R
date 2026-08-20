
# Get slurry mass

level <- level[level_scenario == 'a', ]
level <- merge(tanks[, .(tank, area)], level, by = 'tank')
level[, slurry_mass := 1000 * level * area]

mass <- level[, .(doy, slurry_mass, tank)]
names(mass) <- c('time', 'slurry_mass', 'tank')
# Drop rows with no mass change (do not need daily resolution)
mass <- mass[c(TRUE, diff(slurry_mass) != 0) | time == 1]

temp <- setorder(temp, tank, doy)
# Take only weekly temperatures
temp <- temp[doy %in% seq(1, 365, 7), ]

emis[, doy_start := as.integer(format(date_time_min, '%j'))]
emis[, doy_end := as.integer(format(date_time_max, '%j'))]
# Row 6 seems to have doy error!
emis[, doy_mid := as.integer(format(date_time, '%j'))]
# Change doy_start for single period that spanned years
emis[doy_start == 365, doy_start := 1]

