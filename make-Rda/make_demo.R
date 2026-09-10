# Creates demo inputs

storage_ex_dig_01 <- list(
  type = 'series',
  dat = data.frame(
    time = c(0, 100, 200, 365), 
    slurry_mass = c(4, 0.5, 1.8, 4) * 500 * 1000
  ),
  storage_depth = 4,     
  area = 500,              
  temp_C = 15,
  resid_enrich = 0.64 
)

save(storage_ex_dig_01, file = '../data/storage_ex_dig_01.rda')

storage_ex_pig_01 <- list(
  type = 'regular',
  slurry_prod_rate = 5, 
  storage_depth = 2,     
  area = 0.65,              
  temp_C = 20,
  resid_enrich = 0.6,
  slurry_mass = 6.5,     
  resid_depth = 0.05,      
  empty_int = 28
)

save(storage_ex_pig_01, file = '../data/storage_ex_pig_01.rda')

