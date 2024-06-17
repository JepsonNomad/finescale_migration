#### Define parameters ----
# Sampling requirements
minHr = 12
minFreq = paste0(minHr, " hours")
minDura = "15 weeks"
message(paste0("Min frequency: ", minFreq))
message(paste0("Min duration: ", minDura))

# migrateR parameters
startDateMigR = 304 # Day of year to start migrateR "years"
min.d = 500 # Distance between ranges (meters)
min.r = 21 # Duration of range occupancy
message(paste0("Min delta (for migrants): ", min.d, " meters"))
message(paste0("Min r (for migrants): ", min.r, " days"))

# The width of the window (minus 1) for the actual ssf analysis.
# This is the centerpoint of movement + 1/2 windowWidth and centerpoint - 1/2 windowWidth
# in days, if windStyl=="Fixed". Otherwise the function will take start- and endpoints of
# migration derived from migrateR.
windStyl = "Fixed"
windWide = 20

# Years to consider in the analysis
myYears = c(2003:2022)

# Set a seed for reproducibility
set.seed(80) # Colliseum opened to public

