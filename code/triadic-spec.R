# Specifications for triadic analysis
load('calc/mathrev-years.RData')

# Duration and number of intervals
dur <- 3
int <- 3
# Range of years to consider (data are noticeably incomplete for 2009)
ran <- range(setdiff(years, max(years)))

# Evenly-spaced intervals to include, by last year and index
yrs <- seq(ran[1] - 1 + dur, ran[2], floor((diff(ran) + 1 - int) / (int - 1)))
inds <- sapply(yrs, function(y) y - dur + 1 - ran[1] + 1)

# Remove unnecessary stuff
rm(int, ran)
