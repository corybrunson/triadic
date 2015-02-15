# Compute tables of critical values for tau using all available methods
# Compare with these:
# http://www.york.ac.uk/depts/maths/tables/kendall.pdf
# http://faculty.washington.edu/heagerty/Books/Biostatistics/TABLES/Kendall.pdf
# http://books.google.com/books?id=IMbVyKoZRh8C&pg=PA691#v=onepage&q&f=false
require(tautable)

# Table dimensions
master.n <- 1:100
master.alpha <- as.vector(sapply(1:4, function(e) c(5, 2.5, 1) * 10 ^ (-e)))
# Construct tau table!
inv.crit.vals <- tau.crit.table(n = master.n, alpha = master.alpha,
                                incl.len = TRUE, stat = 'inv')
tau.crit.vals <- cbind(
    n = inv.crit.vals[, 1, drop = TRUE],
    (choose(inv.crit.vals[, 1], 2) - 2 * inv.crit.vals[, -1]) /
        choose(inv.crit.vals[, 1], 2)
)
# Save
save(inv.crit.vals, tau.crit.vals, file = 'calc/tau-tables.RData')

rm(list = ls())
