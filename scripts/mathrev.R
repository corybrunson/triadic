# Load the database as a data frame
classes <- c('numeric', 'integer', rep('character', 3))
mathrev <- data.frame()
for(i in 0:4) mathrev <- rbind(mathrev, read.table(
    paste0('data/distdata1985-2009/distdata',
           1985 + i * 5, '-', 1989 + i * 5, '.csv'),
    header = FALSE, sep = ',', colClasses = classes))
names(mathrev) <- c('index', 'year', 'authors', 'pclass', 'sclass')

# Ensure that only years within range are included
mathrev <- mathrev[mathrev$year >= 1985, ]
mathrev <- mathrev[mathrev$year <= 2009, ]

# Remove publications attributed to the artifact author string '00a00a00a'
mathrev <- mathrev[-grep('00a00a00a', mathrev$authors), ]

# Remove publications attributed to '02c50a20a' (et al)
# (which necessarily have incomplete author lists)
mathrev <- mathrev[-grep('02c50a20a', mathrev$authors), ]

save(mathrev, file = 'calc/mathrev.RData')

years <- sort(unique(mathrev$year))
save(years, file = 'calc/mathrev-years.RData')

rm(list = ls())
