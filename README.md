# triadic

This repo contains code to reproduce the analysis in the paper [*Triadic analysis of affiliation networks*] [1].

[1]: http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=10081620&fileId=S2050124215000387

## Description

The paper describes an attempt to define and evaluate a new clustering coefficient for affiliation networks, establishing along the way several frameworks for understanding the diversity of possibilities and how they can be compared.

## Data

The code in this repo makes use of 24 years of data from [*Mathematical Reviews*] [2], available by request to the director of the American Mathematical Society. The data may require cleaning in order to match the format assumed by `code/mathrev2igraph.R` and `scripts/mathrev.R`. The code also uses data included in the [`bitriad` package] [3], specifically affiliation networks drawn from the studies [*Deep South: A Social Anthropological Study of Caste and Class*] [4] and [*Social Organization of an Urban Grants Economy: A Study of Business Philanthropy and Non-Profit Organizations*] [5].

[2]: http://www.ams.org/mr-database
[3]: https://github.com/corybrunson/bitriad
[4]: http://books.google.com/books?id=Q3b9QTOgLFcC
[5]: http://books.google.com/books?id=fR-LBQAAQBAJ

## Reproduction

The following steps should generate all calculations and visualizations used in the paper:

* Clone this repo to a machine with R installed (version >= 3.0.1).
* From within the `triadic` directory, execute `triadic.R`. Alternatively, open `triadic.R` in a text editor or [RStudio] [6] and execute the commands (which mostly source other files) in order. (The file begins by checking for and, as necessary, installing packages required for the analysis.)

Some of the calculations on *MR* data may take a *long* time (e.g. several hours on a 3.2 GHz Intel Core i5).

[6]: http://www.rstudio.com/

## Feedback

Any questions or suggestions on this code are more than welcome.
