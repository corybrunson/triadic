# triadic

This repo contains code to reproduce the analysis in the paper [*Triadic analysis of affiliation networks*](http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=10081620&fileId=S2050124215000387).

## Description

The paper describes an attempt to define and evaluate a new clustering coefficient for affiliation networks, establishing along the way several frameworks for understanding the diversity of possibilities and how they can be compared.

## Data

The code in this repo makes use of 24 years of data from [*Mathematical Reviews*](http://www.ams.org/mr-database), available by request to the director of the American Mathematical Society. The data may require cleaning in order to match the format assumed by `code/mathrev2igraph.R` and `scripts/mathrev.R`. The code also uses data included in the [`bitriad` package](https://github.com/corybrunson/bitriad), specifically affiliation networks drawn from the studies [*Deep South: A Social Anthropological Study of Caste and Class*](http://books.google.com/books?id=Q3b9QTOgLFcC) and [*Social Organization of an Urban Grants Economy: A Study of Business Philanthropy and Non-Profit Organizations*](http://books.google.com/books?id=fR-LBQAAQBAJ).

## Reproduction

The following steps should generate (subject to differences in random seeds) all calculations and visualizations used in the paper:

* Clone this repo to a machine with R installed (version >= 3.3.2).
* From within the `triadic` directory, execute `triadic.R`. Alternatively, open `triadic.R` in a text editor or [RStudio](http://www.rstudio.com/) and execute the commands (which mostly source other files) in order. (The file begins by checking for and, as necessary, installing packages required for the analysis.)

Some of the calculations on *MR* data may take a *long* time (e.g. several hours on a 3.2 GHz Intel Core i5).

## Feedback

Any questions or suggestions on this code are more than welcome.
