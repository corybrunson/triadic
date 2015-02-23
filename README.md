# triadic

This repo contains code to reproduce the analysis in the paper *Triadic analysis of affiliation networks*.

## Description

The paper describes an attempt to define and evaluate a new clustering coefficient for affiliation networks, establishing along the way several frameworks for understanding the diversity of possibilities and how they can be compared.

## Data

The code in this repo makes use of 24 years of data from [*Mathematical Reviews*] [2], available by request to the director of the American Mathematical Society. If the data obtained is formatted differently, some changes to `code/mathrev2igraph.R` and `scripts/mathrev.R` may be in order. The code also uses data included in the [`bitriad` package] [3], specifically affiliation networks drawn from the studies [*Deep South: A Social Anthropological Study of Caste and Class*] [4] and [*Social Organization of an Urban Grants Economy: A Study of Business Philanthropy and Non-Profit Organizations*] [5].

[2]: http://www.ams.org/mr-database
[3]: https://github.com/corybrunson/bitriad
[4]: http://books.google.com/books?id=Q3b9QTOgLFcC
[5]: http://books.google.com/books?id=fR-LBQAAQBAJ

## Reproduction

To generate all calculations and visualizations used in the paper, clone this repo to a machine with R installed (version >= 3.0.1) and, from within the `triadic` directory, execute `triadic.R`. Alternatively, open `triadic.R` in a text editor or [RStudio] [6] and execute the commands (mostly sourcing other files) in order. Some of the intensive calculations may take some time to complete.

[6]: http://www.rstudio.com/

## Feedback

Any questions or suggestions on this code are more than welcome.
