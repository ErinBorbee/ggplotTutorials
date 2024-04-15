# ggplot Tutorials
This repository contains jupyter notbooks with code tutorials on building different types of graphs in R. These tutorials were originally designed for students in a data visualization course at Texas State University, who had little to no experience with ggplot prior to the course. The majority of these tutorials use the example datasets built into R like the `iris` and `mtcars` datasets. For the ones that use datasets not provided by R the data input is also provided in this repository. 

#### ggplot Introduction
This code goes through the basic syntax of ggplot and how to customize theme and appearance of graphs using boxplots. The code then moves on to show how we can use multiple data layers in a single plot and finishes by showing how to export figures from R using the `ggsave()` function.

#### Distribution plots
Distribution plots include examples of boxplots, violin plots, histograms, and density plots. Within the code for these we also explore more ways of customizing appearance including using preset palettes and custom colors for variables in the graph. 

### Distribution plots & Correlation plots
This set of code concludes distribution plots with ridgeline plots and rain cloud plots and begins correlation plots with scatterplots, regressions, and bubble plots. 

### Ordinations & Evolution plots
This code begins by introducing how to use `facet_grid()` and `facet_wrap()` to separate data into multiple panels by a single variable. We then move into ordinations and demostrate how to plot a PCA from the iris dataset. Finally, the script finishes by demostrating how to plot a line plot and an area plot using the economics dataset.

#### Mapping
This script shows how to plot maps in ggplot using built in map datasets first and then using high resolution datasets from GADM. Using GADM, we plot multiple layers of data, country, state, and county levels. In plotting these maps, we demonstrate how to set your x and y axis ranges to crop your figure to a specific region. Finally, we downloaded river and terrain shapefiles for our geographic area and plotted those on top of our map to create a more detailed image.

Texas river and terrain data source: https://www.twdb.texas.gov/mapping/gisdata.asp 

#### Combining plots & Interactive plotting
These commands show how to use `ggarrange()` to incorpporate multiple graphs into a single figure and label them. The second portion of this script shows how to use the Plotly package to transform your ggplots into interactive plots that provide data or statistical information when hovering over the graph and datapoints. 

