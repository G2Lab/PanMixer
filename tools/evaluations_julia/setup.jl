""" Import package
"""

import Pkg
Pkg.add("CategoricalArrays")
Pkg.add("Distances")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Statistics")
Pkg.add("Plots")
Pkg.add("StatsPlots")
Pkg.add("ArgParse")
Pkg.add("Gadfly")
Pkg.precompile()