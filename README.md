<br>
<div align="center">
    <img src="data/logo.svg?maxAge=0" width="30%">
</div>

This package implements a [Makie](https://makie.juliaplots.org/stable/index.html)
recipe called `escherplot` that plots maps of metabolic models resembling its [namesake
GUI](https://escher.github.io/#/). Its primary purpose is to facilitate the generation of
high quality metabolic maps from within Julia.

## Plot the core metabolism of E. coli with fluxes
Here [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl) is used to estimate fluxes of
reactions using the [metabolic model of iJO1366](http://bigg.ucsd.edu/models/iJO1366), an E. coli
metabolic model. The associated map can be downloaded from the [Escher website](https://escher.github.io/#/). If you want to run these examples, please download the associated models and maps, and 
place them in the `data`directory.
```julia
using Escher, CairoMakie, ColorSchemes
using COBREXA, Tulip
using JSONFBCModels, AbstractFBCModels
using Clustering

# use COBREXA to generate a flux distribution using the associated model
model = load_model(joinpath(pkgdir(Escher), "data", "iJO1366-model.json"))
#=
Bin fluxes for display purposes - assigning colors to edges needs to be done
manually. The binning uses kmeans clustering on logged fluxes due to the large
differences between fluxes.
=#
sol = flux_balance_analysis(model; optimizer=Tulip.Optimizer).fluxes
rids = string.(keys(sol))
fluxes = values(sol)
logged_fluxes = log.(abs.(fluxes) .+ 1e-8)
clusters = kmeans(logged_fluxes', 9)
centers = Dict(j=>i for (i, j) in enumerate(sortperm(clusters.centers'[:])))
order = [centers[i] for i in assignments(clusters)]

rc = Dict(rid => ColorSchemes.RdYlBu_9[10-k] for (rid, k) in zip(rids, order)) # map reaction id to color

f = Figure(size = (1200, 800));
ax = Axis(f[1, 1]);
escherplot!(
    ax,
    joinpath(pkgdir(Escher), "data", "iJO1366-map.json");
    reaction_edge_colors = rc,
)
hidexdecorations!(ax)
hideydecorations!(ax)
f
```
<br>
<div align="center">
    <img src="data/iJO1366-map.svg?maxAge=0" width="80%">
</div>

## Overlay even more data with colors and node/edge sizes
The previous example only highlighted reactions according to their fluxes. `Escher.jl` also
allows you to control the size and colors of the nodes, as well as the size of the reaction
edges. This time we will use a [smaller "core" model](http://bigg.ucsd.edu/models/e_coli_core)
with its associated map.
```julia
using Escher, CairoMakie, ColorSchemes
using COBREXA, Tulip

model = load_model(joinpath(pkgdir(Escher), "data", "core-model.json"))
sol = Dict(string(k) => v for (k,v) in flux_balance_analysis(model; optimizer=Tulip.Optimizer).fluxes)

# Find min and max absolute fluxes for normalization
maxflux = maximum(abs.(values(sol)))
minflux = minimum(abs.(values(sol)))

# Scale width of reaction edges to fluxes
width_interp(x) = 2 + 5 * (abs(x) - minflux) / (maxflux - minflux) # widths between 2 and 5
re = Dict(k => width_interp(v) for (k, v) in sol) # map reaction id to reaction edge width

# Scale color of reaction edges to fluxes (manually binned)
color_interp(x) = begin
    normed_x = (abs(x) - minflux) / (maxflux - minflux)
    if 0 <= normed_x < 0.01
        ColorSchemes.RdYlBu_4[4]
    elseif 0.01 <= normed_x < 0.25
        ColorSchemes.RdYlBu_4[3]
    elseif 0.25 <= normed_x < 0.5
        ColorSchemes.RdYlBu_4[2]
    else
        ColorSchemes.RdYlBu_4[1]
    end
end
rc = Dict(k => color_interp(v) for (k, v) in sol) # map reaction id to reaction edge color

# metabolite node colors
mc = Dict(
    k => ColorSchemes.:Dark2_7[v] for
    (k, v) in zip(AbstractFBCModels.metabolites(model), rand(1:7, AbstractFBCModels.n_metabolites(model)))
)

# metabolite node sizes
ms = Dict(k => v for (k, v) in zip(AbstractFBCModels.metabolites(model), rand(3:10, AbstractFBCModels.n_metabolites(model))))

# Normal Makie plotting features all work (escherplot is a full recipe)
f = Figure(size = (1200, 800));
ax = Axis(f[1, 1]);
escherplot!(
    ax,
    joinpath(pkgdir(Escher), "data", "core-map.json");
    reaction_edge_widths = re,
    reaction_edge_colors = rc,
    metabolite_node_colors = mc,
    metabolite_node_sizes = ms,
)
hidexdecorations!(ax)
hideydecorations!(ax)
f
```
This results in:
<br>
<div align="center">
    <img src="data/map.svg?maxAge=0" width="80%">
</div>

## Attributes
The `escherplot` recipe exposes a number of custom attributes that can be used to modify the
basic metabolic map figure. The names are all self-explanatory but some comments are
provided for clarity.
```julia
metabolite_identifier = "bigg_id"
metabolite_show_text = false
metabolite_text_size = 4
metabolite_primary_node_size = 5 # fallback size
metabolite_secondary_node_size = 3 # fallback size
metabolite_node_sizes = Dict{String,Any}()
metabolite_node_colors = Dict{String,Any}()
metabolite_node_color = :black # fallback color
metabolite_text_color = :black
reaction_identifier = "bigg_id"
reaction_show_text = false
reaction_show_name_instead_of_id = false
reaction_text_size = 4
reaction_text_color = :black
reaction_edge_colors = Dict{String,Any}() # actual color
reaction_edge_color = :black # fallback color
reaction_edge_widths = Dict{String,Any}() # actual edge width
reaction_edge_width = 2.0 # fallback width
reaction_arrow_size = 6
reaction_arrow_head_offset_fraction = 0.5 # between 0 and 1
reaction_directions = Dict{String,Tuple{Dict{String,Number},Symbol}}() # rid => (reaction stoichiometry, :f or :r)
annotation_show_text = false
annotation_text_color = :black
annotation_text_size = 12
```
Note, if `reaction_edge_colors` or `reaction_edge_widths` are supplied but missing an id
that is present in the map, the associated edge will be dotted. 

## More examples
These examples all use the same data as the second example, but demonstrate the use of
different attributes. For brevity it is assumed that the functions below are inserted as
indicated here:
```julia
f = Figure(size = (1200, 800));
ax = Axis(f[1, 1]);
###### PLOT FUNCTION
hidexdecorations!(ax)
hideydecorations!(ax)
f
```

### Basic plot
Basic plot showing only the edges (reactions) and nodes (metabolites).
```julia
escherplot!(ax, joinpath(pkgdir(Escher), "data", "core-map.json"))
```
<br>
<div align="center">
    <img src="data/basic-map.svg?maxAge=0" width="80%">
</div>

### Adding labels
Basic plot showing the labels of the nodes and reactions.
```julia
escherplot!(
    ax,
    joinpath(pkgdir(Escher), "data", "core-map.json");
    metabolite_show_text=true,
    reaction_show_text=true,
    annotation_show_text = true,
)
```
<br>
<div align="center">
    <img src="data/labels-map.svg?maxAge=0" width="80%">
</div>

### Missing data
If incomplete reaction data (edge colors or widths) are supplied the missing reaction edges
are dotted.
```julia
rc = Dict("FBA" => ColorSchemes.RdYlBu_4[4],
    "PFK" => ColorSchemes.RdYlBu_4[3],
    "PEP" => ColorSchemes.RdYlBu_4[2],
    "PYK" => ColorSchemes.RdYlBu_4[1])
escherplot!(
    ax,
    joinpath(pkgdir(Escher), "data", "core-map.json");
    reaction_edge_colors = rc,
)
```
<br>
<div align="center">
    <img src="data/missing-map.svg?maxAge=0" width="80%">
</div>

### Reaction directions 
It is also possible to add direction arrows to reactions through the
`reaction_directions` attribute. It is a dictionary, which maps reaction ids to
reaction stoichiometry of the model used to simulate fluxes, and a symbol
`:forward`, `:backward`, `:bidirectional`. Arrows are then placed on reactions
in the direction relative to the supplied stoichiometry of the reaction.
```julia
rd = Dict(
    "PGM" => (Dict("3pg_c" => 1, "2pg_c" => -1), :backward),
    "PYK" => (Dict("pep_c" => -1, "adp_c" => -1, "h_c" => -1.0, "atp_c" => 1.0, "pyr_c" => 1), :forward),
    "ENO" => (Dict("2pg_c" => 1, "pep_c" => -1, "h2o_c" => -1), :bidirectional),
)

escherplot!(
    ax,
    joinpath(pkgdir(Escher), "data", "core-map.json");
    reaction_directions = rd,
    reaction_arrow_size = 12,
    reaction_show_text = true,
)
hidexdecorations!(ax)
hideydecorations!(ax)
```
<br>
<div align="center">
    <img src="data/directions-map.svg?maxAge=0" width="80%">
</div>

### Map dimensions
Finally, the original map dimensions can be queried by calling `get_resolution`. 

```julia
h, w, x, y = get_resolution(joinpath(pkgdir(Escher), "data", "core-map.json"))
f = Figure(resolution = (x, y))
...
```
