<br>
<div align="center">
    <img src="data/logo.svg?maxAge=0" width="30%">
</div>

This package implements a single [Makie](https://makie.juliaplots.org/stable/index.html)
recipe called `escherplot` that plots maps of metabolic models resembling its [namesake
GUI](https://escher.github.io/#/). It's primary purpose is to facilitate the generation of
high quality metabolic maps from within Julia.

## Example
```julia
using Escher, CairoMakie, ColorSchemes
using COBREXA, Gurobi

# use COBREXA to generate a flux distribution using the associated model
sol = parsimonious_flux_balance_analysis_dict(
    load_model("core-model.json"),
    Gurobi.Optimizer,
)

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

# Normal Makie plotting features all work (escherplot is a full recipe)
f = Figure();
ax = Axis(f[1, 1]);
escherplot!(ax, "core-map.json"; reaction_edge_widths = re, reaction_edge_colors = rc)
hidexdecorations!(ax)
hideydecorations!(ax)
f
CairoMakie.FileIO.save("map.pdf", f) # save figure
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
metabolite_identifier = "bigg_id" # identifier used to extract metabolite ids
metabolite_show_text = false # show the metabolite identifier
metabolite_text_size = 4 # metabolite identifier text size
metabolite_primary_node_size = 5 # primary metabolite node size
metabolite_secondary_node_size = 3 # secondary/co-factor metabolite node size
metabolite_node_color = :black # Color used for all metabolite nodes
metabolite_text_color = :black # Color used for all metabolite text
reaction_identifier = "bigg_id" # identifier used to extract reaction ids
reaction_show_text = false # show the reaction identifier
reaction_show_name_instead_of_id = false # show the reaction name instead of id
reaction_text_size = 4 # reaction identifier text size
reaction_text_color = :black # reaction identifier text color
reaction_edge_colors = Dict{String, Any}() # reaction id => color mapping
reaction_edge_color = :black # fallback color in case reaction id not present in reaction_edge_colors
reaction_edge_widths = Dict{String,Any}() # reaction id => edge size
reaction_edge_width = 2.0 # fallback width in case reaction id not present in reaction_edge_widths
```
Note, if `reaction_edge_colors` or `reaction_edge_widths` are supplied but missing an id
that is present in the map, the associated edge will be dotted.
