module Escher

using JSON, Makie

"""
    _bezier(t, p0, p1, p2, p3)

Plot a Bezier curve with `t` in `[0,1]` and fixed points `p0`, `p1`, `p2`, `p3`.
"""
_bezier(t, p0, p1, p2, p3) =
    (1 - t)^3 .* p0 + 3 * (1 - t)^2 * t .* p1 + 3 * (1 - t) * t^2 .* p2 + t^3 .* p3

"""
    _nodes(
        escher;
        primary_node_size = 5,
        secondary_node_size = 3,
        metabolite_identifier = "bigg_id",
    )

Helper function to extract metabolite node information from `escher`. The kwargs set plot
details.
"""
function _nodes(
    escher;
    primary_node_size = 5,
    secondary_node_size = 3,
    metabolite_identifier = "bigg_id",
)
    nodexs = Float64[]
    nodeys = Float64[]
    labelxs = Float64[]
    labelys = Float64[]
    metabolite_labels = String[]
    node_pos = Dict()
    markersizes = Float64[]
    for (node_id, node) in escher["nodes"]
        if haskey(node, "x") && haskey(node, "y")
            x = node["x"]
            y = -node["y"]
            node_pos[node_id] = (x, y)
            # plotted nodes
            if node["node_type"] == "metabolite"
                if haskey(node, metabolite_identifier)
                    push!(metabolite_labels, node[metabolite_identifier])
                    push!(labelxs, node["label_x"])
                    push!(labelys, -node["label_y"])
                end
                push!(nodexs, x)
                push!(nodeys, y)
                if node["node_is_primary"]
                    push!(markersizes, primary_node_size)
                else
                    push!(markersizes, secondary_node_size)
                end
            end
        end
    end
    return nodexs, nodeys, node_pos, markersizes, labelxs, labelys, metabolite_labels
end

"""
Plot a metabolic map that is compatible with Escher. The only required argument is the
location of the metabolic map in `json` format.

# Example
```
escherplot("core-map.json"; kwargs...)
```
Here `kwargs` are supported attributes, see the `readme` for more information.

# Attributes
```
metabolite_identifier = "bigg_id"
metabolite_show_text = false
metabolite_text_size = 4
metabolite_primary_node_size = 5
metabolite_secondary_node_size = 3
metabolite_node_color = :black
metabolite_text_color = :black
reaction_identifier = "bigg_id"
reaction_show_text = false
reaction_show_name_instead_of_id = false
reaction_text_size = 4
reaction_text_color = :black
reaction_edge_colors = Dict{String,Any}()
reaction_edge_color = :black
reaction_edge_widths = Dict{String,Any}()
reaction_edge_width = 2.0
```
Get or create maps here: `https://escher.github.io/#/`.
"""
@recipe(EscherPlot) do scene
    Attributes(
        metabolite_identifier = "bigg_id",
        metabolite_show_text = false,
        metabolite_text_size = 4,
        metabolite_primary_node_size = 5,
        metabolite_secondary_node_size = 3,
        metabolite_node_color = :black, # could be a vector, but order not fixed TODO
        metabolite_text_color = :black, # could be a vector, but order not fixed TODO
        reaction_identifier = "bigg_id",
        reaction_show_text = false,
        reaction_show_name_instead_of_id = false,
        reaction_text_size = 4,
        reaction_text_color = :black, # could be a vector, but order not fixed TODO
        reaction_edge_colors = Dict{String,Any}(), # actual color
        reaction_edge_color = :black, # fallback color
        reaction_edge_widths = Dict{String,Any}(), # actual edge width
        reaction_edge_width = 2.0, # fallback width
    )
end

function Makie.plot!(ep::EscherPlot{<:Tuple{String}})

    _, escher = JSON.parsefile(to_value(ep[1]))

    metabolite_node_xs,
    metabolite_node_ys,
    metabolite_node_pos,
    metabolite_markersizes,
    metabolite_label_xs,
    metabolite_label_ys,
    metabolite_labels = _nodes(
        escher;
        primary_node_size = to_value(ep.metabolite_primary_node_size),
        secondary_node_size = to_value(ep.metabolite_secondary_node_size),
        metabolite_identifier = to_value(ep.metabolite_identifier),
    )

    reaction_labels = String[]
    reaction_label_xs = Float64[]
    reaction_label_ys = Float64[]
    reaction_id_to_id = Dict(
        rid => rxn[to_value(ep.reaction_identifier)] for (rid, rxn) in escher["reactions"]
    )
    for (rid, rxn) in escher["reactions"]
        # collect reaction label information
        if ep.reaction_show_name_instead_of_id[]
            if haskey(rxn, "name")
                push!(reaction_labels, rxn["name"])
                push!(reaction_label_xs, rxn["label_x"])
                push!(reaction_label_ys, -rxn["label_y"])
            end
        else
            if haskey(rxn, ep.reaction_identifier[])
                push!(reaction_labels, rxn[ep.reaction_identifier[]])
                push!(reaction_label_xs, rxn["label_x"])
                push!(reaction_label_ys, -rxn["label_y"])
            end
        end

        val_rid = reaction_id_to_id[rid]

        lw =
            haskey(to_value(ep.reaction_edge_widths), val_rid) ?
            to_value(ep.reaction_edge_widths)[val_rid] : ep.reaction_edge_width
        c =
            haskey(to_value(ep.reaction_edge_colors), val_rid) ?
            to_value(ep.reaction_edge_colors)[val_rid] : ep.reaction_edge_color

        if isempty(to_value(ep.reaction_edge_widths)) &&
           isempty(to_value(ep.reaction_edge_color_weights))
            no_reaction_data = false
        elseif haskey(to_value(ep.reaction_edge_widths), val_rid) ||
               haskey(to_value(ep.reaction_edge_colors), val_rid)
            no_reaction_data = false
        else
            no_reaction_data = true
        end

        # linestyle = :dot if missing information
        for (_, segment) in rxn["segments"]
            if !isnothing(segment["b1"]) && !isnothing(segment["b2"])
                n1 = metabolite_node_pos[segment["from_node_id"]]
                p0 = [n1[1], n1[2]]
                p1 = [segment["b1"]["x"], -segment["b1"]["y"]]
                p2 = [segment["b2"]["x"], -segment["b2"]["y"]]
                n2 = metabolite_node_pos[segment["to_node_id"]]
                p3 = [n2[1], n2[2]]
                ts = range(0.0, 1.0; length = 10)
                ds = _bezier.(ts, Ref(p0), Ref(p1), Ref(p2), Ref(p3))
                xs = [d[1] for d in ds]
                ys = [d[2] for d in ds]
                if no_reaction_data
                    lines!(ep, xs, ys, linewidth = lw, linestyle = :dot, color = c)
                else
                    lines!(ep, xs, ys, linewidth = lw, color = c)
                end
            else
                n1 = metabolite_node_pos[segment["from_node_id"]]
                n2 = metabolite_node_pos[segment["to_node_id"]]
                if no_reaction_data
                    lines!(
                        ep,
                        [n1[1], n2[1]],
                        [n1[2], n2[2]],
                        linewidth = lw,
                        linestyle = :dot,
                        color = c,
                    )
                else
                    lines!(ep, [n1[1], n2[1]], [n1[2], n2[2]], linewidth = lw, color = c)
                end
            end
        end
    end

    scatter!(
        ep,
        metabolite_node_xs,
        metabolite_node_ys,
        color = ep.metabolite_node_color,
        markersize = metabolite_markersizes,
    )

    if ep.metabolite_show_text[]
        positions = [
            (metabolite_label_xs[i], metabolite_label_ys[i]) for
            i in eachindex(metabolite_labels)
        ]
        text!(
            ep,
            metabolite_labels;
            position = positions,
            textsize = ep.metabolite_text_size,
            color = ep.metabolite_text_color,
        )
    end

    if ep.reaction_show_text[]
        positions = [
            (reaction_label_xs[i], reaction_label_ys[i]) for
            i in eachindex(reaction_labels)
        ]
        text!(
            ep,
            reaction_labels;
            position = positions,
            textsize = ep.reaction_text_size,
            color = ep.reaction_text_color,
        )

    end

    ep
end

export plot!

end # module
