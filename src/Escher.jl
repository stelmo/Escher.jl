module Escher

using JSON, Makie

_bezier(t, p0, p1, p2, p3) =
    (1 - t)^3 .* p0 + 3 * (1 - t)^2 * t .* p1 + 3 * (1 - t) * t^2 .* p2 + t^3 .* p3

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

@recipe(EscherPlot) do scene
    Attributes(
        metabolite_show_text = false,
        metabolite_text_size = 4,
        metabolite_primary_node_size = 5,
        metabolite_secondary_node_size = 3,
        metabolite_node_color = :black,
        metabolite_text_color = :black,
        metabolite_identifier = "bigg_id",
        reaction_show_text = false,
        reaction_text_size = 4,
        reaction_identifier = "bigg_id",
        reaction_edge_width = 2.0,
        reaction_edge_color_weights = Dict{String,Float64}(),
        reaction_edge_width_weights = Dict{String,Float64}(),
    )
end

function Makie.plot!(ep::EscherPlot{<:Tuple{String}})

    _, escher = JSON.parsefile(to_value(ep[1]))

    reaction_identifier = to_value(ep.reaction_identifier)

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
    )

    lw = 0.0
    for (rid, rxn) in escher["reactions"]
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
                if lw == -1.0
                    lines!(ep, xs, ys, linewidth = ep.reaction_edge_width, linestyle = :dot)
                else
                    lines!(ep, xs, ys, linewidth = ep.reaction_edge_width)
                end
            else
                n1 = metabolite_node_pos[segment["from_node_id"]]
                n2 = metabolite_node_pos[segment["to_node_id"]]
                if lw == -1.0
                    lines!(
                        ep,
                        [n1[1], n2[1]],
                        [n1[2], n2[2]],
                        linewidth = ep.reaction_edge_width,
                        linestyle = :dot,
                    )
                else
                    lines!(
                        ep,
                        [n1[1], n2[1]],
                        [n1[2], n2[2]],
                        linewidth = ep.reaction_edge_width,
                    )
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
            i = 1:length(metabolite_labels)
        ]
        text!(
            ep,
            metabolite_labels;
            position = positions,
            textsize = ep.metabolite_text_size,
        )
    end

    ep
end

export plot!


# function plot_metabolism(influxes, map_path;
#     transparent_background=false, colorscheme=ColorSchemes.Reds_9, edge_weights=nothing)

#     # Plot reactions

#     if !isnothing(edge_weights)
#         efluxes = Dict(k => log(abs(v)) for (k, v) in edge_weights if k in ids_in_map && abs(v) > 1e-8)
#         eminflux = minimum(values(efluxes))
#         emaxflux = maximum(values(efluxes))
#     end

#     # Bezier cubic
#     fallback_lw = 2.0
#     for (rid, rxn) in escher["reactions"]

#         if isnothing(edge_weights)
#             lw = 2
#         else
#             if haskey(rxn, "bigg_id") && haskey(efluxes, rxn["bigg_id"])
#                 enum = (efluxes[rxn["bigg_id"]] - eminflux)/(emaxflux - eminflux)
#                 lw = 0.5 + 5*enum
#             else
#                 lw = 0
#             end
#         end

#         if haskey(rxn, "bigg_id") && haskey(fluxes, rxn["bigg_id"])
#             cnum = (fluxes[rxn["bigg_id"]] - minflux)/(maxflux - minflux)
#         else
#             cnum = 0.0
#             lw = 0.0
#         end

#         for (_, segment) in rxn["segments"]

#             if !isnothing(segment["b1"]) && !isnothing(segment["b2"])
#                 n1 = node_pos[segment["from_node_id"]]
#                 p0 = [n1[1], n1[2]]
#                 p1 = [segment["b1"]["x"], -segment["b1"]["y"]]
#                 p2 = [segment["b2"]["x"], -segment["b2"]["y"]]
#                 n2 = node_pos[segment["to_node_id"]]
#                 p3 = [n2[1], n2[2]]
#                 ts = range(0.0, 1.0; length=10)
#                 ds = B3.(ts, Ref(p0), Ref(p1), Ref(p2), Ref(p3))
#                 xs = [d[1] for d in ds]
#                 ys = [d[2] for d in ds]
#                 if lw == 0
#                     lines!(ax, xs, ys, linewidth=fallback_lw, color=get(colorscheme, cnum), linestyle=:dot)
#                 else
#                     lines!(ax, xs, ys, linewidth=lw, color=get(colorscheme, cnum))
#                 end
#             else
#                 n1 = node_pos[segment["from_node_id"]]
#                 n2 = node_pos[segment["to_node_id"]]
#                 if lw == 0
#                     lines!(ax, [n1[1], n2[1]], [n1[2], n2[2]], linewidth=fallback_lw, color=get(colorscheme, cnum), linestyle=:dot)
#                 else
#                     lines!(ax, [n1[1], n2[1]], [n1[2], n2[2]], linewidth=lw, color=get(colorscheme, cnum))
#                 end
#             end
#         end
#     end

#     hidexdecorations!(ax)
#     hideydecorations!(ax)
#     Colorbar(fig[2,1], limits=(0, 1), colormap=colorscheme, label="Flux [AU]", vertical = false)

#     return fig
# end


end # module
