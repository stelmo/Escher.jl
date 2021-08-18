function plot_metabolism(influxes, map_path;
    transparent_background=false, colorscheme=ColorSchemes.Reds_9, edge_weights=nothing)

    _, escher = JSON.parsefile(map_path)

    # Plot metabolites
    if transparent_background
        fig = Figure(backgroundcolor=:transparent);
    else
        fig = Figure();
    end
    ax = Axis(fig[1,1]);

    nodexs = Float64[]
    nodeys = Float64[]
    node_pos = Dict()
    markersizes = Float64[]
    for (node_id, node) in escher["nodes"]
        if haskey(node, "x") && haskey(node, "y")
            x = node["x"]
            y = -node["y"]
            node_pos[node_id] = (x, y)
            # plotted nodes
            if node["node_type"] == "metabolite"
                push!(nodexs, x)
                push!(nodeys, y)
                if node["node_is_primary"]
                    push!(markersizes, 5)
                else
                    push!(markersizes, 3)
                end
            end
        end
    end

    # Plot reactions
    ids_in_map = [rxn["bigg_id"] for rxn in values(escher["reactions"]) if haskey(rxn, "bigg_id")]

    fluxes = Dict(k => log(abs(v)) for (k, v) in influxes if k in ids_in_map && abs(v) > 1e-8)
    minflux = minimum(values(fluxes))
    maxflux = maximum(values(fluxes))

    if !isnothing(edge_weights)
        efluxes = Dict(k => log(abs(v)) for (k, v) in edge_weights if k in ids_in_map && abs(v) > 1e-8)
        eminflux = minimum(values(efluxes))
        emaxflux = maximum(values(efluxes))
    end

    # Bezier cubic
    B3(t, p0, p1, p2, p3) = (1-t)^3 .* p0  + 3 * (1-t)^2 * t .* p1 + 3 * (1-t) * t^2 .* p2 + t^3 .* p3
    fallback_lw = 2.0
    for (rid, rxn) in escher["reactions"]

        if isnothing(edge_weights)
            lw = 2
        else
            if haskey(rxn, "bigg_id") && haskey(efluxes, rxn["bigg_id"])
                enum = (efluxes[rxn["bigg_id"]] - eminflux)/(emaxflux - eminflux)
                lw = 0.5 + 5*enum
            else
                lw = 0
            end
        end

        if haskey(rxn, "bigg_id") && haskey(fluxes, rxn["bigg_id"])
            cnum = (fluxes[rxn["bigg_id"]] - minflux)/(maxflux - minflux)
        else
            cnum = 0.0
            lw = 0.0
        end

        for (_, segment) in rxn["segments"]

            if !isnothing(segment["b1"]) && !isnothing(segment["b2"])
                n1 = node_pos[segment["from_node_id"]]
                p0 = [n1[1], n1[2]]
                p1 = [segment["b1"]["x"], -segment["b1"]["y"]]
                p2 = [segment["b2"]["x"], -segment["b2"]["y"]]
                n2 = node_pos[segment["to_node_id"]]
                p3 = [n2[1], n2[2]]
                ts = range(0.0, 1.0; length=10)
                ds = B3.(ts, Ref(p0), Ref(p1), Ref(p2), Ref(p3))
                xs = [d[1] for d in ds]
                ys = [d[2] for d in ds]
                if lw == 0
                    lines!(ax, xs, ys, linewidth=fallback_lw, color=get(colorscheme, cnum), linestyle=:dot)
                else
                    lines!(ax, xs, ys, linewidth=lw, color=get(colorscheme, cnum))
                end
            else
                n1 = node_pos[segment["from_node_id"]]
                n2 = node_pos[segment["to_node_id"]]
                if lw == 0
                    lines!(ax, [n1[1], n2[1]], [n1[2], n2[2]], linewidth=fallback_lw, color=get(colorscheme, cnum), linestyle=:dot)
                else
                    lines!(ax, [n1[1], n2[1]], [n1[2], n2[2]], linewidth=lw, color=get(colorscheme, cnum))
                end
            end
        end
    end

    scatter!(ax, nodexs, nodeys, color=:black, markersize=markersizes)
    hidexdecorations!(ax)
    hideydecorations!(ax)
    Colorbar(fig[2,1], limits=(0, 1), colormap=colorscheme, label="Flux [AU]", vertical = false)

    return fig
end
