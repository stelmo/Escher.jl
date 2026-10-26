module Escher

using JSON, Makie, DocStringExtensions

"""
$(TYPEDSIGNATURES)

Return a named tuple with the dimensions of the escher map.
Fields are `height`, `width`, `x`, and `y`. If the fields are 
missing in the supplied map, return `missing`. 

# Example
```
h,w,x,y = get_resolution(map_location)

f = Figure(resolution = (x, y))
```
"""
function get_resolution(escher_location::String)
    _, escher = JSON.parsefile(escher_location)
    (
        height = height = get(escher["canvas"], "height", missing),
        width = get(escher["canvas"], "width", missing),
        x = get(escher["canvas"], "x", missing),
        y = get(escher["canvas"], "y", missing),
    )
end

"""
$(TYPEDSIGNATURES)

Plot a Bezier curve with `t` in `[0,1]` and fixed points `p0`, `p1`, `p2`, `p3`.
"""
_bezier(t, p0, p1, p2, p3) =
    (1 - t)^3 .* p0 + 3 * (1 - t)^2 * t .* p1 + 3 * (1 - t) * t^2 .* p2 + t^3 .* p3

"""
$(TYPEDSIGNATURES)

Helper function to extract node information from `escher`. The kwargs set plot
details.
"""
function _nodes(
    escher;
    primary_node_size = 5,
    secondary_node_size = 3,
    metabolite_identifier = "bigg_id",
    metabolite_node_sizes = Dict{String,Any}(),
    metabolite_node_colors = Dict{String,Any}(),
    metabolite_node_color = :black,
)
    nodexs = Float64[]
    nodeys = Float64[]
    labelxs = Float64[]
    labelys = Float64[]
    left_right_aligns = Symbol[]
    top_bottom_aligns = Symbol[]
    justifications = Symbol[]
    metabolite_labels = String[]
    node_pos_lookup = Dict()
    markersizes = Float64[]
    markercolors = []
    metabolite_ids = String[]
    for (node_id, node) in escher["nodes"] #  cycle through all nodes 
        if haskey(node, "x") && haskey(node, "y") # must have location to be meaningful
            #=
            In an escher map, nodes can be either metabolites or markers. Metabolites 
            should get plotted as dots or circles, but markers are used to join reaction 
            line segments to each other. Use node_pos_lookup to map the internal id of 
            nodes to their (x,y) positions if they are markers. Otherwise, the nodes 
            represent metabolites, in which case more plotting information is stripped 
            from the file.
            =#
            x = node["x"]
            y = -node["y"]
            node_pos_lookup[node_id] = (x, y)
            #=
            Strip all the relevant plotting information from the escher file if the 
            node type is a metabolite, and thus should be plotted.
            =#
            if node["node_type"] == "metabolite"
                push!(metabolite_ids, node_id) # internal id 
                push!(nodexs, x)
                push!(nodeys, y)
                push!(
                    markercolors,
                    get(
                        metabolite_node_colors,
                        node[metabolite_identifier],
                        metabolite_node_color,
                    ),
                )
                push!(
                    markersizes,
                    get(
                        metabolite_node_sizes,
                        node[metabolite_identifier],
                        node["node_is_primary"] ? primary_node_size : secondary_node_size,
                    ),
                )
                # add text label of metabolite if it exists
                if haskey(node, metabolite_identifier)
                    nx = node["label_x"]
                    ny = -node["label_y"]
                    push!(metabolite_labels, node[metabolite_identifier])
                    push!(labelxs, nx)
                    push!(labelys, ny)
                    push!(
                        left_right_aligns, 
                        nx < x ? :right : :left, 
                    )
                    push!(
                        justifications, 
                        nx < x ? :left : :right,
                    )
                    push!(
                        top_bottom_aligns,
                        ny < y ? :top : :bottom,
                    )
                end
            end
        end
    end
    metabolites = (
        ids = metabolite_ids, # internal ids
        positions = collect(zip(nodexs, nodeys)),
        colors = markercolors,
        label_positions = collect(zip(labelxs, labelys)),
        alignments = collect(zip(left_right_aligns, top_bottom_aligns)),
        justifications = justifications,
        markersizes = markersizes,
        labels = metabolite_labels,
    )
    return metabolites, node_pos_lookup
end

"""
$(TYPEDSIGNATURES)

Return an named tuple of text annotation mappings.
The fields are `positions`, which maps `x, y` coordinates 
to the field `labels`.
"""
function _text_annotations(escher)
    xs = Float64[]
    ys = Float64[]
    labels = String[]
    for (_, node) in escher["text_labels"]
        push!(xs, get(node, "x", 0.0))
        push!(ys, -get(node, "y", 0.0))
        push!(labels, get(node, "text", ""))
    end
    return (positions = collect(zip(xs, ys)), labels = labels)
end

"""
$(TYPEDSIGNATURES)

A helper function to get edge information from `escher`. Kwargs set 
plot details.
"""
function _edges(
    escher,
    node_pos_lookup,
    metabolite_id_lookup;
    reaction_identifier = "bigg_id", # used in mapping and possibly to display
    reaction_show_name_instead_of_id = false, # use name field in escher, not reaction_identifier flag
    reaction_edge_colors = Dict{String,Any}(), # actual color
    reaction_edge_color = :black, # fallback color
    reaction_edge_widths = Dict{String,Any}(), # actual edge width
    reaction_edge_width = 2.0, # fallback width
)
    # reaction label data
    reaction_labels = String[]
    reaction_label_xs = Float64[]
    reaction_label_ys = Float64[]

    #=  
    A mapping between the human readable reaction id, and the 
    internal id used for each reaction in an escher map.
    =#
    rid_id_lookup = Dict(
        rid => rxn[to_value(reaction_identifier)] for (rid, rxn) in escher["reactions"]
    )

    #=
    A array holding named tuples. Each tuple contains the basic information 
    necessary to plot the reaction.  
    =#
    reactions = []
    for (id, rxn) in escher["reactions"]
        # collect reaction label information
        _reaction_identifier = reaction_show_name_instead_of_id ? "name" : reaction_identifier
        if haskey(rxn, _reaction_identifier)
            push!(reaction_labels, rxn[_reaction_identifier])
            push!(reaction_label_xs, rxn["label_x"])
            push!(reaction_label_ys, -rxn["label_y"])
        end

        rid = rid_id_lookup[id]

        lw = get(reaction_edge_widths, rid, reaction_edge_width)
        c = get(reaction_edge_colors, rid, reaction_edge_color)

        #=
        If reaction data is NOT supplied by the user, all reactions are plotted.
        If reaction data is supplied by the user, dash the lines of reactions missing from these datasets.
        In making :dot looks better than :dash, so use :dot.
        =#
        if isempty(reaction_edge_widths) && isempty(reaction_edge_colors) # no user data
            dot_reaction_line = false
        elseif haskey(reaction_edge_widths, rid) || haskey(reaction_edge_colors, rid) # data supplied, test if rxn has data
            dot_reaction_line = false
        else
            dot_reaction_line = true # data supplied, but no rxn data
        end

        #=
        An array of named tuples, where each tuple contains the building blocks to plot 
        one reaction. The building blocks are the reaction segments in the escher map. 
        The to/from_metabolite field is used to determine the direction of the reaction.
        =#
        segments = []
        for (_, segment) in rxn["segments"]
            n1 = node_pos_lookup[segment["from_node_id"]]
            n2 = node_pos_lookup[segment["to_node_id"]]
            to_metabolite = get(metabolite_id_lookup, segment["to_node_id"], nothing)
            from_metabolite = get(metabolite_id_lookup, segment["from_node_id"], nothing)
            if !isnothing(segment["b1"]) && !isnothing(segment["b2"])
                p0 = [n1[1], n1[2]]
                p1 = [segment["b1"]["x"], -segment["b1"]["y"]]
                p2 = [segment["b2"]["x"], -segment["b2"]["y"]]
                p3 = [n2[1], n2[2]]
                ts = range(0.0, 1.0; length = 200) # higher resolution => looks smoother
                ds = _bezier.(ts, Ref(p0), Ref(p1), Ref(p2), Ref(p3))
                xs = [d[1] for d in ds]
                ys = [d[2] for d in ds]
            else
                xs = range(n1[1], n2[1]; length = 200) # higher resolution => looks smoother
                ys = range(n1[2], n2[2]; length = 200) # higher resolution => looks smoother
            end
            push!(
                segments,
                (
                    positions = collect(zip(xs, ys)),
                    to_metabolite = to_metabolite,
                    from_metabolite = from_metabolite,
                ),
            )
        end

        push!(
            reactions,
            (
                color = c,
                linewidth = lw,
                rid = rid,
                segments = segments,
                dot_reaction_line = dot_reaction_line,
            ),
        )
    end

    reaction_labels = (
        positions = collect(zip(reaction_label_xs, reaction_label_ys)),
        labels = reaction_labels,
    )

    return reaction_labels, reactions
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
reaction_edge_width = 1.0 # fallback width
reaction_arrow_scale_size = 10.0
reaction_arrow_head_offset_fraction = 0.2 # between 0 and 1
reaction_directions = Dict{String,Tuple{Dict{String,Number},Symbol}}() # rid => (reaction stoichiometry, :forward, :backward, :bidirectional)
annotation_show_text = false
annotation_text_color = :black
annotation_text_size = 12
```
Get or create maps here: `https://escher.github.io/#/`.
"""
@recipe(EscherPlot) do scene
    Attributes(
        metabolite_identifier = "bigg_id",
        metabolite_show_text = false,
        metabolite_text_size = 4,
        metabolite_primary_node_size = 5, # fallback size
        metabolite_secondary_node_size = 3, # fallback size
        metabolite_node_sizes = Dict{String,Any}(),
        metabolite_node_colors = Dict{String,Any}(),
        metabolite_node_color = :black, # fallback color
        metabolite_text_color = :black,
        reaction_identifier = "bigg_id",
        reaction_show_text = false,
        reaction_show_name_instead_of_id = false,
        reaction_text_size = 4,
        reaction_text_color = :black,
        reaction_edge_colors = Dict{String,Any}(), # actual color
        reaction_edge_color = :black, # fallback color
        reaction_edge_widths = Dict{String,Any}(), # actual edge width
        reaction_edge_width = 1.0, # fallback width
        reaction_arrow_scale_size = 10.0,
        reaction_arrow_head_offset_fraction = 0.2, # between 0 and 1
        reaction_directions = Dict{String,Tuple{Dict{String,Number},Symbol}}(), # rid => (reaction stoichiometry, :forward, :backward, :bidirectional)
        annotation_show_text = false,
        annotation_text_color = :black,
        annotation_text_size = 12,
    )
end

function Makie.plot!(ep::EscherPlot{<:Tuple{String}})

    _, escher = JSON.parsefile(to_value(ep[1]))

    annotations = _text_annotations(escher)

    metabolites, node_pos_lookup = _nodes(
        escher;
        primary_node_size = to_value(ep.metabolite_primary_node_size),
        secondary_node_size = to_value(ep.metabolite_secondary_node_size),
        metabolite_identifier = to_value(ep.metabolite_identifier),
        metabolite_node_sizes = to_value(ep.metabolite_node_sizes),
        metabolite_node_colors = to_value(ep.metabolite_node_colors),
        metabolite_node_color = to_value(ep.metabolite_node_color),
    )

    # internal id to metabolite id lookup, used for direction
    metabolite_id_lookup = Dict(zip(metabolites.ids, metabolites.labels))

    reaction_labels, reactions = _edges(
        escher,
        node_pos_lookup,
        metabolite_id_lookup;
        reaction_identifier = to_value(ep.reaction_identifier),
        reaction_show_name_instead_of_id = to_value(ep.reaction_show_name_instead_of_id),
        reaction_edge_colors = to_value(ep.reaction_edge_colors),
        reaction_edge_color = to_value(ep.reaction_edge_color),
        reaction_edge_widths = to_value(ep.reaction_edge_widths),
        reaction_edge_width = to_value(ep.reaction_edge_width),
    )

    # Plot reactions
    reaction_directions = to_value(ep.reaction_directions)
    reaction_arrow_scale_size = to_value(ep.reaction_arrow_scale_size)
    reaction_arrow_head_offset_fraction = to_value(ep.reaction_arrow_head_offset_fraction)

    # Default Makie arrow tip lengths/widths
    tip_length = 8
    tip_width = 14

    for reaction in reactions
        rid = reaction.rid
        color = reaction.color
        linewidth = reaction.linewidth
        linestyle = reaction.dot_reaction_line ? :dot : :solid
        # Plot reaction backbone
        for segment in reaction.segments
            lines!(ep, segment.positions; linewidth, linestyle, color)
        end
        #=
        If reaction has direction information supplied, then add that to the
        plot. The convention is to only add arrows to reactions if the user
        specifies the direction. Moreover, arrows only point into the metabolite
        that is produced by the reaction in the direction of flow. No arrows are
        plotted going from the source metabolite. Lastly, only the arrow head is 
        plotted, this makes it easier to control the the plot.
        =#
        if haskey(reaction_directions, rid)
            rs = first(reaction_directions[rid])
            dir = last(reaction_directions[rid])
            if dir == :forward
                target_metabolites = [mid for (mid, x) in rs if x > 0]
            elseif dir == :backward 
                target_metabolites = [mid for (mid, x) in rs if x < 0]
            else # bidirectional is default
                target_metabolites = [mid for (mid, _) in rs]
            end
            for segment in reaction.segments
                to_metabolite = segment.to_metabolite
                from_metabolite = segment.from_metabolite

                #=
                Find index that corresponds to reaction_arrow_head_offset_fraction 
                of the curve length. This is necessary because the scale for curves 
                is not linear, and makes arrow location tricky.
                =#
                if to_metabolite in metabolites.labels &&
                   to_metabolite in target_metabolites
                    reaction_arrow_head_offset = _curve_idx(segment, reaction_arrow_head_offset_fraction)
                    offset_idx = length(segment.positions) - reaction_arrow_head_offset
                    offset_idx2 = offset_idx + 1
                elseif from_metabolite in metabolites.labels &&
                       from_metabolite in target_metabolites
                    reaction_arrow_head_offset = _curve_idx(segment, 1-reaction_arrow_head_offset_fraction)
                    offset_idx = length(segment.positions) - reaction_arrow_head_offset
                    offset_idx2 = offset_idx - 1
                else
                    continue
                end

                points = [Point2f(segment.positions[offset_idx]...)]
                directions = [
                    Point2f(
                        segment.positions[offset_idx2][1] -
                        segment.positions[offset_idx][1],
                        segment.positions[offset_idx2][2] -
                        segment.positions[offset_idx][2],
                    ),
                ]

                #arrows!(
                    #ep,
                    #points,
                    #directions;
                    #arrowcolor = color,
                    #linewidth = 0.0,
                    #arrowsize = reaction_arrow_size,
                    #linestyle,
                #)
                arrows2d!(
                    ep,
                    points,
                    directions;
                    color = color,
                    tiplength = tip_length * reaction_arrow_scale_size,
                    tipwidth = tip_width * reaction_arrow_scale_size,
                    lengthscale = 20,
                    shaftwidth = 0.0
                )
            end
        end
    end

    # Plot reaction labels
    if ep.reaction_show_text[]
        text!(
            ep,
            reaction_labels.labels;
            position = reaction_labels.positions,
            fontsize = ep.reaction_text_size,
            color = ep.reaction_text_color,
            # align = (:center, :center),
            # justification = :left,
        )
    end

    # Plot metabolites 
    scatter!(
        ep,
        metabolites.positions,
        color = metabolites.colors,
        markersize = metabolites.markersizes,
    )

    # Plot metabolite labels
    if ep.metabolite_show_text[]
        met_label_list = zip(
            metabolites.labels,
            metabolites.label_positions,
            metabolites.alignments,
            metabolites.justifications
        )
        for (label, position, alignment, justification) in met_label_list
            text!(
                ep,
                label;
                position = position,
                fontsize = ep.metabolite_text_size,
                color = ep.metabolite_text_color,
                # align = alignment,
                # justification = justification,
            )
        end
    end

    # Plot annotations
    if ep.annotation_show_text[]
        text!(
            ep,
            annotations.labels;
            position = annotations.positions,
            fontsize = ep.annotation_text_size,
            color = ep.annotation_text_color,
        )
    end

    ep
end

"""
$(TYPEDSIGNATURES)

Helper function to approximately measure the length of a curve contained in
`segment`, and return the index of the first point that is `fraction` larger
than the total length
"""
function _curve_idx(segment, fraction)
    curve_length = cumsum([
        _euclidian_dist(p1, p2) for
        (p1, p2) in zip(segment.positions[1:end-1], segment.positions[2:end])
    ])
    findfirst(curve_length .> last(curve_length)*fraction)
end

"""
$(TYPEDSIGNATURES)

A helper function to calculate the Euclidian distance between 
points `p1` and `p2`.
"""
_euclidian_dist(p1, p2) = sqrt(sum((x - y)^2 for (x, y) in zip(p1, p2)))

export plot!, get_resolution

end # module
