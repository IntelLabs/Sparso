@doc """
If a matrix is structure only, it's basically an unweight matrix, whose non-zero values 
may or may not affact computation result.
"""
type StructureOnlyProperty <: MatrixProperty 

    @doc "pass name"
    name                :: AbstractString

    @doc """ set_property_for method"""
    set_property_for    :: Function

    const DEFAULT_PROP_VAL = 0
    const NEG_PROP_VAL = -1

    const prop_spones_pattern = ExprPattern(
        "prop_spones_pattern",
        (:call, GlobalRef(Main, :spones), SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        prop_set_output_action,
        "",
        "",
        (),
        0,
        ()
    )

    const prop_speye_pattern = ExprPattern(
        "prop_speye_pattern",
        (:call, GlobalRef(Main, :speye), Any, Any),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        prop_set_output_action,
        "",
        "",
        (),
        0,
        ()
    )

    @doc """
    Patterns used for discovering matrix structures. 
    """
    prop_propagation_patterns = [
        prop_spones_pattern,
        prop_speye_pattern,
        prop_assign_pattern,
        prop_assign2_pattern,
        prop_apply_type_pattern,
        prop_last_resort_pattern
    ]

    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        region_info         :: RegionInfo,
        mat_property        :: Dict,
    )
        property_map = new_sym_property_map(Int, DEFAULT_PROP_VAL, NEG_PROP_VAL)

        println(region_info.single_defs)

        for k in keys(region_info.depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) 
                if mat_property[k].structure_only != DEFAULT_PROP_VAL
                    property_map[k] = mat_property[k].structure_only
                elseif mat_property[k].constant_structured < 0
                    property_map[k] = NEG_PROP_VAL
                end
            else
                property_map[k] = in(k, region_info.single_defs) ? DEFAULT_PROP_VAL : NEG_PROP_VAL
            end
        end

        for rk in keys(region_info.reverse_depend_map)
            if !haskey(property_map, rk)  && 
                        haskey(mat_property, rk) && mat_property[rk].structure_only != DEFAULT_PROP_VAL
                property_map[rk] = mat_property[rk].structure_only
            end
        end

        #inherit_property(property_map, mat_property, DEFAULT_PROP_VAL, NEG_PROP_VAL, nothing)

        dprintln(1, 1, "\nBefore structure_only analysis:")
        dprint_property_map(1, property_map)

        cnt = propagate_property(property_map, region_info, prop_propagation_patterns, 
            Int, DEFAULT_PROP_VAL, NEG_PROP_VAL,
            nothing,
        )

        dprintln(1, 1, "\nAfter structure_only analysis (", cnt, " iterations):")
        dprint_property_map(1, property_map)

        # copy back result
        # delete!(property_map, :NEGATIVE_PROPERTY)
        for (k, v) in property_map
            if !haskey(mat_property, k)
                mat_property[k] = StructureProxy()
            end
            mat_property[k].structure_only = v
        end
    end 

    @doc """ Constructor """
    function StructureOnlyProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "structure_only"
        return instance 
    end
end
