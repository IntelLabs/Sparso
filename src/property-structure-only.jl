@doc """
If a matrix is structure only, it's basically an unweight matrix, whose non-zero values 
may or may not affact computation result.
"""
type StructureOnlyProperty <: MatrixPropertyPass 

    @doc "pass name"
    name                :: AbstractString

    @doc """ set_property_for method"""
    set_property_for    :: Function

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
        (:call, GlobalRef(Main, :speye), Any),
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
        property_map = new_sym_property_map()

        for k in keys(mat_property)
            property_map[k] = mat_property[k].structure_only
        end

        dprintln(1, 1, "\nBefore structure_only analysis:")
        dprint_property_map(2, property_map)

        cnt = propagate_property(property_map, region_info, 
                        prop_propagation_patterns, nothing)

        dprintln(1, 1, "\nAfter structure_only analysis (", cnt, " iterations):")
        dprint_property_map(2, property_map)
    end 

    @doc """ Constructor """
    function StructureOnlyProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "structure_only"
        return instance 
    end
end
