@doc """
If a matrix is structure only, it's basically an unweight matrix, whose non-zero values 
may or may not affact computation result.
"""
type TransposeProperty <: MatrixPropertyPass 

    @doc "pass name"
    name                :: AbstractString

    @doc """ set_property_for method"""
    set_property_for    :: Function

    const prop_ctranspose_pattern = ExprPattern(
        "prop_ctranspose_pattern",
        (:call, GlobalRef(Main, :ctranspose), SparseMatrixCSC),
        (:NO_SUB_PATTERNS,),
        do_nothing,
        (:NO_CHANGE, ),
        prop_propagate_first_arg,
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
        prop_ctranspose_pattern,
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
            property_map[k] = mat_property[k].transpose_of
        end

        #inherit_property(property_map, mat_property, DEFAULT_PROP_VAL, NEG_PROP_VAL, nothing)

        dprintln(1, 1, "\nBefore transpose_of analysis:")
        dprint_property_map(2, property_map)

        cnt = propagate_property(property_map, region_info, 
                    prop_propagation_patterns, nothing)

        dprintln(1, 1, "\nAfter transpose_of analysis (", cnt, " iterations):")
        dprint_property_map(2, property_map)
    end 

    @doc """ Constructor """
    function TransposeProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "transpose_of"
        return instance 
    end
end
