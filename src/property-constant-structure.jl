@doc """ Set is_constant_structured perperty for mat_property in a region """
type ConstantStructureProperty <: MatrixPropertyPass 

    @doc """ set_property_for method"""
    set_property_for    :: Function
    name                :: AbstractString

    @doc """ 
    Figure out the constant_structured property of all the matrices in a given region.
    """
    function set_property_for(
        region_info         :: RegionInfo,
        mat_property        :: Dict
    )
        # constant value is a subset of constant structure
        for (sym, v) in mat_property
            if v.constant_valued.final_val == PROP_POSITIVE_VAL  
                v.constant_structured.final_val = PROP_POSITIVE_VAL
            end
        end

        property_map = new_sym_property_map()

        # property_map[k] and mat_property[k] share values
        for k in keys(mat_property)
            property_map[k] = mat_property[k].constant_structured
            if mat_property[k].constant_valued.final_val == PROP_POSITIVE_VAL
                property_map[k].predefined = true
                property_map[k].final_val = PROP_POSITIVE_VAL
                push!(property_map[k].vals, PROP_POSITIVE_VAL)
            end
        end

        for k in keys(region_info.depend_map)
            # always prefer predefined values
            if !mat_property[k].constant_structured.predefined &&
                !isempty(region_info.depend_map[k]) &&
                !in(k, region_info.single_defs)
                    property_map[k].final_val = PROP_NEGATIVE_VAL
            end
        end

        for rk in keys(region_info.reverse_depend_map)
            if !haskey(region_info.depend_map, rk) # rk is not in def set
                property_map[rk].final_val = PROP_POSITIVE_VAL
            end
        end

        #for c in constants
        #    if !haskey(property_map, c) || property_map[c] <=0
        #        property_map[c] = 2
        #    end
        #end

        dprint_property_map(2, property_map)

        # propagate non-constant property 
        # 
        working_set = []
        for (k, v) in property_map
            if v.final_val == PROP_NEGATIVE_VAL
                push!(working_set, k)
            end
        end

        while !isempty(working_set) 
            s = shift!(working_set)
            if !haskey(region_info.reverse_depend_map, s)
                continue
            end
            for d in region_info.reverse_depend_map[s]
                if property_map[d].predefined
                    # always favor annotation
                    dprintln(1, 1, "WW fail to overwrites predefined ", d)
                elseif property_map[d].final_val != PROP_NEGATIVE_VAL
                    property_map[d].final_val = PROP_NEGATIVE_VAL
                    push!(working_set, d)
                end
            end
        end

        dprintln(1, 1, "\nAfter non-constant propagation:")
        dprint_property_map(2, property_map)


        # propagate constant property 
        working_set = []
        for (k, v) in property_map
            if v.final_val == PROP_POSITIVE_VAL
                push!(working_set, k)
            end
        end

        while !isempty(working_set) 
            s = shift!(working_set)
            if !haskey(region_info.reverse_depend_map, s)
                continue
            end
            # check every symbol that depends on s
            # 
            for rd in region_info.reverse_depend_map[s]
                if haskey(property_map, rd) && property_map[rd].final_val != nothing
                    continue
                end
                constant = true
                for d in region_info.depend_map[rd]
                    if haskey(property_map, d) && property_map[d].final_val == nothing 
                        constant = false
                    end
                end
                if constant
                    property_map[rd].final_val = PROP_POSITIVE_VAL
                    push!(working_set, rd)
                end
            end
        end

        dprintln(1, 1, "\nAfter constant propagation:")
        dprint_property_map(2, property_map)

        # copy back result
        #for (k, v) in property_map
        #    #if any(t -> symbol_info[k]<:t, MATRIX_RELATED_TYPES) 
        #        if !haskey(mat_property, k)
        #            mat_property[k] = MatrixPropertyValues()
        #        end
        #        mat_property[k].constant_structured = v
        #    #else
        #    #    dprintln(1, 1, "WW skip ", k, " ", symbol_info[k])
        #    #end
        #end
    end 

    @doc """ Constructor """
    function ConstantStructureProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "constant_structured"
        return instance 
    end
end
