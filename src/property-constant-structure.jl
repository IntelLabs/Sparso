@doc """ Set is_constant_structured perperty for mat_property in a region """
type ConstantStructureProperty <: MatrixProperty 

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
        # constant structure is a subset of constant value
        for (sym, v) in mat_property
            if v.constant_valued > 0 && v.constant_structured == 0
                v.constant_structured = v.constant_valued
            end
        end

        # property values: 
        #  -1: not constant 
        #  0: unknow 
        #  1: constant 
        #  2: external(constant)
        #  3: specified by set_matrix_property statement
        #  4: inherited from parent region (constant)
        property_map = new_sym_property_map(Int, 0, -1)

        for k in keys(region_info.depend_map)
            # always prefer predefined values
            if haskey(mat_property, k) 
                if mat_property[k].constant_structured!=0 
                    property_map[k] = mat_property[k].constant_structured 
                end
            else
                property_map[k] = in(k, region_info.single_defs) ? 0 : -1
            end
        end

        for rk in keys(region_info.reverse_depend_map)
            if !haskey(property_map, rk) # rk is not in def set
                property_map[rk] = 2
            end
        end

        #for c in constants
        #    if !haskey(property_map, c) || property_map[c] <=0
        #        property_map[c] = 2
        #    end
        #end

        dprint_property_map(1, property_map)

        # propagate non-constant property 
        # 
        working_set = []
        for (k, v) in property_map
            if v < 0
                push!(working_set, k)
            end
        end

        while !isempty(working_set) 
            s = shift!(working_set)
            if !haskey(region_info.reverse_depend_map, s)
                continue
            end
            for d in region_info.reverse_depend_map[s]
                if property_map[d] >= 0
                    if property_map[d] == 3
                        # always favor annotation
                        dprintln(1, 1, "WW annotation overwrites non_constant ", d)
                    elseif property_map[d] == 4
                        # always favor inherit result
                        dprintln(1, 1, "WW inherit overwrites non_constant ", d)
                    else
                        property_map[d] = -1
                        push!(working_set, d)
                    end
                end
            end
        end

        dprintln(1, 1, "\nAfter non-constant propagation:")
        dprint_property_map(1, property_map)


        # propagate constant property 
        working_set = []
        for (k, v) in property_map
            if v > 0
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
                if haskey(property_map, rd) && property_map[rd] != 0
                    continue
                end
                constant = true
                for d in region_info.depend_map[rd]
                    if haskey(property_map, d) && property_map[d] == 0 
                        constant = false
                    end
                end
                if constant
                    property_map[rd] = 1
                    push!(working_set, rd)
                end
            end
        end

        dprintln(1, 1, "\nAfter constant propagation:")
        dprint_property_map(1, property_map)

        # copy back result
        #delete!(property_map, sym_default_property_key)
        #delete!(property_map, sym_negative_property_key)
        for (k, v) in property_map
            #if any(t -> symbol_info[k]<:t, MATRIX_RELATED_TYPES) 
                if !haskey(mat_property, k)
                    mat_property[k] = StructureProxy()
                end
                mat_property[k].constant_structured = v
            #else
            #    dprintln(1, 1, "WW skip ", k, " ", symbol_info[k])
            #end
        end
    end 

    @doc """ Constructor """
    function ConstantStructureProperty()
        instance = new()
        instance.set_property_for = set_property_for
        instance.name = "constant_structured"
        return instance 
    end
end
