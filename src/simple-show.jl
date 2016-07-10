import Base.show_expr_type, Base.show_linenumber, Base.show_call
# Re-define the function. Otherwise, the show of an expression may be too verbose
function show_expr_type(io::IO, ty)
    return 
end

function show_linenumber(io::IO, file, line)
    return
end

