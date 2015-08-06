# This file contains all the print parameters and routines for warnings,
# errors, and debugging.

module Base
    # Re-define the function. Otherwise, the show of an expression may be 
    # too verbose
    function show_expr_type(io::IO, ty)
        return
    end
end

# Buffer the messages so that indentation can be inserted before sending
# to STDOUT
io_buffer = IOBuffer()

function show(level :: Int, indent :: Int, new_line :: Bool, msgs...)
    if(level >= verbosity)
        message_tuple = msgs[1]
        for msg in message_tuple
            print(io_buffer, msg)
        end
        
        seekstart(io_buffer)
        while !eof(io_buffer)
            [print(" ") for i = 1 : indent * 4]
            s = readline(io_buffer)
            print_unescaped(STDOUT, s)
        end
        if new_line println(STDOUT, "") end

        # Resize the buffer to be empty
        truncate(io_buffer, 0) 
    end
end

# Debug printing
dprint(level, indent, msgs...)   = show(level, indent, false, msgs)
dprintln(level, indent, msgs...) = show(level, indent, true, msgs)
