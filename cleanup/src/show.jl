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

function show_structure(
    node :: Any,
    space :: String
)
    typ = typeof(node)
    if typ == Expr
        println(io_buffer, space, "Expr ", node.head, " [", node.typ, "]")
        args                = node.args
        args_are_statements = false
        if node.head == :lambda
            # args[3] is the body Expr.
            # TODO: LambdaHandling package should identify the body
            args                = node.args[3].args
            args_are_statements = true
        end
        for i in 1 : length(args)
            arg = args[i]
            type_of_arg = typeof(arg)
            if (type_of_arg == Expr && arg.head == :line) 
                println(io_buffer, space, "    ", arg)
            else 
                arg_is_insignificant = (type_of_arg == LabelNode) || 
                    (type_of_arg == LineNumberNode) || (type_of_arg == NewvarNode)
                if args_are_statements && !arg_is_insignificant
                    println(io_buffer, space, "    ", "// ", arg)
                end
                print(io_buffer, space, "    ", i, ": ")
                show_structure(arg, "")
                if args_are_statements
                    println(io_buffer, "")
                end
            end
        end
    else
        println(io_buffer, space, node, " [", typeof(node), "]")
    end
end

function show(structure :: Bool, msg)
    if structure
        show_structure(msg, "")
    else
        print(io_buffer, msg)
    end
end

function show(level :: Int, indent :: Int, new_line :: Bool, structure :: Bool, msgs...)
    if(level >= verbosity)
        message_tuple = msgs[1]
        for msg in message_tuple
            show(structure, msg)
        end
        
        seekstart(io_buffer)
        while !eof(io_buffer)
            [print(" ") for i = 1 : indent * 4]
            s = readline(io_buffer)
            print_unescaped(STDOUT, s)
        end
        if new_line 
            println(STDOUT, "") 
        end

        # Resize the buffer to be empty
        truncate(io_buffer, 0) 
    end
end

# Debug printing
dprint(level, indent, msgs...)   = show(level, indent, false, false, msgs)
dprintln(level, indent, msgs...) = show(level, indent, true, false, msgs)
dsprint(level, indent, msgs...)   = show(level, indent, false, true, msgs)
dsprintln(level, indent, msgs...) = show(level, indent, true, true, msgs)

