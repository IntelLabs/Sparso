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
    node        :: Any,
    space       :: String,
    symbol_info :: Sym2TypeMap, 
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
                show_structure(arg, string(space, "    "), symbol_info)
                if args_are_statements
                    println(io_buffer, "")
                end
            end
        end
    else
        print(io_buffer, space, node, " [", typ, "]")
        typ1 = type_of_ast_node(node, symbol_info)
        if typ != typ1
            println(io_buffer, " [", typ1, "]")
        else
            println(io_buffer, "")
        end
    end
end

function show(
    structure   :: Bool, 
    symbol_info :: Sym2TypeMap, 
    msg
)
    if structure
        show_structure(msg, "", symbol_info)
    else
        if typeof(msg) <: Exception
            print(io_buffer, msg)
            show_backtrace(io_buffer, catch_backtrace())
            println(io_buffer, "")
        elseif typeof(msg) <: Pattern
            flds = fieldnames(msg)
            for i = 1 : length(flds)
                fld = getfield(msg, flds[i])
                if i < length(flds)
                    println(io_buffer, fld)
                else
                    print(io_buffer, fld)
                end
            end
        elseif typeof(msg) <: Vector{Action}
            for action in msg
                println(io_buffer, "Insert action:")
                for stmt in action.new_stmts
                    println(io_buffer, "    ", stmt.expr)
                end
                if typeof(action) <: InsertBeforeStatement
                    println(io_buffer, "before BB ", action.bb.label, " statement")
                    println(io_buffer, "    ", action.bb.statements[action.stmt_idx].expr)
                elseif typeof(action) <: InsertBeforeLoopHead
                    print(io_buffer, action.outside_loop ? "outisde" : "inside", " loop of BBs [")
                    total = 0
                    for i in action.loop.members
                        total = total + 1
                        print(io_buffer, i, 
                            total == length(action.loop.members) ? "" : ", ")
                    end
                    println(io_buffer, "] before loop head BB ", action.loop.head)
                else
                    assert(typeof(action) == InsertOnEdge)
                    println(io_buffer, "on the edge BB ", action.from_bb.label, " --> BB ", action.to_bb.label)
                end
                println(io_buffer, "")
            end
        else
            print(io_buffer, msg)
        end
    end
end

function show(
    level       :: Int, 
    indent      :: Int, 
    new_line    :: Bool, 
    structure   :: Bool,
    symbol_info :: Sym2TypeMap, 
    msgs...
)
    if(level >= verbosity)
        message_tuple = msgs[1]
        for msg in message_tuple
            show(structure, symbol_info, msg)
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
dprint(level, indent, msgs...)   = show(level, indent, false, false, Sym2TypeMap(), msgs)
dprintln(level, indent, msgs...) = show(level, indent, true, false, Sym2TypeMap(), msgs)
dsprint(level, indent, symbol_info, msgs...)   = show(level, indent, false, true, symbol_info, msgs)
dsprintln(level, indent, symbol_info, msgs...) = show(level, indent, true, true, symbol_info, msgs)