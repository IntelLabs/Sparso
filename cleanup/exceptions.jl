# This file describes all the possible exceptions that may be raised during 
# the execution of Sparse Accelerator, as well as some common exception handling.
# Usually, the exception handling is to print out a nice diagnostic message, and
# return to the original user code without optimizing it.

type UndescribedFunction <: Exception
    module_name     :: String
    function_name   :: String
    parameter_types :: Tuple
end

type UnresolvedFunction <: Exception
    head
    args
end

type UnhandledModuleOrFunctionName <: Exception
    arg
end

type UnknownCallFormat <: Exception
    head
    args
    parameter_types :: Tuple
end
