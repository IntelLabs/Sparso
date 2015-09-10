# This file describes all the possible exceptions that may be raised during 
# the execution of Sparse Accelerator, as well as their handlers.
# Usually, an exception handler prints out a nice diagnostic message, and
# returns to the original user code without optimizing it.

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

type UnknownExprDistributivity <: Exception
    head
    args
end

type UnknownASTDistributivity <: Exception
    ast
end

type UnknownTypeToReorder <: Exception
    sym
    typ
end

type AsymmetricMatrixMarketFile <: Exception
    filename
end

type PostPatternReplacementFailure <: Exception
    pattern
end

type MatrixPropertiesUnavail <: Exception
    symexpr
end