module StructureAnalysisDummy

    name = "Dummy"

    function print_action(s)
        dprintln(1, 0, s.args) 
    end

    transfer_rules = [
        ((:call, GlobalRef(Main, :*), Any, Any), print_action),
        ((:call, GlobalRef(Main, :+), Any, Any), print_action),
    ]
end
