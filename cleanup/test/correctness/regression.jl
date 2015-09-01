# This file is used for preventing regressions and keeping the quality of the 
# repository higher every checkin. 
# For the module you write, please write as many unit test cases as possible, and
# add the corresponding entries for them -- so that your code will not be accidently
# affected by other people's checkins!

@doc """"
A pattern to match with a test's output. The pattern is a regular expression.
The comment explains the pattern, e.g. what the pattern is supposed to do, as 
a friendly message when the pattern does not match.
"""
immutable TestPattern
    pattern :: Regex
    comment :: String
end

immutable Test
    name     :: String
    command  :: Cmd
    patterns :: Vector{TestPattern}
end

const sanity_test1 = Test(
    "sanity-test1",
    `julia sanity-test1.jl`,
    [
        TestPattern(r" {4,4}1 tab",
                     "Test tabbed output facility: 1 tab."
        ),
        TestPattern(r" {8,8}2 tabs",
                     "Test tabbed output facility: 2 tabs."
        ),
        TestPattern(r" {12,12}3 tabs",
                     "Test tabbed output facility: 3 tabs."
        )
    ]
)

const sanity_test2 = Test(
    "sanity-test2",
    `julia sanity-test2.jl`,
    [
        TestPattern(r"I do nothing!",
                     "Test if optimization framework adds a new optimization pass."
        )
    ]
)

const sanity_test3 = Test(
    "sanity-test3",
    `julia sanity-test3.jl`,
    [
        TestPattern(r"New AST:(.|\n)*return A::Base.SparseMatrix.SparseMatrixCSC.*x::Array",
                     "Test if optimization framework invokes SparseAccelerator to generate a new AST."
        )
    ]
)

const context_test2 = Test(
    "context-test2",
    `julia context-test2.jl small-diag.mtx`,
    [
        TestPattern(r"Original:(.|\n)*sum of x=-1.577312043410735e-5(.|\n)*rel_err=6.381531131954942e-13(.|\n)*With manual context-sensitive optimization:(.|\n)*sum of x=-1.577312043410735e-5(.|\n)*rel_err=6.381531131217335e-13",
                     "Test pcg_symgs_with_context_opt"
        )
    ]
)

const liveness_test1 = Test(
    "liveness-test1",
    `julia liveness-test1.jl small-diag.mtx`,
    [
        TestPattern(r"Def",
                     "Test liveness for cg."
        )
    ]
)

const tests = [
    sanity_test1,
    sanity_test2,
    sanity_test3,
    context_test2,
    liveness_test1
]

fail = 0
succ = 0
old_stderr = STDERR
old_stdout = STDOUT
for test in tests
    print("Testing ", test.name)
    log  = string(test.name, ".log")
    file = open(log, "w+")
    redirect_stderr(file)
    redirect_stdout(file)
    output = ""
    try
        output = readall(test.command)
    catch ex
    end
    successful = true
    for pattern in test.patterns
        m = match(pattern.pattern, output)
        if m == nothing
            comment = pattern.comment
            write(file, "\n****** Result:\n", output)
            write(file, "\n****** Failed in pattern:\n\tPattern: ",
                  string(pattern.pattern), "\n\tComment: ", comment)
            close(file)
            redirect_stderr(old_stderr)
            redirect_stdout(old_stdout)
            println(": FAIL. See ", log)
            successful = false
            break
        end
    end
    if successful
        succ = succ + 1
        close(file)
        redirect_stderr(old_stderr)
        redirect_stdout(old_stdout)
        rm(log)
        println(": Pass")
    else
        fail = fail + 1
    end
end

println("Total: ", fail + succ)
println("Pass : ", succ)
println("Fail : ", fail)
flush(STDOUT)