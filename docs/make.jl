using StochasticGroundMotionSimulation
using Documenter

makedocs(;
    modules=[StochasticGroundMotionSimulation],
    authors="Peter Stafford <p.stafford@me.com>",
    repo="https://github.com/pstafford/StochasticGroundMotionSimulation.jl/blob/{commit}{path}#L{line}",
    sitename="StochasticGroundMotionSimulation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pstafford.github.io/StochasticGroundMotionSimulation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pstafford/StochasticGroundMotionSimulation.jl",
)
