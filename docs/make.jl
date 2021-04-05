using StochasticGroundMotionSimulation
using Documenter

makedocs(;
    modules=[StochasticGroundMotionSimulation],
    authors="Peter Stafford <p.stafford@me.com>",
    # repo="https://github.com/pstafford/StochasticGroundMotionSimulation.jl/blob/{commit}{path}#L{line}",
    sitename="StochasticGroundMotionSimulation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pstafford.github.io/StochasticGroundMotionSimulation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Model Components" => [
            "Fourier Spectral Parameters" => "fourier_parameters.md",
            "Oscillator Parameters" => "sdof_parameters.md",
            "Random Vibration Parameters" => "random_vibration_parameters.md",
        ],
        "Module Functionality" => [
            "Fourier Components" => [
                "Fourier Amplitude Spectrum" => "fourier_spectrum.md",
                "Source spectrum" => "source_spectrum.md",
                "Path scaling" => "path_scaling.md",
                "Site scaling" => "site_scaling.md",
            ]
        ]

    ],
)

deploydocs(;
    repo="github.com/pstafford/StochasticGroundMotionSimulation.jl.git",
)
