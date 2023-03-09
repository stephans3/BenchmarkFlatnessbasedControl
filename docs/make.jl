using BenchmarkFlatnessbasedControl
using Documenter

DocMeta.setdocmeta!(BenchmarkFlatnessbasedControl, :DocTestSetup, :(using BenchmarkFlatnessbasedControl); recursive=true)

makedocs(;
    modules=[BenchmarkFlatnessbasedControl],
    authors="Stephan Scholz",
    repo="https://github.com/stephans3/BenchmarkFlatnessbasedControl.jl/blob/{commit}{path}#{line}",
    sitename="BenchmarkFlatnessbasedControl.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stephans3.github.io/BenchmarkFlatnessbasedControl.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stephans3/BenchmarkFlatnessbasedControl.jl",
    devbranch="main",
)
