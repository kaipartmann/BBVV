using BBVV
using Documenter

DocMeta.setdocmeta!(BBVV, :DocTestSetup, :(using BBVV); recursive=true)

makedocs(;
    modules=[BBVV],
    authors="Kai Partmann <kai.partmann@gmail.com> and contributors",
    repo="https://github.com/kaipartmann/BBVV.jl/blob/{commit}{path}#{line}",
    sitename="BBVV.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kaipartmann.github.io/BBVV.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kaipartmann/BBVV.jl",
    devbranch="main",
)
