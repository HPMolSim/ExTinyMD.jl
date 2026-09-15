using Documenter, ExTinyMD

DocMeta.setdocmeta!(ExTinyMD, :DocTestSetup, :(using ExTinyMD, StaticArrays);
                    recursive = true)

makedocs(
    sitename = "ExTinyMD.jl",
    modules  = [ExTinyMD],
    format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home"           => "index.md",
        "MD Core"        => "md_core.md",
        "Interactions"   => "interactions.md",
        "Electrostatics" => "electrostatics.md",
        "API Reference"  => "api.md",
    ],
    checkdocs = :exports,
    warnonly  = false,
)

deploydocs(repo = "github.com/HPMolSim/ExTinyMD.jl", devbranch = "main")
