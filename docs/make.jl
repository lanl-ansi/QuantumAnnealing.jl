using Documenter, QuantumAnnealing

makedocs(
    warnonly = Documenter.except(:linkcheck),
    modules = [QuantumAnnealing],
    sitename = "QuantumAnnealing",
    authors = "Zach Morrell, Carleton Coffrin, Marc Vuffray",
    pages = [
        "Home" => "index.md",
        "Library" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/QuantumAnnealing.jl.git",
)