
push!(LOAD_PATH,"./src");
using Documenter, EntanglementDist;

makedocs(format = :html, sitename = "Entanglement Distillation",
    pages = [
        "index.md",
        "intro.md",
        "optimize.md",
        "quantum.md",
        "schemes.md",
        "states.md",
        "examplePPT.md",
        "exampleKEXT.md",
        "exampleSeesaw.md"
        ]
)


