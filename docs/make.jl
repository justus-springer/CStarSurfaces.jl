using Documenter, DocumenterCitations, CStarSurfaces, Oscar

DocMeta.setdocmeta!(CStarSurfaces, :DocTestSetup, :(using CStarSurfaces, Oscar); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "references.bib");
    style = :numeric
)

makedocs(
    plugins = [bib],
    #modules = [CStarSurfaces],
    doctest = true,
    format = Documenter.HTML(prettyurls = false),
    sitename = "CStarSurfaces.jl",
    pages = [
        "Home" => "index.md",
        "Surfaces with torus action" => "surfaces_with_torus_action.md",
        "Normal forms" => "admissible_operations.md",
        "Mori Dream Spaces" => "mori_dream_spaces.md",
        "Database functionality" => "database_functionality.md",
        "Index" => "docs_index.md"
    ]
)

deploydocs(
    repo = "github.com/justus-springer/CStarSurfaces.jl.git",
)
