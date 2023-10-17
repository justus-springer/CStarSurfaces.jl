using Documenter, DocumenterCitations, CStarSurfaces

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
    ]
)

deploydocs(
    repo = "github.com/justus-springer/CStarSurfaces.jl.git",
)
