using Documenter, DocumenterCitations, CStarSurfaces

bib = CitationBibliography(
    joinpath(@__DIR__, "references.bib");
    style = :numeric
)

makedocs(
    plugins = [bib],
    format = Documenter.HTML(prettyurls = false),
    sitename = "CStarSurfaces.jl",
    pages = [
        "Home" => "index.md",
        "General" => [
            "constructors.md"
            "attributes.md"
        ]
    ]
)
