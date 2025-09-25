if length(ARGS) > 1
    println("Usage: julia make.jl <format>")
    println("Format can be 'html', 'thesis'")
    exit(1)
end

format = isempty(ARGS) ? "html" : ARGS[1]

if format ∉ ["html", "thesis"]
    println("Usage: julia make.jl <format>")
    println("Format can be 'html', 'thesis'")
    exit(1)
end

@info "make.jl: Building documentation for format \"$format\""

using Documenter, CStarSurfaces, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

if format == "html"
    makedocs(
        sitename = "CStarSurfaces",
        pages = [
            "CStarSurfaces.jl" => "index.md",
            "``\\mathbb{C}^*``-surfaces" => "cstar_surfaces.md",
            "Local properties" => "local_properties.md",
            "Classifications" => "classifications.md",
        ],
    )

    deploydocs(
        repo = "github.com/justus-springer/CStarSurfaces.jl.git",
    )

elseif format == "thesis"
    makedocs(
        sitename = "CStarSurfaces",
        format = Documenter.LaTeX(platform = "none"),
        pages = [
            "CStarSurfaces.jl" => "index.md",
            "C*-surfaces" => "cstar_surfaces.md",
            "Local properties" => "local_properties.md",
            "Classifications" => "classifications.md",
        ],
        plugins = [bib],
    )

    filename = joinpath(@__DIR__, "build", "CStarSurfaces.tex")

    txt = read(filename, String)

    # Remove bibliography
    start_idx = findfirst("{\\raggedright% @bibliography", txt)
    end_idx = findfirst("}% end @bibliography", txt)
    if start_idx !== nothing && end_idx !== nothing
        txt = txt[1:first(start_idx)-1] * txt[last(end_idx)+1:end]
    end

    # Add chapter heading
    txt = "\\chapter{CStarSurfaces.jl}\n\\label{apx:julia_cstar_surfaces}\n" * txt

    # Add Tex root directive
    txt = "%!TEX root = thesis.tex\n\n" * txt

    # Remvoe line breaks before equations
    txt = replace(txt, r"\n+(\\begin\{equation\*\})" => s"\n\1")

    # Remove math mode for references and add tilde
    txt = replace(txt, r" \\\( (\\ref\{\S+\}) \\\)" => s"~\1")

    # Fix citations
    txt = replace(txt, r"\[\\hyperref\[doc:(\w+)\]\{\d+\}\]" => s"\\cite{\1}")

    write(filename, txt)

else
    println("Invalid format. Please choose 'html' or thesis")
    exit(1)
end
