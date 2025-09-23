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

using Documenter, CStarSurfaces

if format == "html"
    makedocs(
        sitename = "CStarSurfaces",
        pages = [
            "CStarSurfaces.jl" => "index.md",
        ],

        deploydocs(
            repo = "github.com/justus-springer/CStarSurfaces.jl.git",
        )
    )

elseif format == "thesis"
    makedocs(
        sitename = "CStarSurfaces",
        format = Documenter.LaTeX(platform = "none"),
        pages = [
            "CStarSurfaces.jl" => "index.md",
        ],
    )

    filename = joinpath(@__DIR__, "build", "CStarSurfaces.tex")

    txt = read(filename, String)


    # Add chapter heading
    txt = "\\chapter{CStarSurfaces.jl}\n\\label{apx:julia_cstar_surfaces}\n" * txt

    # Add Tex root directive
    txt = "%!TEX root = thesis.tex\n\n" * txt

    # Remove line breaks before examples
    txt = replace(txt, r"\n+(\textbf{Example})" => s"\n\1")

    # Remove math mode for references and add tilde
    txt = replace(txt, r" \\\( (\\ref\{\S+\}) \\\)" => s"~\1")

    write(filename, txt)

else
    println("Invalid format. Please choose 'html' or thesis")
    exit(1)
end
