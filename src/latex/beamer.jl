# Writing beamer slides

"""
	$(SIGNATURES)

Make a file that contains a set of slides. Each slide is just a title with one figure. Input is a vector of tuples. Each contains a title and a file path for the figure file.
Optionally, preprend a common base directory to each path.
"""
function write_figure_slides(io :: IO,  slideV;
    baseDir :: AbstractString = "")
    for slide in slideV
        figPath = joinpath(baseDir, slide[2]);
        write_figure_slide(io, slide[1], figPath);
    end
end


"""
	$(SIGNATURES)

Write lines for a slide that contains one figure to `io`.
"""
function write_figure_slide(io :: IO, title, figPath)
    for line ∈ figure_slide(title, figPath)
        println(io, line);
    end
end

"""
	$(SIGNATURES)

Latex lines for slide with a single figure.
"""
function figure_slide(title, fPath)
    v = Vector{String}();
    push!(v, begin_frame(title));
    for line ∈ include_figure(fPath)
        push!(v, line);
    end
    push!(v, end_frame());
    return v
end

function begin_frame(title)
    return "\\begin{frame}{$title}";
end

function end_frame()
    return "\\end{frame}"
end

function include_figure(fPath)
    return [
        "\\begin{center}",
        "\\includegraphics[width=0.95\\textwidth]{$fPath}",
        "\\end{center}"
    ]    
end

# ------------