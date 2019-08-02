using Test
using EconLH.LatexLH

testDir = joinpath(@__DIR__, "test_files");

@testset "LatexLH" begin
	@testset "CellColor" begin
		x = LatexLH.CellColor("blue", 40);
		s = color_string(x);
		# println(s)
		@test s == "\\cellcolor{blue!40}"

		x = LatexLH.CellColor("green", 0)
		s = color_string(x);
		# println(s)
		@test s == ""
	end


	@testset "Cell" begin
		c = LatexLH.Cell("abc", 1, 'c', LatexLH.CellColor("blue", 40))
		@test c.text == "abc"
		@test cell_string(c) == "\\cellcolor{blue!40}abc"

		c = LatexLH.Cell("abc");
		@test c.text == "abc"

		c = LatexLH.Cell("abc", 'c');
		@test c.text == "abc"
		@test c.align == 'c'
		@test cell_string(c) == "abc"

		c = LatexLH.Cell("abc", 2, 'c', LatexLH.CellColor("blue", 40))
		@test cell_string(c) == "\\multicolumn{2}{c}{\\cellcolor{blue!40}abc}"
	end


	@testset "Table" begin
		nCols = 5;
		tb = LatexLH.Table(nCols, "lSSSS");
		LatexLH.add_row!(tb, " & \\multicolumn{3}{c}{Heading 1} & ");
		@test length(tb.bodyV) == 1

		cellV = [LatexLH.Cell("cell1"),  LatexLH.Cell("cell2", 4, 'c')]
		rowStr = LatexLH.make_row(tb, cellV);
		@test rowStr == "cell1 & \\multicolumn{4}{c}{cell2}"
		LatexLH.add_row!(tb, rowStr);
		@test length(tb.bodyV) == 2

		LatexLH.add_row!(tb, "\\cline{2-5}")
		LatexLH.add_row!(tb, "x1 & 1.2 & 2.3 & 3.4 & 4.5")

		filePath = joinpath(testDir, "TableTest.tex")
		if isfile(filePath)
			rm(filePath)
		end
		write_table(tb, filePath)
		@test isfile(filePath)
	end


	@testset "ParameterTable" begin
		pt = ParameterTable();
		add_row!(pt, "p1", "param 1", "1.23");
		# Not robust. Skipping 2 header rows
		@test nrows(pt.tb) == 3
		@test pt.tb.bodyV[3] == "p1 & param 1 & 1.23"
		add_row!(pt, "p2", "param 2", "2.34");
		write_table(pt, joinpath(testDir, "parameter_table_test.tex"))
	end
end



# ----------
