#! /bin/env julia


#####################################################################
using ArgParse


#####################################################################
mutable struct Model
	#Model                  LogL         AIC      w-AIC        AICc     w-AICc         BIC      w-BIC
	model::String
	scores::Dict
end


#####################################################################
function read_arguments()
	s = ArgParseSettings()
	@add_arg_table s begin
		"-i", "--in", "--infile"
			help = "infile"
			required = true
		"--metric"
			help = "AIC, BIC, AICc"
			required = true
	end
	parsed_args = parse_args(s)
	infile = parsed_args["in"]
	metric = parsed_args["metric"]
	return([infile, metric])
end


#####################################################################
infile, metric = read_arguments()


#####################################################################
model_objs = Vector()

open(infile) do in_fh
	is_log = false
	is_iqtree = false
	global is_start = false
	headers = nothing
	for line in eachline(in_fh)
		if occursin(r"^Model[ ]+LogL", line)
			is_start = true
			is_iqtree = true
		elseif occursin(r"^ No. Model", line)
			is_start = true
			is_log = true
		elseif is_start && occursin(r"^$", line)
			is_start = false
		end

		if is_start
			line = replace(line, r"^\s+" => "")
			line_arr = split(line, r"\s+")
			ori_line_arr = copy(line_arr)
			# discard the 1st element in line_arr
			is_log && deleteat!(line_arr, 1)
			if match(r"^Model", line_arr[1]) !== nothing
				headers = line_arr[2:length(line_arr)]
			else
				if is_log && (! occursin(r"^\d+", ori_line_arr[1]))
					continue
				end
				model = line_arr[1]
				scores = line_arr[2:length(line_arr)]
				scores = filter(x->!occursin(r"^[+-]$",x), scores)
				scores = map(x->parse(Float64,x), scores)
				scores = Dict(zip(headers, scores))
				model_obj = Model(model, scores)
				push!(model_objs, model_obj)
			end
		end
	end
	#println(model_objs)
end


#####################################################################
sorted_model_objs = sort(by=x->x.scores[metric], model_objs)
for model_obj in sorted_model_objs
	println(join([model_obj.model, model_obj.scores[metric]], "\t"))
end


