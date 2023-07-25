# ---------------------------------------------------------------------------------------
# Parsing dates
function exceldatetodate(exceldate::Real)
    t, d = modf(exceldate)
    return Dates.DateTime(1899, 12, 30) + Dates.Day(d) + Dates.Millisecond((floor(t * 86400000)))
end

# ---------------------------------------------------------------------------------------
# Parsing XML summary from ICM+
function parseXML(Parameters, Selections, Results)

    selections = Dict()
    results = Dict()

    for i in 1:length(Selections)
        sel = elements(Selections[i])
        sel_name = Selections[i]["Name"]
        sel_st = sel[1]["StartTime"]
        sel_en = sel[1]["EndTime"]

        sel_st = DateTime(sel_st, dateformat"dd/mm/yyyy HH:MM:SS")
        sel_en = DateTime(sel_en, dateformat"dd/mm/yyyy HH:MM:SS")

        selections = merge!(selections, Dict(sel_name => [sel_st, sel_en]))
    end

    for i = eachindex(Results)
        res = elements(Results[i])
        res_name = Results[i]["Name"]
        name_surrogate = split(res_name, " ")

        # Some parameters name have units e.g. [mmHg] but Dict field naming does not support strings with non-alphanumeric chars
        if isletter(name_surrogate[2][1])
            res_name = name_surrogate[1] * "_" * name_surrogate[2]
        else
            res_name = name_surrogate[1]
        end

        res_val = Results[i]["Value"]
        results = merge!(results, Dict(res_name => parse(Float64, res_val)))
    end
    return results, selections
end

# ---------------------------------------------------------------------------------------
# Core function for obtaining data from ICM+ hdf5 files
function readCSF(filename)
    Data = Dict{String,Any}()
    recording_start_time = 0.0
    recording_end_time = 0.0

    fid = h5open(filename, "r")

    # Get attributes outside of datasets/groups
    completeAttrFlg = false
    try # For some reason some data exported to hdf5 does not have these attributes so need to calculate ad hoc
        recording_end_time = read_attribute(fid, "dataEndTime")
        recording_start_time = read_attribute(fid, "dataStartTime")
        dur = read_attribute(fid, "duration")
        dur = split(dur[1], " ")
        completeAttrFlg = true
    catch # Obtain these from ICP datetime & parse later on when it is loaded
        recording_start_time = [DateTime(2013)]
        recording_end_time = [DateTime(2013)]   # Pre-allocate so the type matches
    endŃ

    # Begin handling data in XML string - some files do not have infusion test output
    xml_obj = fid["aux/ICM+/icmtests"]
    xml_data = read(xml_obj)
    xml_string = String(xml_data[1])
    data = parsexml(xml_string)
    icm_data = root(data)
    vars = elements(icm_data)

    # Find which test is CSF infusion 
    ToolName = elements(parentelement(vars[1]))
    ToolName = elements(ToolName[1])
    for tn = 1:length(ToolName)
        tool = ToolName[tn]["Name"]
        tool == "CSF Infusion Test" ? (global ToolID = tn) : 0
    end

    SingleAnalysis = elements(vars[2])
    SingleAnalysis = SingleAnalysis[ToolID]
    SingleAnalysis = elements(SingleAnalysis)
    # Variables = elements(SingleAnalysis[1]) # Not currently used
    Selections = elements(SingleAnalysis[2])
    Parameters = elements(SingleAnalysis[3])
    Results = elements(SingleAnalysis[4])

    results, selections = parseXML(Parameters, Selections, Results)

    # This way of dereferencing is pain, there has to be a more elegant solution
    # Probably using eval()
    t_series_ds = fid["summaries/minutes/"]
    t_series = t_series_ds[:]
    numsamples = length(t_series)

    # Pre-allocate ?needed
    ICP = zeros(numsamples)
    AMP = zeros(numsamples)
    timestamp = zeros(numsamples)
    P0 = zeros(numsamples)
    AMP_P = zeros(numsamples)

    # Dereferencing named tuple...
    for i in 1:numsamples
        ICP[i] = t_series[i].ICP
        AMP[i] = t_series[i].AMP

        timestamp[i] = t_series[i].datetime
        # Some files do not have P_0 saved?
        try
            P0[i] = t_series[i].P0
            AMP_P[i] = t_series[i].AMP_P
        catch
        end
    end

    if !completeAttrFlg # If for some reason the recording start and end times are not saved, obtain them from timestamp data
        recording_start_time[1] = exceldatetodate(timestamp[1])
        recording_end_time[1] = exceldatetodate(timestamp[end])
        start_time = recording_start_time[1]
        end_time = recording_end_time[1]
    else
        start_time = DateTime(recording_start_time[1], dateformat"yyyy/mm/dd HH:MM:SS")
        end_time = DateTime(recording_end_time[1], dateformat"yyyy/mm/dd HH:MM:SS")
    end

    # Inconsistent naming
    try
        Data["E"] = results["Elasticity"]
    catch
        Data["E"] = results["Elastance"]
    end

    Data["I_b"] = results["CSF_production"]
    Data["Rcsf"] = results["Rcsf"]

    # Inconsistent naming
    try
        Data["P_0"] = results["Pss"]
    catch
        Data["P_0"] = results["pss"]
    end

    Data["ICP"] = ICP
    Data["AMP"] = AMP
    Data["P_p"] = results["ICP_plateau"]
    Data["P_b"] = results["ICP_baseline"]
    Data["T"] = [0:numsamples-1...] * 1 / 6

    Data["baseline_start_frame"] = round(Int, (selections["Baseline"][1] - start_time).value / 10000)
    Data["baseline_end_frame"] = round(Int, (selections["Baseline"][2] - start_time).value / 10000)
    Data["infusion_start_frame"] = round(Int, (selections["Infusion"][1] - start_time).value / 10000)
    Data["infusion_end_frame"] = round(Int, (selections["Infusion"][2] - start_time).value / 10000)
    Data["transition_start_frame"] = round(Int, (selections["Transition"][1] - start_time).value / 10000)
    Data["transition_end_frame"] = round(Int, (selections["Transition"][2] - start_time).value / 10000)
    Data["plateau_start"] = round(Int, (selections["Plateau"][1] - start_time).value / 10000)
    Data["plateau_end"] = round(Int, (selections["Plateau"][2] - start_time).value / 10000)
    Data["rec_dur_s"] = (end_time - start_time).value / 1000
    Data["start_time"] = start_time
    Data["end_time"] = end_time
    Data["I_inf"] = parse(Float64, Parameters[1]["Value"])
    Parameters[3]["Value"] == "yes" ? Data["one_needle"] = 1 : Data["one_needle"] = 0
    Data["one_needle"] == 1 ? Data["Rn"] = results["Needle_resist."] : Data["Rn"] = 0

    return Data
end

# ---------------------------------------------------------------------------------------
# Plot infusion study trace
function plotmodel(I_b, E, P_0, Pss, μ, σ, cscheme, fmodel)

    icp = Data["ICP"]
    infstart = Data["infusion_start_frame"]
    infend = Data["infusion_end_frame"]
    P_b = Data["P_b"]
    plateau_end = Data["plateau_end"]
    plateau_end > infend ? endidx = plateau_end : endidx = infend

    Rout = (P_b - P_0) / I_b

    p = plot()
    color = p[:background_color]

    if sum([color.r, color.g, color.b]) > 1.5
        grayval = 0.8
        dashcolor = :black
    else
        grayval = 2.0
        dashcolor = :white
    end
    # sum([color.r, color.g, color.b]) > 1.5 ? (grayval = 0.7, dashcolor=:black) : (grayval = 2.0, dashcolor=:white)
    newvals = [color.r, color.g, color.b, color.alpha] .*= [grayval, grayval, grayval, 1.0]
    newcol = RGBA{Float64}(newvals[1], newvals[2], newvals[3], newvals[4])

    vspan!([infstart, infend], fillcolor=newcol, alpha=0.2, linewidth=1, linestyle=:dash, linecolor=dashcolor, label="Infusion")

    xtks = LinRange(0, infend, 10)
    xtklabels = round.(collect(xtks) ./ 6, digits=1)

    mean([color.r, color.g, color.b]) > 0.5 ? icp_col = RGBA{Float64}(0.0, 0.0, 0.0) : icp_col = RGBA{Float64}(1.0, 1.0, 1.0)

    plot!(icp, lw=2, grid=false, xticks=(xtks, xtklabels), linecolor=icp_col, legend=:outertopright, label="ICP", ylims=[minimum(icp) * 0.9, maximum(icp) * 1.1], xlims=(firstindex(icp), endidx))

    Pm = zeros(endidx)
    Pmodel, rmserr = calc_model_plot(Data["I_b"], Data["E"], Data["P_0"], Data["P_0"])
    Pm[firstindex(Pm):infend] .= Pmodel
    Pm[infend+1:end] .= Pm[infend]
    Pm[firstindex(Pm):infstart] .= P_b
    plot!(Pm, label="Gradient descent", lw=3)

    Pm = zeros(endidx)
    Pmodel, rmserr = calc_model_plot(I_b, E, P_0, P_0)
    Pm[firstindex(Pm):infend] .= Pmodel
    Pm[infend+1:end] .= Pm[infend]
    Pm[firstindex(Pm):infstart] .= P_b
    num_iter = 1
    w, ci = getCI(μ, σ, num_iter)
    plot!(Pm, ribbon=w, fillalpha=0.1, label="Bayesian", lw=3)

    title!("I_b = $(round(I_b,digits=2))\n" * "Rcsf = $(round((Rout),digits=2))\n" * "E = $(round((E),digits=2))\n" * "P_0 = $(round((P_0),digits=2))\n" * "error = $rmserr", titlealign=:left, titlefontsize=8, xlabel="Time [min]", ylabel="ICP [mmHg]")
end

# ---------------------------------------------------------------------------------------
# Ad hoc calc of the model curve
function calc_model_plot(I_b, E, P_0, Pss)
    infstart = Data["infusion_start_frame"]
    infend = Data["infusion_end_frame"]
    I_inf = Data["I_inf"]
    Rn = Data["Rn"]
    ΔP = Data["P_b"] - P_0
    icp = Data["ICP"]
    It = I_b + I_inf
    Pm = zeros(infend) .+ Data["P_b"]
    errorVal = 0.0

    for i = infstart:infend
        t = (i - infstart) / 6
        y = It * ΔP / (I_b + (I_inf * exp(-E * It * t))) + P_0 + (I_inf * Rn)
        Pm[i] = y
        # errorVal += (icp[i] - y)^2
    end

    icp_trace = Data["ICP"][infstart:infend]
    Pm_inf = Pm[infstart:infend]

    # rmserr = 100 * sqrt(errorVal) / length(Pm) / abs(mean(icp[infstart:infend]))
    rmserr = rmsd(Pm_inf, icp_trace, normalize=true)
    return Pm, rmserr
end

# ---------------------------------------------------------------------------------------
# Obtain confidence intervals for ribbon plot of uncertainty using boostrapping
function getCI(μ, σ, num_iter)
    infstart = Data["infusion_start_frame"]
    infend = Data["infusion_end_frame"]
    icp = Data["ICP"][infstart:infend]
    Pmodel = zeros(infend)
    model_err = zeros(num_iter)
    numvars = length(μ)
    θ = zeros(numvars, num_iter)

    for i = 1:numvars
        d = Normal(μ[i], σ[i])
        θ[i, :] = rand(d, num_iter)
    end

    if numvars == 3
        θ[1, :] = (Data["P_b"] .- θ[3, :]) ./ θ[1, :]
        for j = 1:num_iter
            θᵢ = θ[:, j]
            Pmodel = calc_model_plot(θᵢ[1], θᵢ[2], θᵢ[3], θᵢ[3])[1]
            Pmodel = Pmodel[infstart:end]
            model_err[j] = mean(abs.(Pmodel .- icp))
        end
    elseif numvars == 2
        θ[1, :] = (Data["P_b"] .- Data["P0_static"]) ./ θ[1, :]
        for j = 1:num_iter
            θᵢ = θ[:, j]
            Pmodel = calc_model_plot(θᵢ[1], θᵢ[2], Data["P0_static"], Data["P0_static"])[1]
            Pmodel = Pmodel[infstart:end]
            model_err[j] = mean(abs.(Pmodel .- icp))
        end
    else
        θ[1, :] = (Data["P_b"] .- θ[4, :]) ./ θ[1, :]
        for j = 1:num_iter
            θᵢ = θ[:, j]
            Pmodel = calc_model_plot(θᵢ[1], θᵢ[2], θᵢ[3], θᵢ[4])[1]
            Pmodel = Pmodel[infstart:end]
            model_err[j] = mean(abs.(Pmodel .- icp))
        end
    end

    x̂ = mean(model_err)
    s = std(model_err)
    z = 0.95
    n = num_iter

    ci_low = x̂ + z * s / sqrt(n)
    ci_high = x̂ + z * s / sqrt(n)
    y1 = Pmodel .- ci_low
    y2 = Pmodel .+ ci_high
    w = (y2 .- y1) ./ 2

    return w, ci_low
end

# ---------------------------------------------------------------------------------------
# CSF dynamics model
function dyn_model(Rcsf, E, P_0)
    I_b = (Data["P_b"] - P_0) / Rcsf # CSF formation rate

    infstart = Data["infusion_start_frame"]
    infend = Data["infusion_end_frame"]
    I_inf = Data["I_inf"]

    Rn = Data["Rn"] # Needle resistance (one-needle)
    ΔP = Data["P_b"] - P_0
    It = I_b + I_inf
    # Pm = zeros(infend - infstart + 1)
    # Pm = zeros(eltype(Data["ICP"]), infend - infstart + 1)
    Pm = Vector{Any}(undef, infend - infstart + 1)

    for i = infstart:infend
        t = (i - infstart) / 6
        Pm[i-infstart+1] = It * ΔP / (I_b + (I_inf * exp(-E * It * t))) + P_0 + (I_inf * Rn)
    end

    return Pm
end

# ---------------------------------------------------------------------------------------
# Define the Turing.jl model
@model function curve_fitting(data, ::Type{T}=Vector{Float64}) where {T}

    # Priors
    E ~ truncated(Normal(params_means[2], params_stddevs[2]), lower_bounds[2], upper_bounds[2])
    P0 ~ truncated(Normal(params_means[3], params_stddevs[3]), lower_bounds[3], upper_bounds[3])
    Rout ~ truncated(Normal(params_means[1], params_stddevs[1]), lower_bounds[1], upper_bounds[1])

    # Accounting for production rate range:
    # 1) Soft Constraint (Potential) - penalty function
    # 2) custom distribution ranges for P0 and Rout based on Ib bounds
    # 3) optimise Ib and derive Rout/P0 from this (but Rout is robust, and Ib is sensitive)
    # 4) post-processing of the chain/posterior e.g. rejection sampling

    # Soft Constraint (Potential)
    Ib = (Data["P_b"] - P0) / Rout
    min_value = 0.01
    max_value = 1.00
    acceptable_Ib_range = (min_value, max_value) # Define range
    is_acceptable = (acceptable_Ib_range[1] <= Ib <= acceptable_Ib_range[2])
    constraint_penalty = is_acceptable ? 0.0 : Inf # Infinite penalty if outside range
    # Turing.@addlogprob! -constraint_penalty # Subtract penalty from log probability

    # ReLU-like penalty
    # relu_penalty = max(0, abs(Ib - mean(acceptable_Ib_range)) - 0.5 * (max_value - min_value))
    # Turing.@addlogprob! -relu_penalty # Subtract penalty from log probability (change to -relu_penalty for ReLU)
    
    # Sigmoid-like penalty
    # sigmoid_penalty = 1.0 / (1.0 + exp(-10 * (abs(Ib - mean(acceptable_Ib_range)) - 0.5 * (max_value - min_value))))
    # Turing.@addlogprob! -sigmoid_penalty

    # E ~ truncated(Laplace(params_means[2], params_stddevs[2]/sqrt(2)), lower_bounds[2], upper_bounds[2])
    # P0 ~ truncated(Laplace(params_means[3], params_stddevs[3]/sqrt(2)), lower_bounds[3], upper_bounds[3])
    # Rout ~ truncated(Laplace(params_means[1], params_stddevs[1]/sqrt(2)), lower_bounds[1], upper_bounds[1])

    # Standard deviation for the likelihood, i.e. how much to trust the data
    # Higher values -> less trust
    sigma ~ truncated(Normal(0, 1), 0, Inf)

    # Compute model predictions
    predicted = dyn_model(Rout, E, P0)

    # Likelihood (Gaussian)
    for i in 1:length(data)
        data[i] ~ Normal(predicted[i], sigma)
    end
end