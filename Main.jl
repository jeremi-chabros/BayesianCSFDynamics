using Pkg
Pkg.activate(".")
using Turing, CSV, EzXML, Glob, Dates, StatsPlots, HDF5, Measures, DataFrames, StatsBase, Logging
Logging.global_logger(Logging.SimpleLogger(stdout, Logging.Error))

include("Utils.jl")
# theme(:dracula; palette=palette(:seaborn_colorblind))
theme(:dracula)
# theme(:wong2)

results = DataFrame(
    recording_name=Vector{String31}(),
    Rout_GD=Float64[],
    E_GD=Float64[],
    P0_GD=Float64[],
    Ib_GD=Float64[],
    Rout_mean=Float64[],
    E_mean=Float64[],
    P0_mean=Float64[],
    Ib_mean=Float64[],
    Rout_std=Float64[],
    E_std=Float64[],
    P0_std=Float64[],
    Ib_std=Float64[],
    Rout_MAP=Float64[],
    E_MAP=Float64[],
    P0_MAP=Float64[],
    Ib_MAP=Float64[],
    Rhat=Float64[],
    NRMSE_Bayes_mean=Float64[],
    NRMSE_Bayes_MAP=Float64[],
    NRMSE_GD=Float64[],
    acceptance_rate=Float64[],
    execution_time=Float64[])

# Load data
datapath = "/Users/jjc/CSF/Recordings/"
path = pwd();
savepath = "/Users/jjc/CSF/"
files = glob("*.hdf5", datapath)
filenames_trunc = Vector{String}()
for i = eachindex(files)
    push!(filenames_trunc, files[i][length(datapath)+1:end-5])
end

ppath = "/Users/jjc/BayesCSF/ModelPlots/"
processed_files = glob("*.png", ppath)

# contains(processed_files, "")
processed_trunc = Vector{String}()
for i = eachindex(processed_files)
    push!(processed_trunc, processed_files[i][length(ppath)+1:end-4])
end

error_log = []
for fid = eachindex(files)
    # for fid = 13:15

    results_file_exists = isfile("Results.csv")

    println("File $fid/$(length(files))")
    filename = String127(filenames_trunc[fid])

    if results_file_exists
        existing_data = DataFrame(CSV.File("Results.csv"))
        filename in existing_data.recording_name ? continue : 0
    end
    # if sum(occursin.(processed_trunc, filenames_trunc[fid])) > 0
    #     continue
    # end

    try

        time_taken = @elapsed begin
            global Data = readCSF(files[fid])
            # global Data = readCSF("/Users/jjc/CSF/Recordings/Inf_20180419114929_INF1.hdf5")
            icp = Data["ICP"][Data["infusion_start_frame"]:Data["infusion_end_frame"]]
            icp = icp[.~isnan.(icp)]

            # Priors and bounds
            global params_means = [10.45, 0.33, 7.5]
            global params_stddevs = [2.03, 0.08, 1.5]
            global lower_bounds = [0.01, 0.01, -10.0]
            global upper_bounds = [50.0, 1.0, Data["P_b"]]
            global Ib_max = 1.00
            global Ib_min = 0.01

            # Turing.jl settings
            sampler = NUTS()
            # sampler = MH()
            num_samples = 5_000
            num_chains = 4

            # Sampling
            # x0 = [0.33, 7.5, 10.45, 1.0] # starting point, but prevents convergence in some files so not used currently
            # chain = sample(curve_fitting(icp), sampler, MCMCThreads(), num_samples, num_chains, init_params=Iterators.repeated(x0))
            chain = sample(curve_fitting(icp), sampler, MCMCThreads(), num_samples, num_chains)

            # Populate results
            res_summary = DataFrame(summarize(chain))
            E_mean, P0_mean, Rout_mean = res_summary.mean[1:3]
            Ib_chain = (Data["P_b"] .- chain[:P0]) ./ chain[:Rout]
            Ib_mean = (Data["P_b"] - P0_mean) / Rout_mean
            df = DataFrame(chain)
            acceptance_rate = round(mean(df.acceptance_rate), digits=2)
            idxoptim = findfirst(df.lp .== maximum(df.lp))
            pms = df[idxoptim, [:E, :Rout, :P0]]
            Rout_MAP = pms.Rout
            E_MAP = pms.E
            P0_MAP = pms.P0
            Ib_MAP = (Data["P_b"] - P0_MAP) / Rout_MAP

            Rhats = DataFrame(rhat(chain))
            Rhat = maximum(Rhats.rhat) # take the worst-converged chain statistic

            NRMSE_Bayes_mean = round(calc_model_plot(Ib_mean, E_mean, P0_mean, P0_mean)[2], digits=3)
            NRMSE_Bayes_MAP = round(calc_model_plot(Ib_MAP, E_MAP, P0_MAP, P0_MAP)[2], digits=3)
            NRMSE_GD = round(calc_model_plot(Data["I_b"], Data["E"], Data["P_0"], Data["P_0"])[2], digits=3)

            Rout_std = std(chain[:Rout])
            E_std = std(chain[:E])
            P0_std = std(chain[:P0])
            Ib_std = std(Ib_chain)

            Rout_mean, E_mean, P0_mean, Ib_mean, Rout_MAP, E_MAP, P0_MAP, Ib_MAP, Rhat, Rout_std, E_std, P0_std, Ib_std = round.([Rout_mean, E_mean, P0_mean, Ib_mean, Rout_MAP, E_MAP, P0_MAP, Ib_MAP, Rhat, Rout_std, E_std, P0_std, Ib_std], digits=2)

            Rout_GD = Data["Rcsf"]
            E_GD = Data["E"]
            P0_GD = Data["P_0"]
            Ib_GD = Data["I_b"]
        end

        execution_time = round(time_taken, digits=1)

        push!(results,
            (filename,
                Rout_GD,
                E_GD,
                P0_GD,
                Ib_GD,
                Rout_mean,
                E_mean,
                P0_mean,
                Ib_mean,
                Rout_std,
                E_std,
                P0_std,
                Ib_std,
                Rout_MAP,
                E_MAP,
                P0_MAP,
                Ib_MAP,
                Rhat,
                NRMSE_Bayes_mean,
                NRMSE_Bayes_MAP,
                NRMSE_GD,
                acceptance_rate,
                execution_time
            ), promote=true)

        # Create and save plots
        plotmodel(Ib_mean, E_mean, P0_mean, P0_mean, zeros(3), zeros(3), "dark", "")

        title!(
            "Resistance to CSF outflow = $(round(Rout_mean,digits=2)) ± $(round(std(chain[:Rout]), digits=2)) [mmHg/mL/min]\n" *
            "Elasticity coefficient = $(round(E_mean,digits=2)) ± $(round(std(chain[:E]), digits=2)) [1/mL]\n" *
            "Reference pressure = $(round(P0_mean,digits=2)) ± $(round(std(chain[:P0]), digits=2)) [mmHg]\n" *
            "CSF production rate = $(round(Ib_mean,digits=2)) ± $(round(std(Ib_chain), digits=2)) [mL/min]\n" *
            "Error (Bayesian) = $NRMSE_Bayes_mean\n" *
            "Error (Gradient descent) = $NRMSE_GD\n",
            grid=true,
            size=(700, 500),
            dpi=300,
            margin=10mm,
            legend=:topleft
        )

        savefig("ModelPlots/$filename.png")

        autocorplot(chain, dpi=300)
        savefig("ACFPlots/$filename.png")

        plot(chain, dpi=300)
        savefig("ChainPlots/$filename.png")

        Plots.CURRENT_PLOT.nullableplot = nothing # Clear plots so they don't accumulate

        if results_file_exists
            existing_data = DataFrame(CSV.File("Results.csv"))
            check_exists = existing_data.recording_name .== filename
            res_exists = length(findall(check_exists)) > 0
            res_idx = findfirst(check_exists)
            if res_exists
                existing_data[res_idx, :] = results[1, :]
            else
                push!(existing_data, results[1, :], promote=true)
            end
            CSV.write("Results.csv", existing_data, append=false)
        else
            CSV.write("Results.csv", results, writeheader=true)
        end

    catch e
        # push!(error_log, "Error at file '$(filenames_trunc[fid])': $e")
        bt = catch_backtrace()
        stack_trace = sprint(Base.show_backtrace, bt)
        # push!(error_log, "Error at file $fid '$(filenames_trunc[fid])': $e\n$stack_trace")

        open("/Users/jjc/BayesCSF/ErrorLogs/ErroLogFile$fid.txt", "w") do file
            write(file, stack_trace)
        end
        continue
    end
    empty!(results)
end

