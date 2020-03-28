using Statistics;
using Plots; gr()
using Plots.PlotMeasures;
using Test;
using LaTeXStrings;
using StatsPlots;
using Distributions;

function Canonical_MonteCarlo(ρ::Type, T::Type, R_Cut::Type = 3.) where {Type <: Real}
    ##################################### CONFIGURATIONAL STEPS #############################
    println("\t\tCANONICAL MONTE CARLO")
    MC_Relaxation_Steps = 20_000;
    MC_Equilibrium_Steps = 250_000;
    MC_Measurement = 10;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;
    
    ##################################### VARIABLE INITIALIZATION ###########################
    L = 20.;
    V = L ^ 3.;
    N = convert(Int64, round(V*ρ))
    σ, λ = 1.0, 1.5;
    Beta, μ_Id, N_Bins = 1.0 / T, T * log(ρ), 150;
    Displacement, N_Displacement, N_Displacement_Accepted, N_Measurements = 0.1*ρ^(-1.0/3.0), 0, 0, 0;
    Energy_Array, μ_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Mean_Energy_Array, Mean_μ_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    STD_Energy_Array, STD_μ_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    g_r = zeros(Float64, N_Bins)
    ####################################### OUTPUT ROUTE, INITIAL POSITIONS AND INITIAL ENERGY ###########################
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/Density_$(round(ρ, digits = 2))/"
    mkpath(Output_Route)
    x, y, z = InitialPositions(N, L);
    Energy = Total_Energy_Calculation(L, x, y, z);
    ################################################# SIMULATION CYCLES ###########################################
    println("- Starting Simulation Cycles")
    for k = 1:MC_Steps
        @test all(Array(x) .<= L / 2.) && all(Array(x) .>= -L / 2.)
        @test all(Array(y) .<= L / 2.) && all(Array(y) .>= -L / 2.)
        @test all(Array(z) .<= L / 2.) && all(Array(z) .>= -L / 2.)
        ####################################### PRINTS SIMULATION PROGRESS TO SCREEN ########################################
        if k < MC_Relaxation_Steps && k % .01MC_Relaxation_Steps == 0
            println("$(convert(Int64, 100k / MC_Relaxation_Steps))% RELAXATION: [ρ* = $ρ, L* = $L, T* = $T]")
            println("N = $N Particles")
            println("U* / N = $(round(Energy / length(x), digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 2))%)\tRejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100 - 100N_Displacement_Accepted/N_Displacement, digits = 2))%)\n")
            N_Displacement, N_Displacement_Accepted = 0, 0;
        end

        if k > MC_Relaxation_Steps && k % .01MC_Equilibrium_Steps == 0
            println("$(convert(Int64, ceil(100(k - MC_Relaxation_Steps) / MC_Equilibrium_Steps)))% EQUILIBRIUM: [ρ* = $ρ, L* = $L, T* = $T, $(N_Measurements + 1) Measurements]")
            println("N = $N Particles")
            println("U* / N = $(round(Energy / length(x), digits = 6))")
            println("μ* = $(round(μ_Array[N_Measurements], digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 2))%)\tRejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100 - 100N_Displacement_Accepted/N_Displacement, digits = 2))%)\n")
            N_Displacement, N_Displacement_Accepted = 0, 0;
        end

        k == MC_Relaxation_Steps ? println("- FINISHED RELAXATION STEPS\n") : nothing
        @inbounds for i = 1:N
            N_Displacement += 1;
            Energy, N_Displacement_Accepted = Movement(i, L, Beta, Displacement, Energy, N_Displacement_Accepted, x, y, z)
        end
        ##################################################### MEASUREMENTS ##################################################
        if k % MC_Measurement == 0
            if k > MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                μ_Ex = TestParticleInsertion(L, Beta, x, y, z)
                μ_Array[N_Measurements] = μ_Id + μ_Ex;
                Mean_Energy_Array[N_Measurements] = mean(Energy_Array[1:N_Measurements]);
                Mean_μ_Array[N_Measurements] = mean(μ_Array[1:N_Measurements])
                if N_Measurements > 1
                    STD_Energy_Array[N_Measurements] = std(Energy_Array[1:N_Measurements]);
                    STD_μ_Array[N_Measurements] = std(μ_Array[1:N_Measurements]);
                end
                g_r += RadialDistributionFunction(N_Bins, L, length(x) / V, x, y, z)
            end
            1. * N_Displacement_Accepted / N_Displacement > 0.55 ? Displacement *= 1.05 : Displacement *= 0.95
            Displacement < 0.05 ? Displacement = 0.05 : nothing
            Displacement > L / 4. ? Displacement = L / 4. : nothing
        end  
    end
    ####################################################### END OF SIMULATION CYCLES #############################################
    ########################################################### SUMMARY FILE ##############################################
    println("< U* / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 6)) ± $(round(STD_Energy_Array[N_Measurements], digits = 6))")
    println("< μ* > = $(round(Mean_μ_Array[N_Measurements], digits = 6)) ± $(round(STD_μ_Array[N_Measurements], digits = 6))")
    Summary_File = open("$Output_Route/Summary.dat", "w+")
    println(Summary_File, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   INPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, "Density = $ρ\tN = $N\nL = $L\tV = $V\tT = $T\n$MC_Relaxation_Steps Relaxation Steps.\n$MC_Equilibrium_Steps Equilibrium Steps.\tMeasurements every $MC_Measurement steps.")
    println(Summary_File, "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   OUTPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, "< U* / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 6)) ± $(round(STD_Energy_Array[N_Measurements], digits = 6))")
    println(Summary_File, "< μ* > = $(round(Mean_μ_Array[N_Measurements], digits = 6)) ± $(round(STD_μ_Array[N_Measurements], digits = 6))")
    close(Summary_File)
    ################################################################ OUTPUT ##############################################################
    Energy_File = open("$Output_Route/Energy.dat", "w+");
    μ_File = open("$Output_Route/ChemicalPotential.dat", "w+");
    println(Energy_File, "Step\tEnergy\tMean_Energy\tStd_Energy")
    println(μ_File, "Step\tChemicalPotential\tMean_ChemicalPotential\tStd_ChemicalPotential")
    for i = 1:N_Measurements
        println(Energy_File, "$i\t$(round(Energy_Array[i], digits = 6))\t$(round(Mean_Energy_Array[i], digits = 6))\t$(round(STD_Energy_Array[i], digits = 6))")
        println(μ_File, "$i\t$(round(μ_Array[i], digits = 6))\t$(round(Mean_μ_Array[i], digits = 6))\t$(round(STD_μ_Array[i], digits = 6))")
    end
    close(Energy_File)
    close(μ_File)
    ##################################################### RADIAL DISTRIBUTION FUNCTION ########################################################
    g_r ./= N_Measurements;
    Delta, r = R_Cut / N_Bins, zeros(Float64, N_Bins);
    g_r_File = open("$Output_Route/RadialDistribution.dat", "w+");
    println(g_r_File, "r\tg_r\n")
    @inbounds for i = 1:N_Bins
        r[i] = round((i + 0.5)*Delta, digits = 6);
        println(g_r_File, "$(r[i])\t$(round(g_r[i], digits = 6))")
    end
    close(g_r_File)
    
    ###################################################### PLOTS ##############################################
    Radial_Distribution_Plot = plot(r, g_r, xlabel = L"r^*", ylabel = L"g(r^*)", xlim = (0, 3), xticks = 0:0.5:3, width = 3, guidefontsize = 20, tickfontsize = 18, bottom_margin = 7mm, left_margin = 10mm, right_margin = 3mm, legend = false, size = [1920, 1080], dpi = 300)
    hline!([1.0], color = :black, width = 2, linestyle = :dash)
    savefig(Radial_Distribution_Plot, "$Output_Route/RadialDistribution")

    Energy_Plot = plot(Energy_Array, xlabel = "N. Measurements", ylabel = L"U^* / N", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!(Energy_Plot, [Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    Energy_Histogram =  histogram(Energy_Array, normalize = true, xlabel = L"U^* / N", ylabel = "Normalized Frequency", legend = false, framestyle = :box, bins = 20, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(Mean_Energy_Array[N_Measurements], STD_Energy_Array[N_Measurements]), width = 3, linecolor = :black)

    Mean_Energy_Plot = plot(Mean_Energy_Array, ribbon = STD_Energy_Array, xlabel = "N. Measurements", ylabel = L"\langle U^* / N \rangle", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    Energy_Plots = plot(Energy_Plot, Energy_Histogram, Mean_Energy_Plot, layout = (@layout [a{0.3h} ; b c]))
    savefig(Energy_Plots, "$Output_Route/Energy_Plots")

    μ_Plot = plot(μ_Array, xlabel = "N. Measurements", ylabel = L"\mu^*", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!(μ_Plot, [Mean_μ_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    μ_Histogram =  histogram(μ_Array, normalize = true, xlabel = L"\mu^*", ylabel = " Normalized Frequency", legend = false, framestyle = :box, bins = 20, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([Mean_μ_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(Mean_μ_Array[N_Measurements], STD_μ_Array[N_Measurements]), width = 3, linecolor = :black)

    Mean_μ_Plot = plot(Mean_μ_Array, ribbon = STD_μ_Array, xlabel = "N. Measurements", ylabel = L"\langle \mu^* \rangle", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([Mean_μ_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

    μ_Plots = plot(μ_Plot, μ_Histogram, Mean_μ_Plot, layout = (@layout [a{0.3h} ; b c]))
    savefig(μ_Plots, "$Output_Route/ChemicalPotential_Plots")
end

function InitialPositions(N::Int64, L::Type) where {Type <: Real}
    x, y, z, i = zeros(Float64, N), zeros(Float64, N), zeros(Float64, N), 1;
    while i <= N
        x[i], y[i], z[i] = L * (rand() - 0.5), L * (rand() - 0.5), L * (rand() - 0.5);
        Overlap = false;
        for j = 1:i
            if i != j
                Delta_x, Delta_y, Delta_z = x[i] - x[j], y[i] - y[j], z[i] - z[j];
                Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
                Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
                Delta_z = PeriodicBoundaryConditions!(L, Delta_z);   
                r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
                if r2 < 1.
                    Overlap = true
                    break
                end
            end
        end
        Overlap == true ? nothing : i += 1
    end
    println("\n- Initial Configuration Established\n")
    return x, y, z
end

function Total_Energy_Calculation(L::Type, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, σ::Type = 1., λ::Type = 1.5, R_Cut::Type = 3.) where {Type <: Real}
    Energy = 0.;
    @inbounds for i = 1:length(x) - 1, j = i + 1:length(x)
        Delta_x, Delta_y, Delta_z = x[j] - x[i], y[j] - y[i], z[j] - z[i];
        Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
        Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
        Delta_z = PeriodicBoundaryConditions!(L, Delta_z);
        r = sqrt(Delta_x^2. + Delta_y^2. + Delta_z^2.);
        Energy += U(r);
    end
    return Energy
end

function Movement(i::Int64, L::Type, Beta::Float64, Displacement::Float64, Energy::Float64, N_Displacement_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    Energy_Old = Energy_Calculation(L, x[i], y[i], z[i], x, y, z);
    x_Old, y_Old, z_Old = x[i], y[i], z[i];
    x[i] += Displacement * (rand() - 0.5);
    y[i] += Displacement * (rand() - 0.5);
    z[i] += Displacement * (rand() - 0.5);
    x[i] = PeriodicBoundaryConditions!(L, x[i]);
    y[i] = PeriodicBoundaryConditions!(L, y[i]);
    z[i] = PeriodicBoundaryConditions!(L, z[i]);
    Energy_New = Energy_Calculation(L, x[i], y[i], z[i], x, y, z);
    Delta_E = Energy_New - Energy_Old;
    if rand() < exp(-Beta * Delta_E)
        N_Displacement_Accepted += 1;
        Energy += Delta_E;
    else
        x[i], y[i], z[i] = x_Old, y_Old, z_Old;
    end
    return Energy, N_Displacement_Accepted
end

function Energy_Calculation(L::Type, rx::Float64, ry::Float64, rz::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    Energy = 0.;
    @inbounds for i = 1:length(x)
        Delta_x, Delta_y, Delta_z = rx - x[i], ry - y[i], rz - z[i];
        Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
        Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
        Delta_z = PeriodicBoundaryConditions!(L, Delta_z);
        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
        r != 0. ? Energy += U(r) : nothing
    end
    return Energy
end

function PeriodicBoundaryConditions!(L::Type, x::Float64) where {Type <: Real}
    return x - L * round(x / L)
end

function U(r::Float64, σ::Type = 1., λ::Type = 1.5, e::Type = 1.) where {Type <: Real}
    r <= σ ? (return Inf) : r <= λ ? (return -e) : (return 0)
end

function U_LJ(r::Float64, σ::Type = 1., e::Type = 1.)
    return 4e * (r^(-12.) - r^(-6))
end

function TestParticleInsertion(L::Type, Beta::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, Insertions::Int64 = 5000, R_Cut::Type = 3.) where {Type <: Real}
    μ_Sum = 0.;
    for i = 1:Insertions
        Energy = 0;
        x_μ, y_μ, z_μ = L * (rand() - 0.5), L * (rand() - 0.5), L * (rand() - 0.5);
        for j = 1:length(x)
            Delta_x, Delta_y, Delta_z = x_μ - x[j], y_μ - y[j], z_μ - z[j];
            Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
            Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
            Delta_z = PeriodicBoundaryConditions!(L, Delta_z);
            r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
            Energy += U(r);
        end
        μ_Sum += exp(-Beta * Energy);
    end
    return -log(μ_Sum / Insertions) / Beta
end

function RadialDistributionFunction(N_Bins::Int64, L::Type, Density::Type, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Type = 3.) where {Type <: Real}
    Delta = R_Cut / N_Bins;
    g_r = zeros(Float64, N_Bins);
    @inbounds for i = 1:length(x) - 1, j = i + 1:length(x)
        Delta_x, Delta_y, Delta_z = x[i] - x[j], y[i] - y[j], z[i] - z[j];
        Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
        Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
        Delta_z = PeriodicBoundaryConditions!(L, Delta_z);
        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
        if r <= R_Cut
            l = convert(Int64, ceil(r / Delta));
            g_r[l] += 2;
        end
    end
    @inbounds for l = 1:N_Bins
        g_r[l] /= (length(x) * 4 * pi * (l * Delta)^2 * Delta * Density)
    end
    return g_r
end

Canonical_MonteCarlo(parse(Float64, ARGS[1]), parse(Float64, ARGS[2]))