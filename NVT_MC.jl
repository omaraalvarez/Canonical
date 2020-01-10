using Statistics;
using Plots;
using Test;

function MonteCarlo(ρ::Float64, L::Float64, T::Float64, R_Cut::Float64 = 3.0)
    println("NVT MONTE CARLO")
    """ CONFIGURATIONAL STEPS """
    MC_Relaxation_Steps = 100_000;
    MC_Equilibrium_Steps = 250_000;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;
    MC_Measurement = 10;
    """ VARIABLE INITIALIZATION """
    V = L ^ 3.;
    N = convert(Int64, ceil(V*ρ))
    σ_p, λ_p = 0.5, 1.5;
    Beta = 1.0/T;
    Displacement, N_Displacement, N_Displacement_Accepted = 0.1*ρ^(-1.0/3.0), 0, 0;
    Energy_Sum, μ_Sum, μ_Ex_Sum, N_Measurements = 0., 0., 0., 0;
    Energy_Array, μ_Array, μ_Ex_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Average_Energy_Array, Average_μ_Array, Average_μ_Ex_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    STD_Energy_Array, STD_μ_Array, STD_μ_Ex_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    N_Bins = 150;
    g_r = zeros(Float64, N_Bins)
    """ OUTPUT FILES """
    Output_Route = pwd() * "/Output/Density_$(round(ρ, digits = 2))/T_$(round(T, digits = 2))"
    mkpath(Output_Route)
    Average_Energy_File = open("$Output_Route/Average_Energy.dat", "w");
    println(Average_Energy_File, "#\t< E / N >")
    Energy_File = open("$Output_Route/Energy.dat", "w");
    println(Energy_File, "#\tE / N ")
    Average_Chemical_File = open("$Output_Route/AverageChemicalPotential.dat", "w");
    println(Average_Chemical_File, "#\t< mu >")
    Chemical_File = open("$Output_Route/ChemicalPotential.dat", "w");
    println(Chemical_File, "#\tmu ")
    Average_ExcessChemical_File = open("$Output_Route/AverageExcessChemicalPotential.dat", "w");
    println(Average_Chemical_File, "#\t< mu_ex >")
    ExcessChemical_File = open("$Output_Route/ExcessChemicalPotential.dat", "w");
    println(Chemical_File, "#\tmu_ex ")
    """ INITIAL POSITIONS   """
    x, y, z = InitialPositions(N, L);
    μ, μ_Ex = 0., 0.;
    Energy = Total_Energy_Calculation(N, L, R_Cut, x, y, z);
    """ SIMULATION CYCLES """
    for k = 1:MC_Steps

        @test all(Array(x) .<= L / 2.) && all(Array(x) .>= -L / 2.)
        @test all(Array(y) .<= L / 2.) && all(Array(y) .>= -L / 2.)
        @test all(Array(z) .<= L / 2.) && all(Array(z) .>= -L / 2.)

        """     PRINTS PROGRESS TO SCREEN   """
        if k < MC_Relaxation_Steps && k % .01MC_Relaxation_Steps == 0
            println("$(convert(Int64, 100k / MC_Relaxation_Steps))% Relaxation. ($N Particles)")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 4))%)")
            println("   Rejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100 - 100N_Displacement_Accepted/N_Displacement, digits = 4))%)")
            println("")
        end

        if k > MC_Relaxation_Steps && k % .01MC_Equilibrium_Steps == 0
            println("$(convert(Int64, floor(100(k - MC_Relaxation_Steps) / MC_Equilibrium_Steps)))% Equilibrium ($N_Measurements Measurements). ($N Particles)")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("μ_ex = $(round(μ_Ex, digits = 6))")
            println("μ = $(round(μ, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 4))%)")
            println("   Rejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100 - 100N_Displacement_Accepted/N_Displacement, digits = 4))%)")
            println("")
        end

        k == MC_Relaxation_Steps ? println("~~~    STARTING MEASUREMENT STEPS    ~~~") : nothing
        @inbounds for i = 1:N
            N_Displacement += 1;
            Energy, N_Displacement_Accepted = Movement(L, Beta, Displacement, Energy, N_Displacement_Accepted, R_Cut, x, y, z, i)
        end

        if k % MC_Measurement == 0
            if k > MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                println(Energy_File, "$N_Measurements\t$(round(Energy_Array[N_Measurements], digits = 6))")
                Energy_Sum += Energy / length(x);
                Average_Energy_Array[N_Measurements] = Energy_Sum / N_Measurements;
                if N_Measurements > 1
                    STD_Energy_Array[N_Measurements] = std(Energy_Array[1:N_Measurements]);
                end
                println(Average_Energy_File, "$N_Measurements\t$(round(Average_Energy_Array[N_Measurements], digits = 6))\t$(round(STD_Energy_Array[N_Measurements], digits = 6))")

                μ_Ex = WidomInsertion(L, Beta, R_Cut, x, y, z)
                μ_Ex = - T * log(μ_Ex)
                μ_Ex_Array[N_Measurements] = μ_Ex;
                println(ExcessChemical_File, "$N_Measurements\t$(round(μ_Ex_Array[N_Measurements], digits = 6))")
                μ_Ex_Sum += μ_Ex;
                Average_μ_Ex_Array[N_Measurements] = μ_Ex_Sum / N_Measurements;
                if N_Measurements > 1
                    STD_μ_Ex_Array[N_Measurements] = std(μ_Ex_Array[1:N_Measurements]);
                end
                println(Average_ExcessChemical_File, "$N_Measurements\t$(round(Average_μ_Ex_Array[N_Measurements], digits = 6))\t$(round(STD_μ_Ex_Array[N_Measurements], digits = 6))")

                μ = T * log(ρ) + μ_Ex;
                μ_Array[N_Measurements] = μ;
                println(Chemical_File, "$N_Measurements\t$(round(μ_Array[N_Measurements], digits = 6))")
                μ_Sum += μ;
                Average_μ_Array[N_Measurements] = μ_Sum / N_Measurements;
                if N_Measurements > 1
                    STD_μ_Array[N_Measurements] = std(μ_Array[1:N_Measurements]);
                end
                println(Average_Chemical_File, "$N_Measurements\t$(round(Average_μ_Array[N_Measurements], digits = 6))\t$(round(STD_μ_Array[N_Measurements], digits = 6))")

                g_r += Distribution(N_Bins, L, length(x) / V, x, y, z, R_Cut)
            end
            1. * N_Displacement_Accepted / N_Displacement > 0.55 ? Displacement *= 1.05 : Displacement *= 0.95
            Displacement < 0.05 ? Displacement = 0.05 : nothing
            Displacement > L / 4. ? Displacement = L / 4. : nothing
            N_Displacement, N_Displacement_Accepted = 0, 0;
        end  

    end
    close(Average_Energy_File)
    close(Energy_File)
    close(Average_Chemical_File)
    close(Chemical_File)
    close(Average_ExcessChemical_File)
    close(ExcessChemical_File)

    g_r /= N_Measurements;
    Delta = R_Cut / N_Bins;
    r = zeros(Float64, N_Bins - 1)
    g_r_File = open("$Output_Route/Radial_Distribution.dat", "w");
    println(Average_Energy_File, "#r\t#Normalized Density\n")
    @inbounds for i = 1:N_Bins - 1
        r[i] = round((i + 0.5)*Delta, digits = 6)
        println(g_r_File, "$(r[i])\t$(round(g_r[i], digits = 6))")
    end
    g_r = g_r[1 : N_Bins - 1]
    close(g_r_File)
    Radial_Distribution_Plot = plot(r, g_r, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Distance [r]", ylabel = "Normalized Density", width = 3, size = [1200, 800])
    hline!([1.0], color = :black, width = 2, linestyle = :dash)
    savefig(Radial_Distribution_Plot, "$Output_Route/RadialDistribution")

    Summary_File = open("$Output_Route/Summary.dat", "w")

    println("< E / N > = $(round(mean(Energy_Array[1:end - 1]), digits = 6)) ± $(round(std(Energy_Array[1:end - 1]), digits = 6))")
    println(Summary_File, "< E / N > = $(round(mean(Energy_Array[1:end - 1]), digits = 6)) ± $(round(std(Energy_Array[1:end - 1]), digits = 6))")
    Energy_Plot = plot(Energy_Array[1:end - 1], guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Measurements", ylabel = "Energy [Unitless]", width = 2, size = [1200, 800])
    hline!([mean(Energy_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Energy_Plot, "$Output_Route/Energy")
    Energy_Histogram = histogram(Energy_Array[convert(Int64, floor(MC_Relaxation_Steps/MC_Measurement)):end - 1], guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, bins = 20, legend = false, xlabel = "Energy [Unitless]", ylabel = "Frequency", size = [1200, 800])
    vline!([mean(Energy_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Energy_Histogram, "$Output_Route/Energy_Histogram")
    Average_Energy_Plot = plot(Average_Energy_Array[1:end - 1], ribbon = STD_Energy_Array, fillalpha = 0.2, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Measurements", ylabel = "< Energy > [Unitless]", width = 3, size = [1200, 800])
    hline!([mean(Energy_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Average_Energy_Plot, "$Output_Route/Average_Energy")

    println("< mu_Ex > = $(round(mean(μ_Ex_Array[1:end - 1]), digits = 6)) ± $(round(std(μ_Ex_Array[1:end - 1]), digits = 6))")
    println(Summary_File, "< mu_Ex > = $(round(mean(μ_Ex_Array[1:end - 1]), digits = 6)) ± $(round(std(μ_Ex_Array[1:end - 1]), digits = 6))")
    ExcessChemical_Plot = plot(μ_Ex_Array[1:end - 1], guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Measurements", ylabel = "Excess Chemical Potential [Unitless]", width = 2, size = [1200, 800])
    hline!([mean(μ_Ex_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(ExcessChemical_Plot, "$Output_Route/ExcessChemicalPotential")
    ExcessChemical_Histogram = histogram(μ_Ex_Array[convert(Int64, floor(MC_Relaxation_Steps/MC_Measurement)):end - 1], guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, bins = 20, legend = false, xlabel = "Excess Chemical Potential [Unitless]", ylabel = "Frequency", size = [1200, 800])
    vline!([mean(μ_Ex_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(ExcessChemical_Histogram, "$Output_Route/ExcessChemicalPotential_Histogram")
    Average_ExcessChemical_Plot = plot(Average_μ_Ex_Array[1:end - 1], ribbon = STD_μ_Ex_Array, fillalpha = 0.2, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Measurements", ylabel = "< Excess Chemical Potential > [Unitless]", width = 3, size = [1200, 800])
    hline!([mean(μ_Ex_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Average_ExcessChemical_Plot, "$Output_Route/Average_ExcessChemicalPotential")

    println("< mu > = $(round(mean(μ_Array[1:end - 1]), digits = 6)) ± $(round(std(μ_Array[1:end - 1]), digits = 6))")
    println(Summary_File, "< mu > = $(round(mean(μ_Array[1:end - 1]), digits = 6)) ± $(round(std(μ_Array[1:end - 1]), digits = 6))")
    Chemical_Plot = plot(μ_Array[1:end - 1], guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Measurements", ylabel = "Total Chemical Potential [Unitless]", width = 2, size = [1200, 800])
    hline!([mean(μ_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Chemical_Plot, "$Output_Route/ChemicalPotential")
    Chemical_Histogram = histogram(μ_Array[convert(Int64, floor(MC_Relaxation_Steps/MC_Measurement)):end - 1], guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, bins = 20, legend = false, xlabel = "Chemical Potential [Unitless]", ylabel = "Frequency", size = [1200, 800])
    vline!([mean(μ_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Chemical_Histogram, "$Output_Route/ChemicalPotential_Histogram")
    Average_Chemical_Plot = plot(Average_μ_Array[1:end - 1], ribbon = STD_μ_Array, fillalpha = 0.2, guidefontsize = 14, tickfontsize = 10, widen = true, dpi = 300, legend = false, xlabel = "Measurements", ylabel = "< Chemical Potential > [Unitless]", width = 3, size = [1200, 800])
    hline!([mean(μ_Array[1:end - 1])], color = :black, width = 2, linestyle = :dash)
    savefig(Average_Chemical_Plot, "$Output_Route/Average_ChemicalPotential")

    close(Summary_File)
    return mean(μ_Array[1:end - 1]), std(μ_Array[1:end - 1])
end

function InitialPositions(N::Int64, L::Float64)
    x, y, z = zeros(Float64, N), zeros(Float64, N), zeros(Float64, N);
    i = 1
    while i <= N
        x[i], y[i], z[i] = (L - 0.5) * (rand() - 0.5), (L - 0.5) * (rand() - 0.5), (L - 0.5) * (rand() - 0.5);
        overlap = false
        for j = 1:i
            Delta_x = x[i] - x[j];
            Delta_x = PeriodicBoundaryConditions(L, Delta_x);
            Delta_y = y[i] - y[j];
            Delta_y = PeriodicBoundaryConditions(L, Delta_y);
            Delta_z = z[i] - z[j];
            Delta_z = PeriodicBoundaryConditions(L, Delta_z);
            r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
            if r2 < 1. && r2 != 0
                overlap = true
                break
            end
        end
        overlap == true ? nothing : i += 1
    end
    return x, y, z
end

function Movement(L::Float64, Beta::Float64, Displacement::Float64, Energy::Float64, N_Displacement_Accepted::Int64, R_Cut::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, j::Int64)
    Energy_Old = Energy_Calculation(L, R_Cut, x[j], y[j], z[j], x, y, z)
    x_Old, y_Old, z_Old = x[j], y[j], z[j];
    x[j] += Displacement * (rand() - 0.5);
    x[j] = PeriodicBoundaryConditions(L, x[j]);
    y[j] += Displacement * (rand() - 0.5);
    y[j] = PeriodicBoundaryConditions(L, y[j]);
    z[j] += Displacement * (rand() - 0.5);
    z[j] = PeriodicBoundaryConditions(L, z[j]);
    Energy_New = Energy_Calculation(L, R_Cut, x[j], y[j], z[j], x, y, z);
    Delta_E = Energy_New - Energy_Old;
    if rand() < exp(-Beta * Delta_E)
        N_Displacement_Accepted += 1;
        Energy += Delta_E;
    else
        x[j] = x_Old;
        y[j] = y_Old;
        z[j] = z_Old;
    end
    return Energy, N_Displacement_Accepted
end

function Energy_Calculation(L::Float64, R_Cut::Float64, rx::Float64, ry::Float64, rz::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    Energy = 0.;
    @inbounds for i = 1:length(x)
        Delta_x = rx - x[i];
        Delta_x = PeriodicBoundaryConditions(L, Delta_x);
        Delta_y = ry - y[i];
        Delta_y = PeriodicBoundaryConditions(L, Delta_y);
        Delta_z = rz - z[i];
        Delta_z = PeriodicBoundaryConditions(L, Delta_z);
        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
        if r2 != 0.
            if r2 < R_Cut^2
                Energy += u_SquareWell(r2, 0.5, 1.5);
            end
        end
    end
    return Energy
end

function Total_Energy_Calculation(N::Int64, L::Float64, R_Cut::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    Energy = 0.;
    @inbounds for i = 1:N - 1, j = i + 1:N
        Delta_x = x[j] - x[i];
        Delta_x = PeriodicBoundaryConditions(L, Delta_x);
        Delta_y = y[j] - y[i];
        Delta_y = PeriodicBoundaryConditions(L, Delta_y);
        Delta_z = z[j] - z[i];
        Delta_z = PeriodicBoundaryConditions(L, Delta_z);
        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
        if r2 < R_Cut^2.
            Energy += u_SquareWell(r2, 0.5, 1.5);
        end
    end
    return Energy
end

function PeriodicBoundaryConditions(L::Float64, x::Float64)
    if x < - L / 2.
        x += L;
    elseif x > L / 2.
        x -= L;
    end
    return x
end

function u_SquareWell(r2::Float64, σ::Float64, λ::Float64, e::Float64 = 1.)
    if r2 <= (2σ)^2
        return Inf
    elseif r2 <= λ^2
        return -e
    else
        return 0
    end
end

function Distribution(N_Bins::Int64, L::Float64, Density::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, R_Cut::Float64)
    Delta = L / (2N_Bins);
    Delta = R_Cut / N_Bins;
    g_r = zeros(Float64, N_Bins);
    @inbounds for i = 1:length(x) - 1
        @inbounds for j = i + 1:length(x)
            Delta_x = x[i] - x[j];
            Delta_x = PeriodicBoundaryConditions(L, Delta_x);
            Delta_y = y[i] - y[j];
            Delta_y = PeriodicBoundaryConditions(L, Delta_y);
            Delta_z = z[i] - z[j];
            Delta_z = PeriodicBoundaryConditions(L, Delta_z);
            r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
            if r2 < R_Cut^2
                l = convert(Int64, floor(sqrt(r2) / Delta));
                g_r[l] += 2;
            end
        end
    end
    @inbounds for l = 1:N_Bins
        g_r[l] /= (length(x) * 4 * pi * (l * Delta)^2 * Delta * Density)
    end
    return g_r
end

function WidomInsertion(L::Float64, Beta::Float64, R_Cut::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, Insertions::Int64 = 5000)
    μ_Sum = 0.;
    for i = 1:Insertions
        x_μ = L * (rand() - 0.5);
        y_μ = L * (rand() - 0.5);
        z_μ = L * (rand() - 0.5);
        Energy = 0;
        for j = 1:length(x)
            Delta_x = x_μ - x[j];
            Delta_x = PeriodicBoundaryConditions(L, Delta_x);
            Delta_y = y_μ - y[j];
            Delta_y = PeriodicBoundaryConditions(L, Delta_y);
            Delta_z = z_μ - z[j];
            Delta_z = PeriodicBoundaryConditions(L, Delta_z);
            r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
            if r2 < R_Cut^2
                Energy += u_SquareWell(r2, 0.5, 1.5)
            end
        end
        μ_Sum += exp(-Beta * Energy)
    end
    return μ_Sum / Insertions
end

function Cycle()
    T_Array = [3.];
    μ_Mean = zeros(Float64, length(T_Array));
    μ_Std = zeros(Float64, length(T_Array));
    i = 1;
    #Route = pwd();
    #Chemical_File = open("$Route/CriticalDensity_0.437.dat", "w+")
    #println(Chemical_File, "T mu mu_error")
    for T in T_Array
        μ_Mean[i], μ_Std[i] = MonteCarlo(0.6, 10., T)
        #println(Chemical_File, "$T $(μ_Mean[i]) $(μ_Std[i])")
        i += 1;
    end
    #close(Chemical_File)
    #Chemical_Plot = plot(T_Array, μ_Mean, yerror = μ_Std, title = "Critical Density = 0.3", legend = false, xlabel = "Temperature", ylabel= "< Chemical Potential > [Unitless]", width = 3, size = [1200, 800]) 
    #savefig(Chemical_Plot, "$Route/CriticalDensity_0.3")
end

Cycle()