using Plots
using Statistics
using StatsBase
import Base.length
# histograms

"""Define an histogram"""
struct Histo1d
	edges::Vector{Float64}
	weights::Vector{Float64}
	centers::Vector{Float64}
end

"""
Take an object of type Histogram and return a Histo1d
"""
get_histo1d(h::Histogram) = Histo1d(edges(h), h.weights, centers(h))

"""
Returns a Histo1d struct

""" 
function h1d(x::Vector{T}, nbins::Integer,
             xmin::T=typemin(T), xmax::T=typemax(T), norm=false) where T

    h = hist1d(x, nbins, xmin, xmax, norm)
    Histo1d(edges(h), h.weights, centers(h))
end

h1d(x::Vector{Float64}, bins::Vector{Float64}, norm=false) = get_histo1d(hist1d(x, bins, norm))


"""plot a h1d histogram"""
function plot_h1dx(h; xs="edges value", ys ="weights values",
                   yl=(1.0,10^4), xl=(0.0,100.0),
                   ylog=false, xlog=false,
                   ylim=false, xlim=false,
                   markersize=3, fwl=false,
                   label="", legend=false)

    xg = h.centers
    yg = h.weights

    marker = (:circle, markersize, 0.6, :black)

    xaxis =(xs)
    if xlim && xlog
    xaxis = (xs, xl, :log)
    elseif xlim
    xaxis = (xs, xl)
    elseif xlog
    xaxis = (xs, :log)
    end

    yaxis =(ys)
    if ylim && ylog
    yaxis = (ys, yl, :log)
    elseif ylim
    yaxis = (ys, yl)
    elseif ylog
    yaxis = (ys, :log)
    end


    p = scatter(xg, yg, yerr=sqrt.(abs.(yg)), 
    xaxis =xaxis, yaxis = yaxis, 
    markersize=markersize, label=label, legend=legend)
    if fwl
    p = plot!(p, xg, yg, 
        xaxis =xaxis, yaxis = yaxis, fmt=:png,
        linewidth=1, label=label, legend=legend)
    end

    p
    end


function plot_h1dx(h1, h2; xs="edges value", ys ="weights values",
        yl=(1.0,10^4), xl=(0.0,100.0),
        ylog=false, xlog=false,
        ylim=false, xlim=false,
        markersize=3, fwl=false,
        label1="", label2="", legend=false)

    xg1 = h1.centers
    yg1 = h1.weights

    xg2 = h2.centers
    yg2 = h2.weights

    marker = (:circle, markersize, 0.6, :black)

    xaxis =(xs)
    if xlim && xlog
    xaxis = (xs, xl, :log)
    elseif xlim
    xaxis = (xs, xl)
    elseif xlog
    xaxis = (xs, :log)
    end

    yaxis =(ys)
    if ylim && ylog
    yaxis = (ys, yl, :log)
    elseif ylim
    yaxis = (ys, yl)
    elseif ylog
    yaxis = (ys, :log)
    end


    p = scatter(xg1, yg1, yerr=sqrt.(abs.(yg1)), 
    xaxis =xaxis, yaxis = yaxis, 
    markersize=markersize, label=label1, legend=legend)
    p = scatter!(p, xg2, yg2, yerr=sqrt.(abs.(yg2)), 
    xaxis =xaxis, yaxis = yaxis, 
    markersize=markersize, label=label2, legend=legend)
    if fwl
    p = plot!(p, xg1, yg1, 
    xaxis =xaxis, yaxis = yaxis, fmt=:png,label=label1,
    linewidth=1, legend=legend)
    p = plot!(p, xg2, yg2, label=label2,
    xaxis =xaxis, yaxis = yaxis, fmt=:png,
    linewidth=1, legend=legend)
    end

    p
    end


"""
    digitize(x, bins)

    Return the indices of the bins to which each value in input array belongs.
"""
function digitize(x::Vector{<:Real}, bins::LinRange{<:Real})
    return searchsortedlast.(Ref(bins), x)
end
# Notice that in this case the bins are expected to be sorted
function digitize(x::Vector{<:Real}, bins::Vector{<:Real})
    return searchsortedlast.(Ref(bins), x)
end


"""
    hist1d(x::Vector{T}, nbins::Integer, xmin::T=typemin(T), xmax::T=typemax(T))
    hist1d(x::Vector{T}, xs::String, nbins::Integer,
                    xmin::T=typemin(T), xmax::T=typemax(T); datap = true, fwl=true)
    hist1d(h::Histogram, xs::String; datap = true, markersize=3, fwl=false)

return a 1d histogram and its corresponding graphics (plots)

"""
function hist1d(x::Vector{<:Real}, bins::Vector{<:Real}, norm=false)
    h = fit(Histogram, x, bins)
    if norm
        h = StatsBase.normalize(h, mode=:density)
    end
    return h
end


function hist1d(x::Vector{T}, nbins::Integer,
                xmin::T=typemin(T), xmax::T=typemax(T), norm=false) where T
    
    dx = (xmax - xmin) /nbins

    edges = [xmin + i * dx for i in 0:nbins]
    h = fit(Histogram, x, edges)
    if norm
        h = StatsBase.normalize(h, mode=:density)
    end
    h

    #xx   = in_range(x, xmin, xmax)
    #bins = collect(LinRange(xmin, xmax, nbins + 1))
    #h    = hist1d(xx, bins, norm)
    #return h
end


"""
    edges(h::Histogram)

edges of the histogram
"""
function edges(h::Histogram, index::Int64=1)
    collect(h.edges[index])
end


"""
    centers(h::Histogram)

return the centres of the histogram bins.
"""
centers(h::Histogram) = centers(edges(h))


"""
    centers(edges::Vector{<:Real})

Calculate the bin centers from a vector of bin edges
"""
function centers(edges::Vector{<:Real})
    edges[1:end-1] + .-(@view(edges[2:end]), @view(edges[1:end-1])) / 2
end


"""
    hist_weights

return histogram weights of a set of data.
"""
function hist_weights(edges::Vector{<:Real})
  function get_weights(y::SubArray{<:Real})
    histo = fit(Histogram, y, edges)
    return histo.weights
  end
  return get_weights
end

"""
    hist1d(x::Vector{T}, xs::String, nbins::Integer,
           xmin::T=typemin(T), xmax::T=typemin(T);
           norm::Bool=false, datap::Bool=true,
           markersize::Int64=3, fwl::Bool=false,
           label::String="", legend::Bool=false) where T <: Real
return 1D histogram and plot axes for data x and number of bins nbins.
"""
function hist1d(x::Vector{T}, xs::String, nbins::Integer,
                xmin::T=typemin(T), xmax::T=typemax(T);
                norm::Bool=false, datap::Bool=true,
                markersize::Int64=3, fwl::Bool=false,
                label::String="", legend::Bool=false) where T <: Real

    h = hist1d(x, nbins, xmin, xmax, norm)
    return h, hist1d(h, xs, datap=datap, markersize=markersize, fwl=fwl,
                     label=label, legend=legend)
end


"""
    hist1d(h::Histogram, xs::String;
           datap::Bool=true, markersize::Int64=3, fwl::Bool=false,
           label::String="", legend::Bool=false)
Plot 1D histogram from Histogram object and options.
"""
function hist1d(h::Histogram, xs::String;
                datap::Bool=true, markersize::Int64=3, fwl::Bool=false,
                label::String="", legend::Bool=false)

    if datap
        xg = centers(h)
        yg = eltype(xg).(h.weights)
        p  = scatter(xg, yg, yerr=sqrt.(yg),
                     markersize=markersize, label=label, legend=legend)
        if fwl
            p = plot!(p, xg, yg, yerr=sqrt.(yg), fmt=:png,
                      linewidth=1, label=label, legend=legend)
        end
    else
        p = plot(h, xlabel=xs, fmt=:png, yl="frequency")
    end
    xlabel!(xs)
    ylabel!("frequency")

    return p
end


"""
    plot2hist1d(h1::Histogram, h2::Histogram, xs::String;
                markersize::Int64=2, norm::Bool=false)
Plot two 1D histograms on the same axes.
"""
function plot2hist1d(h1::Histogram, h2::Histogram, xs::String;
                     markersize::Int64=2, norm::Bool=false)

    if norm
        h1 = StatsBase.normalize(h1, mode=:density)
        h2 = StatsBase.normalize(h2, mode=:density)
    end

    xg1 = centers(h1)
    yg1 = eltype(xg1).(h1.weights)
    xg2 = centers(h2)
    yg2 = eltype(xg2).(h2.weights)

    p1 = scatter(     xg1, yg1, yerr=sqrt.(yg1), fmt=:png, markersize=markersize)
    p1 = scatter!(p1, xg2, yg2, yerr=sqrt.(yg2), fmt=:png, markersize=markersize)
    xlabel!(xs)
    ylabel!("frequency")

    return p1
end


"""
    hist2d(x::Vector{T}, y::Vector{T}, nbins::Integer,
           xl::String, yl::String,
           xmin::T=typemin(T), xmax::T=typemax(T),
           ymin::T=typemin(T), ymax::T=typemax(T)) where T <: Real

ggiven data x, y, labels xl and yl and limits for each axis.
"""
function hist2d(x::Vector{T},y::Vector{T}, nbins::Integer,
                xl::String, yl::String,
                xmin::T=typemin(T), xmax::T=typemax(T),
                ymin::T=typemin(T), ymax::T=typemax(T); title="") where T <: Real
    function select_data(xval::T, yval::T)
        range_bound(xmin, xmax, OpenBound)(xval) .&
        range_bound(ymin, ymax, OpenBound)(yval)
    end
    mask = select_data.(x, y)
    data = (y[mask], x[mask])
    h    = fit(Histogram, data, nbins=nbins)
    hm   = heatmap(xe, ye, h.weights)
    xlabel!(xl)
    ylabel!(yl)
    if title != ""
        title!(title)
    end
    return h, hm
end


"""
    fmean_std

returns the weighted mean and std of an array given a minimum
bin weight for consideration.
"""
function fmean_std(x::Vector{<:Real}, min_prop::Float64=0.1)
    mask_func = weights -> weights .> min_prop * maximum(weights)
    function filt_mean(w::SubArray{<:Real})
        mask = mask_func(w)
        return [mean_and_std(x[mask], FrequencyWeights(w[mask]), corrected=true)]
    end
    return filt_mean
end


"""
    p1df(x, y, nbins; ybin_width, ymin, ymax, min_proportion)

return a profile DataFrame. This is a DF in which the variable y is
histogrammed as a function of the average of variable x in each bin.
The keyword arguments allow for a filtering in the y variable by
range and min_proportion of the maximum in the bin.
TODO: Protect against fine binning that results in zeros.
"""
function p1df(x::Vector{<:Real}, y::Vector{<:Real}, nbins::Integer;
              ybin_width::Real=0.1, ymin::Real=minimum(y), ymax::Real=maximum(y),
              min_proportion::Real=0.0)
    df           = DataFrame(:y => y)
    x_upper      = nextfloat(maximum(x))
    bin_edges    = LinRange(minimum(x), x_upper, nbins + 1)
    df[!, "bin"] = digitize(x, bin_edges)
    bin_centers  = centers(collect(bin_edges))
    bin_width    = .-(@view(bin_edges[2:end]), @view(bin_edges[1:end-1]))
    ## Bin in y so filtering can be done on bin content
    y_upper      = nextfloat(ymax)
    y_bins       = collect(ymin:ybin_width:y_upper)
    y_centers    = centers(y_bins)
    y_weights    = combine(groupby(df, :bin),
                           :y => hist_weights(y_bins) => :yweights,
                           ungroup = false)
    bin_stats    = fmean_std(y_centers, min_proportion)
    ndf = combine(y_weights,
                  :yweights => bin_stats => [:y_mean, :y_std])
    ndf = ndf[:, [:y_mean, :y_std]]
    ndf[!, :x_mean] = bin_centers
    ndf[!, :x_std]  = bin_width / 2

    p1 = plot(ndf.x_mean,ndf.y_mean, yerror=ndf.y_std, fmt = :png,
              shape = :circle, color = :black, legend=false)
    return ndf, p1
end
