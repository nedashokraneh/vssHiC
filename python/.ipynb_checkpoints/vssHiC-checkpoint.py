import cooler
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from random import sample

def make_rep_df(cools, chromosomes):
    reps_df = pd.DataFrame()
    offsets = {}
    for i, chrom1 in enumerate(chromosomes):
        for chrom2 in chromosomes[i:]:
            mats = {}
            for rep in [1, 2]:
                mats[rep] = cools[rep].matrix(balance = False).fetch(chrom1, chrom2)
            offsets[chrom1, chrom2] = reps_df.shape[0]
            if chrom1 == chrom2:
                ut_idx = np.triu_indices(mats[1].shape[0], 1)
                reps_df = pd.concat([reps_df, pd.DataFrame({"rep1": mats[1][ut_idx], "rep2": mats[2][ut_idx]})])
            else:
                reps_df = pd.concat([reps_df, pd.DataFrame({"rep1": mats[1].flatten(), "rep2": mats[2].flatten()})])
    return reps_df, offsets

def get_mean_sd(reps_df, aux_type = "rep1", group_type = "bin", zero_bin = True, bin_size = 100):
    if aux_type == "rep1":
        base_type = "rep2"
    else:
        base_type = "rep1"
    nz_rows = ((reps_df["rep1"] + reps_df["rep2"]) != 0)
    order = np.argsort(reps_df[nz_rows][aux_type].values)
    ordered_nz_df = reps_df[nz_rows].iloc[order,:]
    if zero_bin:
        zero_counts = ordered_nz_df[ordered_nz_df[aux_type] == 0][base_type].values
        zero_mean = np.mean(zero_counts)
        zero_sd = np.sqrt(np.var(zero_counts))
        ordered_nz_df = ordered_nz_df[ordered_nz_df[aux_type] != 0]
    ordered_nz_df["index"] = (np.arange(ordered_nz_df.shape[0]) / bin_size).astype(int)
    mean_var_df = ordered_nz_df.groupby("index").agg({base_type: ["mean", "var"]})
    mean_var_df.columns = mean_var_df.columns.droplevel()
    mean_var_df["var"] = np.sqrt(mean_var_df["var"])
    mean_var_df.rename(columns = {"var": "stdev"}, inplace = True)
    if zero_bin:
        mean_var_df = pd.concat([pd.DataFrame({"mean": [zero_mean], "stdev": [zero_sd]}), mean_var_df])
    return mean_var_df

def fit_mean_sd(mean_sd_df, log_scale = True):
    mean_sd_df.drop_duplicates(inplace = True)
    if log_scale:
        r_x = robjects.FloatVector(np.log(mean_sd_df["mean"].values + 1))
        r_y = robjects.FloatVector(np.log(mean_sd_df["stdev"].values + 1))
    else:
        r_x = robjects.FloatVector(mean_sd_df["mean"].values)
        r_y = robjects.FloatVector(mean_sd_df["stdev"].values)
    r_smooth_spline = robjects.r['smooth.spline'] #extract R function
    kwargs = {"x": r_x, "y": r_y, "lambda":  0.1}
    spline_xy = r_smooth_spline(**kwargs)
    #spline_xy = r_smooth_spline(x = r_x, y = r_y, lambda = r_lambda)
    pred_sd = np.array(robjects.r['predict'](spline_xy, r_x).rx2('y'))
    if log_scale:
        mean_sd_df["pred_sd"] = np.exp(pred_sd) - 1
    else:
        mean_sd_df["pred_sd"] = pred_sd
    return mean_sd_df, spline_xy

def vss_transform(reps_df, spline, log_scale = True, log_converge = False):
    max_ct = np.max(reps_df)
    xg = np.sinh(np.arange(np.arcsinh(0), np.arcsinh(max_ct), (np.arcsinh(max_ct) - np.arcsinh(0)) / 1000))
    if log_scale:
        r_xg = robjects.FloatVector(np.log(xg + 1))
    else:
        r_xg = robjects.FloatVector(xg)
    xg_sd = np.array(robjects.r['predict'](spline, r_xg).rx2('y'))
    if log_scale:
        xg_sd = np.exp(xg_sd) - 1
    inv_xg_sd = 1 / xg_sd
    splf_x = robjects.FloatVector(np.arcsinh((xg[1:] + xg[:-1]) / 2))
    cum_inv_xg_sd = np.cumsum((xg[1:] - xg[:-1]) * ((inv_xg_sd[1:] + inv_xg_sd[:-1]) / 2))
    splf_y = robjects.FloatVector(cum_inv_xg_sd)
    r_splinefun = robjects.r['splinefun'] #extract R function
    trans_splinefun = r_splinefun(x = splf_x, y = splf_y)
    vss_reps_df = pd.DataFrame({"rep1": np.array(trans_splinefun(robjects.FloatVector(np.arcsinh(reps_df["rep1"].values)))),
                                "rep2": np.array(trans_splinefun(robjects.FloatVector(np.arcsinh(reps_df["rep2"].values))))})
    if log_converge:
        h1 = np.quantile(np.maximum(reps_df["rep1"].values, reps_df["rep2"].values), 0.95)
        h2 = np.quantile(np.maximum(reps_df["rep1"].values, reps_df["rep2"].values), 0.9999)
        eta = (np.log2(h2) - np.log2(h1)) / (np.array(trans_splinefun(robjects.FloatVector([np.arcsinh(h2)]))) - 
                                             np.array(trans_splinefun(robjects.FloatVector([np.arcsinh(h1)]))))
        xi = np.log2(h1) - eta * np.array(trans_splinefun(robjects.FloatVector([np.arcsinh(h1)])))
        vss_reps_df = vss_reps_df * eta[0] + xi[0]
    return vss_reps_df

def vss(reps_df, mean_sd_df, log_scale = True, log_converge = False):
    if log_scale:
        r_x = robjects.FloatVector(np.log(mean_sd_df["mean"].values + 1))
        r_y = robjects.FloatVector(np.log(mean_sd_df["stdev"].values + 1))
    else:
        r_x = robjects.FloatVector(mean_sd_df["mean"].values)
        r_y = robjects.FloatVector(mean_sd_df["stdev"].values)
    r_smooth_spline = robjects.r['smooth.spline'] #extract R function
    spline_xy = r_smooth_spline(x = r_x, y = r_y)
    pred_sd = np.array(robjects.r['predict'](spline_xy, r_x).rx2('y'))
    if log_scale:
        mean_sd_df["pred_sd"] = np.exp(pred_sd) - 1
    else:
        mean_sd_df["pred_sd"] = pred_sd
    max_ct = np.max(reps_df)
    xg = np.sinh(np.arange(np.arcsinh(0), np.arcsinh(max_ct), (np.arcsinh(max_ct) - np.arcsinh(0)) / 1000))
    if log_scale:
        r_xg = robjects.FloatVector(np.log(xg + 1))
    else:
        r_xg = robjects.FloatVector(xg)
    xg_sd = np.array(robjects.r['predict'](spline_xy, r_xg).rx2('y'))
    if log_scale:
        xg_sd = np.exp(xg_sd) - 1
    inv_xg_sd = 1 / xg_sd
    splf_x = robjects.FloatVector(np.arcsinh((xg[1:] + xg[:-1]) / 2))
    cum_inv_xg_sd = np.cumsum((xg[1:] - xg[:-1]) * ((inv_xg_sd[1:] + inv_xg_sd[:-1]) / 2))
    splf_y = robjects.FloatVector(cum_inv_xg_sd)
    r_splinefun = robjects.r['splinefun'] #extract R function
    trans_splinefun = r_splinefun(x = splf_x, y = splf_y)
    vss_reps_df = pd.DataFrame({"rep1": np.array(trans_splinefun(robjects.FloatVector(np.arcsinh(reps_df["rep1"].values)))),
                                "rep2": np.array(trans_splinefun(robjects.FloatVector(np.arcsinh(reps_df["rep2"].values))))})
    if log_converge:
        h1 = np.quantile(np.maximum(reps_df["rep1"].values, reps_df["rep2"].values), 0.95)
        h2 = np.quantile(np.maximum(reps_df["rep1"].values, reps_df["rep2"].values), 0.9999)
        eta = (np.log2(h2) - np.log2(h1)) / (np.array(trans_splinefun(robjects.FloatVector([np.arcsinh(h2)]))) - 
                                             np.array(trans_splinefun(robjects.FloatVector([np.arcsinh(h1)]))))
        xi = np.log2(h1) - eta * np.array(trans_splinefun(robjects.FloatVector([np.arcsinh(h1)])))
        vss_reps_df = vss_reps_df * eta[0] + xi[0]
    return vss_reps_df, mean_sd_df
    
    return mean_sd_df, spline_xy