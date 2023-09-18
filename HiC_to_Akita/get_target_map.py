from cooltools.lib.numutils import observed_over_expected, adaptive_coarsegrain
from cooltools.lib.numutils import interpolate_bad_singletons, set_diag, interp_nan
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
import numpy as np

def get_target(chrom,start,end,genome_hic_cool):
    pool_width = 2048
    crop_bp = 65536
    diagonal_offset = 2
    kernel_stddev = 1
    clip = 2
    seq_len_nt = end - start
    seq_len_pool = seq_len_nt // pool_width
    crop_start = crop_bp // pool_width
    crop_end = seq_len_pool - crop_start
    seq_len_crop = seq_len_pool - 2*crop_start
    triu_tup = np.triu_indices(seq_len_crop, diagonal_offset)
    seq_len_nodiag = seq_len_crop - diagonal_offset
    seq_len_hic = seq_len_nodiag*(seq_len_nodiag + 1) // 2
    ## based on Aktia paper, this line should be "kernel = Gaussian2DKernel(x_stddev=kernel_stddev,x_size=5)";
    ## however, on basenji github, it is "kernel = Gaussian2DKernel(x_stddev=kernel_stddev)", 
    ## I tested both and found there was no significant difference, but you can retest it
    kernel = Gaussian2DKernel(x_stddev=kernel_stddev)
    
    mseq_str = '%s:%d-%d' % (chrom, start, end)
    seq_hic_raw = genome_hic_cool.matrix(balance=True).fetch(mseq_str)
    seq_hic_nan = np.isnan(seq_hic_raw)
    num_filtered_bins = np.sum(np.sum(seq_hic_nan,axis=0) == len(seq_hic_nan))
    if num_filtered_bins > (.5*len(seq_hic_nan)):
        print("WARNING: %s >50% bins filtered, check:  %s. " % (genome_hic_file, mseq_str))
    
    clipval = np.nanmedian(np.diag(seq_hic_raw,diagonal_offset))
    for i in range(-diagonal_offset+1,diagonal_offset):
        set_diag(seq_hic_raw, clipval, i)
    seq_hic_raw = np.clip(seq_hic_raw, 0, clipval)
    seq_hic_raw[seq_hic_nan] = np.nan
    
    seq_hic_smoothed = adaptive_coarsegrain(
                            seq_hic_raw,
                            genome_hic_cool.matrix(balance=False).fetch(mseq_str),
                            cutoff=2, max_levels=8)
    seq_hic_nan = np.isnan(seq_hic_smoothed)
    
    seq_hic_obsexp = observed_over_expected(seq_hic_smoothed, ~seq_hic_nan)[0]
    seq_hic_obsexp = np.log(seq_hic_obsexp)
    seq_hic_obsexp = np.clip(seq_hic_obsexp, -clip, clip)
    seq_hic_obsexp = interp_nan(seq_hic_obsexp)
    for i in range(-diagonal_offset+1, diagonal_offset): set_diag(seq_hic_obsexp, 0,i)
    seq_hic = convolve(seq_hic_obsexp, kernel)
    
    seq_hic = seq_hic[crop_start:crop_end,:]
    seq_hic = seq_hic[:,crop_start:crop_end]
    seq_hic = seq_hic[triu_tup].astype('float16')
    
    return seq_hic
