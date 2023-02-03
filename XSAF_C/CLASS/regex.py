float_regex = r"(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
bool_regex = r"(?:[FT]{1})"
nc_spaces_regex = r"(?:\s*)"
equal_reqex = fr"{nc_spaces_regex}=){nc_spaces_regex}"
redshift_regex = r"(transfer_interp_matterpower = F)"
redshift_rm_regex = (
    fr"\ntransfer_redshift\(\d*\){nc_spaces_regex}={nc_spaces_regex}{float_regex}"
)

transfer_file_rm_regex = (
    fr"\ntransfer_filename\(\d*\){nc_spaces_regex}={nc_spaces_regex}transfer_out\d*.dat"
)
regexs = {
    "do_nonlinear": fr"(\ndo_nonlinear{equal_reqex} {float_regex}",
    "w": fr"(\nw{equal_reqex} {float_regex}",
    "wa": fr"(\nwa{equal_reqex} {float_regex}",
    "ombh2": fr"(\nombh2{equal_reqex} {float_regex}",
    "omch2": fr"(\nomch2{equal_reqex} {float_regex}",
    "hubble": fr"(\nhubble{equal_reqex} {float_regex}",
    "scalar_amp": fr"(\nscalar_amp\(1\){equal_reqex} {float_regex}",
    "scalar_spectral_index": (
        fr"(\nscalar_spectral_index\(1\){equal_reqex} {float_regex}"
    ),
    "transfer_num_redshifts": fr"(\ntransfer_num_redshifts{equal_reqex} {float_regex}",
    "output_root": fr"(\noutput_root{equal_reqex} (?:.*)\n",
}
