"""Software parameters for pipeline"""


def get_trimmomatic_parameters(software_parameters):
    trimmomatic_params = (
        f'ILLUMINACLIP {software_parameters["trimmo-ILLUMINACLIP"]} -phred33'
        if software_parameters["trimmo-ILLUMINACLIP"] is not None
        else "-phred33"
    )
    trimmomatic_params += (
        f' HEADCROP:{software_parameters["trimmo-HEADCROP"]}'
        if software_parameters["trimmo-HEADCROP"] is not None
        else ""
    )
    trimmomatic_params += (
        f' CROP:{software_parameters["trimmo-CROP"]}'
        if software_parameters["trimmo-CROP"] is not None
        else ""
    )
    trimmomatic_params += (
        f' SLIDINGWINDOW:{software_parameters["trimmo-SLIDINGWINDOW1"]}:{software_parameters["trimmo-SLIDINGWINDOW2"]}'
        if software_parameters["trimmo-SLIDINGWINDOW2"] is not None
        else ""
    )
    trimmomatic_params += (
        f' LEADING:{software_parameters["trimmo-LEADING"]}'
        if software_parameters["trimmo-LEADING"] is not None
        else ""
    )
    trimmomatic_params += (
        f' TRAILING:{software_parameters["trimmo-TRAILING"]}'
        if software_parameters["trimmo-TRAILING"] is not None
        else ""
    )
    trimmomatic_params += (
        f' MINLEN:{software_parameters["trimmo-MINLEN"]}'
        if software_parameters["trimmo-MINLEN"] is not None
        else ""
    )
    trimmomatic_params += " TOPHRED33"
    return trimmomatic_params


def mask_regions_parameters(software_parameters):
    mask = (
        f'-s {software_parameters["MASK_SITES"]} '
        if software_parameters["MASK_SITES"] is not None
        else ""
    )
    mask += (
        f'-r {software_parameters["MASK_REGIONS"]} '
        if software_parameters["MASK_REGIONS"] is not None
        else ""
    )
    mask += (
        f'-b {software_parameters["MASK_F_BEGINNING"]} '
        if software_parameters["MASK_F_BEGINNING"] is not None
        else ""
    )
    mask += (
        f'-e {software_parameters["MASK_F_END"]} '
        if software_parameters["MASK_F_END"] is not None
        else ""
    )
    return mask


def get_snippy_parameters(software_parameters, primer):
    snippy_parameters = (
        f' --mapqual {software_parameters["mapqual"]}'
        if software_parameters["mapqual"]
        else ""
    )
    snippy_parameters += (
        f' --mincov {software_parameters["mincov"]}'
        if software_parameters["mincov"]
        else ""
    )
    snippy_parameters += (
        f' --minfrac {software_parameters["minfrac"] }'
        if software_parameters["minfrac"]
        else ""
    )
    snippy_parameters += (
        f' --primer ../../../{primer}'
        if primer
        else ""
    )
    return snippy_parameters


def get_nanofilt_parameters(software_parameters):
    nanofilt_parameters = (
        f' -q {software_parameters["QUALITY"]}'
        if software_parameters["QUALITY"] is not None
        else ""
    )
    nanofilt_parameters += (
        f' -l {software_parameters["LENGTH"]}'
        if software_parameters["LENGTH"] is not None
        else ""
    )
    nanofilt_parameters += (
        f' --headcrop {software_parameters["HEADCROP"]}'
        if software_parameters["HEADCROP"] is not None
        else ""
    )
    nanofilt_parameters += (
        f' --tailcrop {software_parameters["TAILCROP"]}'
        if software_parameters["TAILCROP"] is not None
        else ""
    )
    return nanofilt_parameters
