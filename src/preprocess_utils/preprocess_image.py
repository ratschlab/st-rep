import numpy as np
import pyvips

format_to_dtype = {
    'uchar': np.uint8,
    'char': np.int8,
    'ushort': np.uint16,
    'short': np.int16,
    'uint': np.uint32,
    'int': np.int32,
    'float': np.float32,
    'double': np.float64,
    'complex': np.complex64,
    'dpcomplex': np.complex128,
}


def get_low_res_image(image_path, downsample_factor):
    image = pyvips.Image.new_from_file(image_path, access='sequential')
    image_low_res = image.resize(1 / downsample_factor)
    image_low_res_arr = np.ndarray(buffer=image_low_res.write_to_memory(),
                                   dtype=format_to_dtype[image_low_res.format],
                                   shape=[image_low_res.height, image_low_res.width, image_low_res.bands])
    return image_low_res_arr
