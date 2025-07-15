import os
import glob
from landlab import load_params


def load_params_txt(ext=".txt"):
    """
    Finds a .txt file in the current working directory 
    and loads it using the load_params function from landlab.

    Args:
        ext (str): extension for the param file desired (string), e.g. 'txt' (default)

    Returns:
        params (dict): Parameters loaded from the .txt file.
    
    Raises:
        FileNotFoundError: If no .txt files are found.
        TypeError: If the loaded parameters are not a dictionary.
    """
    if not isinstance(ext, str):
        raise TypeError("The 'ext' argument must be a string (e.g., 'txt').")
    
    ext = ext if ext.startswith(".") else f".{ext}"
    print('ext', ext)

        
    cwd = os.getcwd(); ## current working directory
    print('cwd', cwd)
    paramfile = glob.glob(os.path.join(cwd, f"*{ext}"))[0]; ## find the params file (of whatever extension, e.g. txt) in the current folder
    print('pf', paramfile)
    
    ## No file found?? --> Raise an error 
    if not paramfile:
        raise FileNotFoundError("No {} parameter file found in the current directory.".format(ext))

    params = load_params(paramfile); ## load the params file (using landlab load_params)
    
    ## File not loaded as desired format --> Raise an error
    if not isinstance(params, dict):
        raise TypeError("Expected 'params' to be a dictionary.")

    return params





