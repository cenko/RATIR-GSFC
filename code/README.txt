To add a new instrument:

- Only edit "specific_instruments.py"

- Create instrument with inherited class "instrument"
    - initialize instrument with instrument name and
      number of cameras
    - Define the functions that were defined in the
      instrument class as abstract methods:
    
        has_cam_bias
        has_cam_dark
        change_header_keywords
        slice
        is_cam_split
        get_cam_sat
        get_cam_gain
        get_filter
        possible_filters
        get_centered_filter
        original_file_format
        change_file_names
     
       Information of what each of these functions
       should do (input/outputs) are commented in
       the instrument class
    - In the change_header_keywords must have or create
      the following keywords for future pipeline use:
      
        NAXIS1
        NAXIS2
        EXPTIME
        FILTER
        CRPIX1
        CRPIX2
        CTYPE1
        CTYPE2
        CD1_1
        CD2_1
        CD1_2
        CD2_2
        SATURATE
        GAIN
        TARGNAME
        PIXSCALE
        WAVELENG
        
        (PV*_* if there are distortion parameters)

- Add instrument with instrument name and class name to
  instrument_dict